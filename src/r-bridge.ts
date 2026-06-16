/**
 * r-bridge.ts — Manages a persistent R child process running bridge.R.
 *
 * Protocol: newline-delimited JSON on stdin (requests) / stdout (responses).
 * Each request gets a UUID; responses are matched back by id.
 * R objects never cross the bridge — only string ref keys.
 */

import { spawn, type ChildProcess } from "node:child_process";
import { randomUUID } from "node:crypto";
import { resolve, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { existsSync } from "node:fs";

const __dirname = dirname(fileURLToPath(import.meta.url));

interface PendingRequest {
  resolve: (value: unknown) => void;
  reject: (reason: Error) => void;
  timer: ReturnType<typeof setTimeout>;
}

interface BridgeResponse {
  id: string;
  result?: unknown;
  error?: string;
}

export interface RBridgeOptions {
  /** Path to bridge.R (default: src/bridge.R relative to this file) */
  bridgePath?: string;
  /** Rscript binary name or path */
  rscriptBin?: string;
  /** Per-request timeout in ms (default: 300_000 = 5 min, training can be slow) */
  timeoutMs?: number;
  /** Extra env vars passed to the R process */
  env?: Record<string, string>;
}

export class RBridge {
  private proc: ChildProcess | null = null;
  private pending = new Map<string, PendingRequest>();
  private buffer = "";
  private readyPromise: Promise<Record<string, unknown>>;
  private readyResolve!: (v: Record<string, unknown>) => void;
  private readyReject!: (e: Error) => void;
  private opts: Required<RBridgeOptions>;

  constructor(options: RBridgeOptions = {}) {
    this.opts = {
      bridgePath: options.bridgePath ?? resolve(__dirname, "..", "src", "bridge.R"),
      rscriptBin: options.rscriptBin ?? "Rscript",
      timeoutMs: options.timeoutMs ?? 600_000,
      env: options.env ?? {},
    };

    this.readyPromise = new Promise((res, rej) => {
      this.readyResolve = res;
      this.readyReject = rej;
    });
  }

  /** Spawn the R process and wait for the ready signal. */
  async start(): Promise<Record<string, unknown>> {
    if (this.proc) throw new Error("R bridge already started");

    if (!existsSync(this.opts.bridgePath)) {
      throw new Error(
        `bridge.R not found at: ${this.opts.bridgePath}\n` +
        `  Expected location relative to dist/: ../src/bridge.R`
      );
    }

    const env = { ...process.env, ...this.opts.env };

    this.proc = spawn(this.opts.rscriptBin, [this.opts.bridgePath], {
      stdio: ["pipe", "pipe", "pipe"],
      env,
    });

    this.proc.stdout!.on("data", (chunk: Buffer) => this.onData(chunk));

    this.proc.stderr!.on("data", (chunk: Buffer) => {
      // Forward R logs to MCP server stderr so they show in terminal
      process.stderr.write(chunk);
    });

    this.proc.on("exit", (code, signal) => {
      const msg = `R process exited (code=${code}, signal=${signal})`;
      process.stderr.write(`[r-bridge] ${msg}\n`);
      // Reject all pending requests
      for (const [id, req] of this.pending) {
        req.reject(new Error(msg));
        clearTimeout(req.timer);
        this.pending.delete(id);
      }
      this.readyReject(new Error(msg));
      this.proc = null;
    });

    this.proc.on("error", (err) => {
      this.readyReject(
        new Error(`Failed to spawn Rscript: ${err.message}. Is R installed?`)
      );
    });

    // Wait for the __ready__ signal from bridge.R
    const readyTimeout = setTimeout(() => {
      this.readyReject(
        new Error(
          "R bridge did not become ready within 600s (Julia init may have hung)"
        )
      );
    }, 600_000);

    const info = await this.readyPromise;
    clearTimeout(readyTimeout);
    return info;
  }

  /** Send a command to R and await the result. */
  async call(method: string, params: Record<string, unknown> = {}): Promise<unknown> {
    if (!this.proc?.stdin?.writable) {
      throw new Error("R bridge is not running — call start() first");
    }

    const id = randomUUID();
    const request = JSON.stringify({ id, method, params }) + "\n";

    return new Promise<unknown>((resolve, reject) => {
      const timer = setTimeout(() => {
        this.pending.delete(id);
        reject(new Error(`R bridge timeout after ${this.opts.timeoutMs}ms for method '${method}'`));
      }, this.opts.timeoutMs);

      this.pending.set(id, { resolve, reject, timer });
      this.proc!.stdin!.write(request);
    });
  }

  /** Gracefully shut down the R process. */
  async stop(): Promise<void> {
    if (!this.proc) return;
    this.proc.stdin!.end();
    // Give R 5s to clean up, then kill
    await new Promise<void>((resolve) => {
      const forceKill = setTimeout(() => {
        this.proc?.kill("SIGKILL");
        resolve();
      }, 5_000);
      this.proc!.on("exit", () => {
        clearTimeout(forceKill);
        resolve();
      });
    });
  }

  // ── private ────────────────────────────────────────────────────────────────

  private onData(chunk: Buffer): void {
    this.buffer += chunk.toString();

    // Process complete lines
    let newlineIdx: number;
    while ((newlineIdx = this.buffer.indexOf("\n")) !== -1) {
      const line = this.buffer.slice(0, newlineIdx).trim();
      this.buffer = this.buffer.slice(newlineIdx + 1);
      if (!line) continue;

      let parsed: BridgeResponse;
      try {
        parsed = JSON.parse(line) as BridgeResponse;
      } catch {
        // Not JSON — R startup noise that leaked to stdout
        process.stderr.write(`[r-bridge] non-JSON stdout: ${line}\n`);
        continue;
      }

      // Handle the ready signal separately
      if (parsed.id === "__ready__") {
        this.readyResolve(parsed.result as Record<string, unknown>);
        continue;
      }

      const req = this.pending.get(parsed.id);
      if (!req) {
        process.stderr.write(`[r-bridge] orphan response id=${parsed.id}\n`);
        continue;
      }

      this.pending.delete(parsed.id);
      clearTimeout(req.timer);

      if (parsed.error) {
        req.reject(new Error(parsed.error));
      } else {
        req.resolve(parsed.result);
      }
    }
  }
}
