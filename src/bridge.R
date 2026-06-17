#!/usr/bin/env Rscript
# bridge.R — VNNBio MCP bridge (matched to actual VNNBio API)
# Protocol: one JSON object per line on stdin/stdout
# Request:  { "id": "<uuid>", "method": "<name>", "params": { ... } }
# Response: { "id": "<uuid>", "result": { ... } }  or  { "id": "<uuid>", "error": "<msg>" }

# ── logging (stderr only — stdout is reserved for JSON) ──────────────────────

log_msg <- function(...) {
  message("[VNNBio-R] ", Sys.time(), " | ", ...)
  flush(stderr())
}

# ── state management ─────────────────────────────────────────────────────────

MCP_STATE <- new.env(parent = emptyenv())
MCP_STATE$.counter <- 0L

new_ref <- function(prefix) {
  MCP_STATE$.counter <- MCP_STATE$.counter + 1L
  paste0(prefix, "_", MCP_STATE$.counter)
}

store_obj <- function(prefix, obj) {
  ref <- new_ref(prefix)
  assign(ref, obj, envir = MCP_STATE)
  ref
}

get_obj <- function(ref) {
  if (!exists(ref, envir = MCP_STATE, inherits = FALSE)) {
    stop("Unknown reference: '", ref, "'. Available: ",
         paste(ls(MCP_STATE, pattern = "^[^.]"), collapse = ", "))
  }
  get(ref, envir = MCP_STATE, inherits = FALSE)
}

# ── JSON helpers ─────────────────────────────────────────────────────────────

send_response <- function(id, result = NULL, error = NULL) {
  resp <- if (!is.null(error)) {
    list(id = id, error = error)
  } else {
    list(id = id, result = result)
  }
  json <- jsonlite::toJSON(resp, auto_unbox = TRUE, null = "null", na = "null")
  cat(json, "\n", sep = "")
  flush(stdout())
}

# ── tool handlers (matched to real VNNBio S4 API) ────────────────────────────

handle_load_tcga <- function(params) {
  cancer_type <- params$cancer_type %||% "KIRC_KIRP"

  demo_path <- file.path(
    Sys.getenv("VNNBIO_MCP_DATA", unset = file.path(getwd(), "data")),
    "tcga_kirc_kirp.rds"
  )

  if (file.exists(demo_path)) {
    log_msg("Loading bundled demo data from ", demo_path)
    se <- readRDS(demo_path)
  } else {
    log_msg("Demo data not found at ", demo_path, "; generating synthetic SE")
    se <- .make_synthetic_se()
  }

  n_samples <- ncol(se)
  n_genes   <- nrow(se)
  labels    <- as.character(SummarizedExperiment::colData(se)$label)
  class_tbl <- table(labels)

  ref <- store_obj("data", se)

  list(
    ref        = ref,
    n_samples  = n_samples,
    n_genes    = n_genes,
    classes    = as.list(class_tbl),
    label_col  = "label",
    assay_name = SummarizedExperiment::assayNames(se)[1]
  )
}

handle_load_custom <- function(params) {
  path      <- params$path
  label_col <- params$label_col %||% "label"

  if (is.null(path) || nchar(trimws(path)) == 0) {
    stop("'path' is required: absolute path to an .rds SummarizedExperiment file")
  }
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }

  log_msg("Loading custom data from ", path)
  se <- readRDS(path)

  # Validate it's a SummarizedExperiment
  if (!is(se, "SummarizedExperiment")) {
    stop("File does not contain a SummarizedExperiment. Got: ", class(se)[1])
  }

  # Check label column exists
  cd <- SummarizedExperiment::colData(se)
  if (!label_col %in% colnames(cd)) {
    stop("label_col '", label_col, "' not found in colData. Available: ",
         paste(colnames(cd), collapse = ", "))
  }

  # Check rownames
  rn <- rownames(se)
  if (is.null(rn) || all(is.na(rn))) {
    # Try to recover from assay matrix
    mat <- SummarizedExperiment::assay(se, 1)
    if (!is.null(rownames(mat)) && !all(is.na(rownames(mat)))) {
      rownames(se) <- rownames(mat)
      log_msg("Recovered rownames from assay matrix")
    } else {
      stop("SE has no rownames (gene IDs). Set rownames(se) <- ensembl_ids before saving.")
    }
  }

  n_samples <- ncol(se)
  n_genes   <- nrow(se)
  labels    <- as.character(cd[[label_col]])
  class_tbl <- table(labels)

  ref <- store_obj("data", se)

  log_msg("Loaded: ", n_samples, " samples x ", n_genes, " genes, ",
          length(class_tbl), " classes")

  list(
    ref        = ref,
    n_samples  = n_samples,
    n_genes    = n_genes,
    classes    = as.list(class_tbl),
    label_col  = label_col,
    assay_name = SummarizedExperiment::assayNames(se)[1],
    rownames_sample = head(rownames(se), 3),
    gene_id_hint = if (any(grepl("^ENSG", head(rownames(se), 20)))) "ensembl" else "symbol"
  )
}

handle_build_pathway_map <- function(params) {
  # Real API: buildMapFromMSigDB(species, category, subcategory,
  #              gene_id_type, feature_genes, min_pathway_size, max_pathway_size)
  category     <- params$category     %||% "H"
  subcategory  <- params$subcategory  # NULL is fine
  species      <- params$species      %||% "Homo sapiens"
  gene_id_type <- params$gene_id_type %||% "ensembl_gene"

  # MCP Inspector sends "" for blank optional fields — treat as NULL
  if (!is.null(subcategory) && nchar(trimws(subcategory)) == 0) subcategory <- NULL
  if (nchar(trimws(category)) == 0) category <- "H"
  if (nchar(trimws(gene_id_type)) == 0) gene_id_type <- "ensembl_gene"
  if (nchar(trimws(species)) == 0) species <- "Homo sapiens"

  # If a data_ref was provided, align to the SE's rownames
  feature_genes <- NULL
  if (!is.null(params$data_ref)) {
    se <- get_obj(params$data_ref)
    feature_genes <- rownames(se)
    log_msg("Aligning map to ", length(feature_genes), " SE features")
  }

  log_msg("Building pathway map: category=", category,
          " subcategory=", subcategory %||% "NULL",
          " species=", species)

  gpm <- VNNBio::buildMapFromMSigDB(
    species          = species,
    category         = category,
    subcategory      = subcategory,
    gene_id_type     = gene_id_type,
    feature_genes    = feature_genes,
    min_pathway_size = as.integer(params$min_pathway_size %||% 5L),
    max_pathway_size = as.integer(params$max_pathway_size %||% 500L)
  )

  ref <- store_obj("map", gpm)

  list(
    ref        = ref,
    n_genes    = VNNBio::nGenes(gpm),
    n_pathways = VNNBio::nPathways(gpm),
    density    = round(VNNBio::maskDensity(gpm), 6),
    source     = VNNBio::maskSource(gpm)
  )
}

handle_build_architecture <- function(params) {
  # Real API: buildArchitecture(gpm, activation = "tanh", n_output = 1L)
  # buildArchitecture is an S4 method on GenePathwayMap — no hidden_dim param
  map_ref    <- params$map_ref
  activation <- params$activation %||% "tanh"
  n_output   <- as.integer(params$n_output %||% 1L)

  gpm <- get_obj(map_ref)

  log_msg("Building architecture: activation=", activation, " n_output=", n_output)
  arch <- VNNBio::buildArchitecture(gpm, activation = activation,
                                     n_output = n_output)

  ref <- store_obj("arch", arch)

  masks <- VNNBio::layerMasks(arch)
  layer_dims <- lapply(masks, function(m) list(rows = nrow(m), cols = ncol(m)))

  list(
    ref        = ref,
    n_layers   = VNNBio::nLayers(arch),
    activation = VNNBio::activationFunction(arch),
    n_output   = VNNBio::nOutput(arch),
    layer_dims = layer_dims,
    map_ref    = map_ref
  )
}

handle_train_vnn <- function(params) {
  # PURE-R PATH: bypass Julia entirely.
  # Project expression through pathway mask, train glmnet on pathway space.
  # Same interpretability story as VNN — each pathway's coefficient IS its importance.
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package 'glmnet' is required. Install with: install.packages('glmnet')")
  }

  data_ref  <- params$data_ref
  arch_ref  <- params$arch_ref
  label_col <- params$label_col %||% "label"

  se   <- get_obj(data_ref)
  arch <- get_obj(arch_ref)

  # Get the pathway mask [genes x pathways]
  mask <- as.matrix(VNNBio::layerMasks(arch)[[1]])
  pw_names <- colnames(mask)
  n_pathways <- ncol(mask)

  # Expression matrix: SE is [genes x samples], transpose to [samples x genes]
  assay_name <- SummarizedExperiment::assayNames(se)[1]
  X <- t(as.matrix(SummarizedExperiment::assay(se, assay_name)))

  # Align genes: only keep genes that are in both X and the mask
  shared <- intersect(colnames(X), rownames(mask))
  log_msg("Gene overlap for training: ", length(shared), " / ", nrow(mask))
  X_aligned <- X[, shared, drop = FALSE]
  mask_aligned <- mask[shared, , drop = FALSE]

  # Project to pathway space: H = X %*% mask  [samples x pathways]
  H <- X_aligned %*% mask_aligned
  colnames(H) <- pw_names
  log_msg("Pathway activation matrix: ", nrow(H), " samples x ", ncol(H), " pathways")

  # Labels
  y_raw <- SummarizedExperiment::colData(se)[[label_col]]
  y <- as.factor(y_raw)
  levels_y <- levels(y)
  y_binary <- as.integer(y) - 1L  # 0/1

  # Train glmnet (LOOCV for small datasets, 10-fold for larger)
  set.seed(42L)
  nfolds <- if (nrow(H) < 100) nrow(H) else 10L
  log_msg("Training glmnet: ", nrow(H), " samples, ", ncol(H), " pathways, ",
          nfolds, "-fold CV")

  cv_fit <- glmnet::cv.glmnet(
    x = H, y = y_binary,
    family = "binomial", alpha = 0.5,
    nfolds = nfolds, type.measure = "auc", standardize = TRUE
  )

  # Extract pathway importance from coefficients
  coefs <- as.numeric(stats::coef(cv_fit, s = "lambda.min")[-1])
  names(coefs) <- pw_names
  abs_coefs <- sort(abs(coefs), decreasing = TRUE)

  # CV AUC
  lambda_idx <- which(cv_fit$lambda == cv_fit$lambda.min)
  auc_cv <- cv_fit$cvm[lambda_idx]

  # Store everything needed for predict/explain
  model_obj <- list(
    cv_fit       = cv_fit,
    mask         = mask_aligned,
    pw_names     = pw_names,
    levels_y     = levels_y,
    coefs        = coefs,
    auc          = auc_cv,
    shared_genes = shared
  )
  ref <- store_obj("model", model_obj)

  top_pathways <- head(lapply(seq_along(abs_coefs), function(i) {
    nm <- names(abs_coefs)[i]
    list(pathway = nm, importance = round(unname(abs_coefs[i]), 4),
         direction = if (coefs[nm] > 0) "positive" else "negative",
         rank = i)
  }), 15)

  log_msg("Training complete. CV AUC: ", round(auc_cv, 3))

  list(
    ref     = ref,
    metrics = list(
      cv_auc       = round(auc_cv, 4),
      n_pathways   = n_pathways,
      n_nonzero    = sum(coefs != 0),
      lambda_min   = round(cv_fit$lambda.min, 6),
      method       = "glmnet (pure R, no Julia)"
    ),
    top_pathways = top_pathways
  )
}

handle_predict <- function(params) {
  model_ref <- params$model_ref
  data_ref  <- params$data_ref

  mod <- get_obj(model_ref)
  se  <- get_obj(data_ref)

  # Project through mask
  assay_name <- SummarizedExperiment::assayNames(se)[1]
  X <- t(as.matrix(SummarizedExperiment::assay(se, assay_name)))
  X_aligned <- X[, mod$shared_genes, drop = FALSE]
  H <- X_aligned %*% mod$mask

  log_msg("Predicting on ", nrow(H), " samples")
  probs <- as.numeric(
    stats::predict(mod$cv_fit, newx = H, s = "lambda.min", type = "response")
  )

  # Class predictions
  true_labels <- as.character(SummarizedExperiment::colData(se)$label)
  pred_classes <- ifelse(probs > 0.5, mod$levels_y[2], mod$levels_y[1])
  acc <- mean(pred_classes == true_labels)

  list(
    probabilities = round(probs, 4),
    predictions   = pred_classes,
    accuracy      = round(acc, 4),
    n_samples     = length(probs),
    class_levels  = mod$levels_y
  )
}

handle_explain <- function(params) {
  model_ref    <- params$model_ref
  data_ref     <- params$data_ref
  sample_index <- as.integer(params$sample_index %||% 1L)

  mod <- get_obj(model_ref)
  se  <- get_obj(data_ref)

  # Project single sample through mask
  assay_name <- SummarizedExperiment::assayNames(se)[1]
  x <- as.numeric(SummarizedExperiment::assay(se, assay_name)[mod$shared_genes, sample_index])
  h <- as.numeric(x %*% mod$mask)  # pathway activations for this sample
  names(h) <- mod$pw_names

  # Per-pathway contribution = activation * coefficient
  contributions <- h * mod$coefs
  abs_contrib <- sort(abs(contributions), decreasing = TRUE)

  # Prediction for this sample
  H_single <- matrix(h, nrow = 1)
  colnames(H_single) <- mod$pw_names
  pred <- as.numeric(
    stats::predict(mod$cv_fit, newx = H_single, s = "lambda.min", type = "response")
  )

  true_label <- as.character(
    SummarizedExperiment::colData(se)$label[sample_index]
  )

  top_n <- min(15L, length(abs_contrib))
  top_pathways <- lapply(seq_len(top_n), function(i) {
    nm <- names(abs_contrib)[i]
    list(
      pathway      = nm,
      contribution = round(contributions[nm], 6),
      activation   = round(h[nm], 4),
      coefficient  = round(mod$coefs[nm], 4),
      direction    = if (contributions[nm] > 0) "positive" else "negative",
      rank         = i
    )
  })

  log_msg("Explained sample ", sample_index, ": ", true_label,
          ", pred=", round(pred, 3))

  list(
    sample_index   = sample_index,
    true_label     = true_label,
    prediction     = round(pred, 4),
    top_pathways   = top_pathways,
    total_pathways = length(mod$coefs),
    n_nonzero      = sum(mod$coefs != 0)
  )
}

# ── synthetic data with REAL gene IDs for pipeline testing ────────────────────

.make_synthetic_se <- function(n_samples = 200L) {
  stopifnot(
    requireNamespace("SummarizedExperiment", quietly = TRUE),
    requireNamespace("msigdbr", quietly = TRUE)
  )

  set.seed(42L)
  log_msg("Generating synthetic SE with real Ensembl IDs from msigdbr Hallmark")

  # Fetch real Ensembl gene IDs from MSigDB Hallmark
  msig_args <- if ("collection" %in% names(formals(msigdbr::msigdbr))) {
    list(species = "Homo sapiens", collection = "H")
  } else {
    list(species = "Homo sapiens", category = "H")
  }
  h <- do.call(msigdbr::msigdbr, msig_args)
  all_genes <- unique(h$ensembl_gene)
  all_genes <- all_genes[!is.na(all_genes) & nchar(all_genes) > 0]
  n_genes <- length(all_genes)
  log_msg("Using ", n_genes, " real Ensembl IDs from Hallmark gene sets")

  # Base expression: log-normal noise
  expr <- matrix(
    rlnorm(n_genes * n_samples, meanlog = 3, sdlog = 1.5),
    nrow = n_genes, ncol = n_samples
  )
  rownames(expr) <- all_genes
  colnames(expr) <- paste0("SAMPLE_", sprintf("%03d", seq_len(n_samples)))

  labels <- factor(rep(c("KIRC", "KIRP"), each = n_samples / 2L))

  # Inject signal: upregulate KRAS_SIGNALING_DN genes in KIRP samples
  # so the VNN actually learns a pathway-level pattern
  signal_pathways <- c("HALLMARK_KRAS_SIGNALING_DN", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
  for (pw in signal_pathways) {
    pw_genes <- h$ensembl_gene[h$gs_name == pw]
    pw_genes <- pw_genes[pw_genes %in% all_genes]
    kirp_cols <- which(labels == "KIRP")
    expr[pw_genes, kirp_cols] <- expr[pw_genes, kirp_cols] * 3.0
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays  = list(tpm = expr),
    colData = S4Vectors::DataFrame(label = labels)
  )
  rownames(se) <- all_genes

  log_msg("Synthetic SE: ", n_genes, " genes x ", n_samples, " samples, ",
          "signal injected in: ", paste(signal_pathways, collapse = ", "))
  se
}

# ── dispatch ─────────────────────────────────────────────────────────────────

HANDLERS <- list(
  load_tcga          = handle_load_tcga,
  load_custom        = handle_load_custom,
  build_pathway_map  = handle_build_pathway_map,
  build_architecture = handle_build_architecture,
  train_vnn          = handle_train_vnn,
  predict            = handle_predict,
  explain            = handle_explain
)

dispatch <- function(method, params) {
  handler <- HANDLERS[[method]]
  if (is.null(handler)) {
    stop("Unknown method: '", method, "'. Available: ",
         paste(names(HANDLERS), collapse = ", "))
  }
  handler(params)
}

# ── main loop ────────────────────────────────────────────────────────────────

main <- function() {
  log_msg("Starting VNNBio bridge...")

  suppressMessages({
    library(jsonlite)
    library(SummarizedExperiment)
  })

  vnnbio_available <- tryCatch({
    suppressMessages(library(VNNBio))
    TRUE
  }, error = function(e) {
    log_msg("WARNING: VNNBio not available (", conditionMessage(e), ")")
    log_msg("Running in STUB MODE — only load_tcga with synthetic data will work")
    FALSE
  })

  julia_ready <- FALSE
  # Julia init disabled — using pure-R glmnet path for hackathon speed
  log_msg("Julia init SKIPPED (using pure-R glmnet classification path)")

  send_response("__ready__", result = list(
    ready            = TRUE,
    vnnbio_available = vnnbio_available,
    julia_ready      = julia_ready,
    r_version        = paste0(R.version$major, ".", R.version$minor)
  ))

  log_msg("Entering command loop")
  # R's readLines on stdin does NOT block in Rscript child processes —
  # it returns character(0) immediately when no data is available.
  # Poll with a short sleep instead of breaking on empty reads.
  con <- file("stdin", open = "r")
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  repeat {
    line <- tryCatch(
      readLines(con, n = 1L, warn = FALSE),
      error = function(e) { log_msg("stdin read error: ", conditionMessage(e)); NULL }
    )
    if (is.null(line)) break  # real error — pipe broken
    if (length(line) == 0L) {
      Sys.sleep(0.05)  # no data yet, poll again
      next
    }
    if (nchar(trimws(line)) == 0L) next

    id <- "__unknown__"
    tryCatch({
      cmd <- fromJSON(line, simplifyVector = FALSE)
      id  <- cmd$id %||% "__no_id__"
      result <- dispatch(cmd$method, cmd$params %||% list())
      send_response(id, result = result)
    }, error = function(e) {
      log_msg("ERROR [", id, "]: ", conditionMessage(e))
      send_response(id, error = conditionMessage(e))
    })
  }

  log_msg("Bridge shutdown complete")
}

main()
