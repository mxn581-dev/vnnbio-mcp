suppressMessages({library(ggplot2); library(SummarizedExperiment)})
dir.create("plots", showWarnings = FALSE)

# ------------------------------------------------------------------
# PLOT 1 — Shapley pathway attribution for patient #229 (Died, P_died=0.75)
# Diverging "waterfall" of signed contributions
# ------------------------------------------------------------------
sh <- data.frame(
  pathway = c("HYPOXIA","P53_PATHWAY","PI3K_AKT_MTOR_SIGNALING","COMPLEMENT",
              "APICAL_JUNCTION","NOTCH_SIGNALING","GLYCOLYSIS","IL2_STAT5_SIGNALING",
              "APICAL_SURFACE","CHOLESTEROL_HOMEOSTASIS","REACTIVE_OXYGEN_SPECIES",
              "EPITHELIAL_MESENCHYMAL_TRANSITION","TNFA_SIGNALING_VIA_NFKB",
              "PROTEIN_SECRETION","SPERMATOGENESIS"),
  contribution = c(-38.8241,-14.6558,12.7247,10.853,10.294,10.0729,9.908,7.7724,
                   6.7027,-6.2804,6.1701,-5.3914,-4.4688,3.8104,-3.2464))
sh$dir <- ifelse(sh$contribution < 0, "Pushes toward DEATH", "Pushes toward SURVIVAL")
sh$pathway <- factor(sh$pathway, levels = sh$pathway[order(sh$contribution)])

p1 <- ggplot(sh, aes(x = contribution, y = pathway, fill = dir)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey30") +
  geom_text(aes(label = sprintf("%+.1f", contribution),
                hjust = ifelse(contribution < 0, 1.1, -0.1)), size = 3) +
  scale_fill_manual(values = c("Pushes toward DEATH" = "#c0392b",
                               "Pushes toward SURVIVAL" = "#2471a3"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.18)) +
  labs(title = "Why the model predicted DEATH for patient GSM1692079 (#229)",
       subtitle = "True label: Died   |   P(death) = 0.75   |   Shapley pathway attribution (Hallmark)",
       x = "Signed Shapley contribution  (logit scale)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
ggsave("plots/1_patient229_shapley.png", p1, width = 9.5, height = 6.2, dpi = 140)

# ------------------------------------------------------------------
# PLOT 2 — Cohort-level pathway importance (model coefficients)
# ------------------------------------------------------------------
co <- data.frame(
  pathway = c("NOTCH_SIGNALING","APICAL_SURFACE","HYPOXIA","PI3K_AKT_MTOR_SIGNALING",
              "REACTIVE_OXYGEN_SPECIES","CHOLESTEROL_HOMEOSTASIS","P53_PATHWAY",
              "APICAL_JUNCTION","GLYCOLYSIS","SPERMATOGENESIS","COMPLEMENT",
              "EPITHELIAL_MESENCHYMAL_TRANSITION","IL2_STAT5_SIGNALING",
              "PEROXISOME","PROTEIN_SECRETION"),
  importance = c(0.0949,0.0609,0.0505,0.028,0.0247,0.0219,0.0188,0.0185,0.0149,
                 0.014,0.0128,0.0126,0.0098,0.0084,0.0083),
  direction = c("positive","positive","negative","positive","positive","negative",
                "negative","positive","positive","negative","positive","negative",
                "positive","positive","positive"))
co$assoc <- ifelse(co$direction == "negative", "Higher -> DEATH", "Higher -> SURVIVAL")
co$pathway <- factor(co$pathway, levels = co$pathway[order(co$importance)])

p2 <- ggplot(co, aes(x = importance, y = pathway, fill = assoc)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("Higher -> DEATH" = "#c0392b",
                               "Higher -> SURVIVAL" = "#2471a3"), name = NULL) +
  labs(title = "Cohort-wide Hallmark pathway importance",
       subtitle = "VNN trained on 479 ICU patients (114 died / 365 survived)  |  CV AUC = 0.69",
       x = "|model coefficient|  (importance)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
ggsave("plots/2_cohort_pathway_importance.png", p2, width = 9.5, height = 6, dpi = 140)

# ------------------------------------------------------------------
# PLOT 3 — Predicted survival probability by true outcome (with patient #229 marked)
# ------------------------------------------------------------------
se <- readRDS("data/sepsis_mortality.rds")
lab <- as.character(colData(se)$label)
p <- c(0.8904,0.9094,0.7635,0.9405,0.9014,0.5378,0.5751,0.755,0.6486,0.6013,0.6876,0.865,0.5751,0.5322,0.6088,0.7086,0.7501,0.9023,0.8883,0.8863,0.9003,0.7721,0.9412,0.8697,0.8401,0.5713,0.8092,0.8564,0.7292,0.8286,0.7452,0.7373,0.9716,0.6833,0.7871,0.8403,0.8582,0.7971,0.2777,0.7517,0.8532,0.765,0.5151,0.7906,0.7476,0.765,0.7794,0.6291,0.904,0.4742,0.6271,0.5263,0.7816,0.5284,0.6638,0.8924,0.8671,0.79,0.7663,0.7306,0.9142,0.8681,0.6763,0.3999,0.8354,0.8174,0.5728,0.6865,0.9221,0.5413,0.3313,0.7574,0.8554,0.7211,0.7605,0.7094,0.8492,0.7706,0.7994,0.8076,0.3775,0.6724,0.7977,0.8342,0.6835,0.63,0.8963,0.8295,0.8864,0.3932,0.8781,0.9085,0.6294,0.6169,0.6278,0.5359,0.7655,0.8975,0.7117,0.746,0.753,0.7657,0.7808,0.8439,0.7138,0.792,0.826,0.9486,0.709,0.7297,0.6906,0.8562,0.5472,0.8264,0.756,0.4419,0.8454,0.8822,0.84,0.8622,0.9127,0.8761,0.884,0.5022,0.8711,0.7267,0.808,0.9198,0.885,0.9193,0.9005,0.7397,0.7521,0.8519,0.8819,0.7811,0.9507,0.8688,0.7004,0.6444,0.9487,0.8167,0.4072,0.9323,0.5347,0.6416,0.7313,0.8702,0.7826,0.7534,0.8699,0.7799,0.794,0.8271,0.9337,0.9305,0.838,0.8229,0.8133,0.42,0.8048,0.2545,0.6941,0.8689,0.8438,0.6678,0.6857,0.7283,0.4804,0.7805,0.8382,0.6706,0.915,0.8652,0.415,0.7551,0.9375,0.7468,0.4784,0.8672,0.9253,0.8761,0.893,0.898,0.7137,0.838,0.8714,0.826,0.3132,0.7744,0.7969,0.6221,0.7976,0.7905,0.7057,0.4877,0.8757,0.7433,0.783,0.7117,0.8673,0.9048,0.8758,0.6091,0.9003,0.8651,0.7883,0.6254,0.7535,0.5254,0.7338,0.7925,0.8503,0.6781,0.3841,0.8725,0.7101,0.582,0.844,0.8773,0.2751,0.9233,0.8831,0.9557,0.6114,0.8939,0.828,0.9722,0.2533,0.4104,0.719,0.7592,0.7991,0.9421,0.7219,0.9422,0.8514,0.7816,0.7844,0.6596,0.9845,0.7946,0.8271,0.8621,0.4151,0.7546,0.587,0.7584,0.528,0.8754,0.7757,0.7983,0.4993,0.9302,0.8206,0.5456,0.9035,0.7834,0.8574,0.9196,0.8741,0.8144,0.726,0.734,0.7843,0.708,0.6959,0.8839,0.8939,0.8461,0.9342,0.8482,0.8073,0.4498,0.6188,0.6424,0.8108,0.8217,0.8436,0.8756,0.5549,0.8445,0.8128,0.8714,0.5916,0.945,0.8897,0.9113,0.531,0.5686,0.5684,0.653,0.7433,0.8639,0.625,0.6038,0.785,0.7113,0.7071,0.5508,0.7404,0.8874,0.779,0.913,0.7793,0.9215,0.6856,0.9001,0.7039,0.9358,0.7406,0.8187,0.6766,0.8557,0.6707,0.8282,0.7499,0.5522,0.8648,0.8288,0.8268,0.8663,0.7481,0.5571,0.7416,0.8863,0.9169,0.9167,0.7468,0.8089,0.7287,0.8928,0.7738,0.6565,0.8074,0.8999,0.918,0.8956,0.3037,0.8756,0.6946,0.8042,0.8676,0.8495,0.8753,0.8163,0.901,0.6679,0.8472,0.7045,0.8185,0.7841,0.7811,0.8296,0.7683,0.7228,0.9164,0.7236,0.8735,0.7133,0.8282,0.9374,0.811,0.6181,0.9251,0.5532,0.689,0.7884,0.7495,0.6723,0.7066,0.9198,0.6219,0.5563,0.8427,0.4954,0.8124,0.7945,0.8901,0.5264,0.8614,0.8584,0.6789,0.7883,0.7559,0.7892,0.7195,0.8132,0.5921,0.9682,0.7854,0.7478,0.5794,0.6021,0.434,0.8681,0.8221,0.8293,0.7586,0.9593,0.8541,0.7512,0.7794,0.606,0.5838,0.6619,0.7901,0.7549,0.8414,0.9574,0.7254,0.61,0.9309,0.8663,0.4732,0.4408,0.9481,0.825,0.7753,0.8971,0.668,0.8257,0.8605,0.5644,0.8509,0.9014,0.9483,0.9274,0.7957,0.6508,0.6221,0.9068,0.8059,0.9274,0.6831,0.8976,0.9495,0.7984,0.9066,0.8144,0.874,0.8218,0.8854,0.7525,0.7513,0.8665,0.623,0.9337,0.8631,0.8332,0.7331,0.8709,0.7864,0.6215,0.7449,0.7923,0.8055,0.8143,0.8316,0.7534,0.9179,0.5802,0.4527,0.8386,0.7333,0.9051,0.5702,0.4866,0.8839,0.6013,0.6866,0.7975,0.914,0.8395,0.8602,0.7913,0.8233,0.8392,0.8812)
df <- data.frame(p_survival = p, outcome = lab)
pt <- df$p_survival[229]

p3 <- ggplot(df, aes(x = p_survival, fill = outcome)) +
  geom_histogram(binwidth = 0.04, color = "white", position = "identity", alpha = 0.6) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.5, y = Inf, label = "decision threshold", vjust = 1.5,
           hjust = -0.05, size = 3, color = "grey30") +
  geom_vline(xintercept = pt, color = "#c0392b", linewidth = 0.9) +
  annotate("text", x = pt, y = Inf, label = "patient #229", vjust = 3, hjust = 1.1,
           size = 3.2, color = "#c0392b", fontface = "bold") +
  scale_fill_manual(values = c("Died" = "#c0392b", "Survived" = "#2471a3"), name = "True outcome") +
  labs(title = "Predicted survival probability separates outcomes",
       subtitle = "Patient #229 sits deep in the low-survival region the model learned",
       x = "Model P(survival)", y = "Number of patients") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave("plots/3_probability_distribution.png", p3, width = 9.5, height = 5.2, dpi = 140)

cat("Saved 3 plots to plots/\n")
