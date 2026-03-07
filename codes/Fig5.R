# Fig5.R — Bulk RNA-seq: Volcano plots & GSEA (Adipose and Liver)
#
# Fig5c: Adipose volcano plot
# Fig5d: Liver volcano plot
# Fig5e: Adipose GSEA bubble plot
# Fig5f: Liver GSEA bubble plot
# Fig5g: GSEA upregulated terms Venn diagram (Adipose vs Liver)
# Fig5h: Shared upregulated pathway dot plot
#
# Run this script from the repository root.
# Outputs are written to figures/

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(fgsea)
  library(msigdbr)
  library(ggVennDiagram)
})

set.seed(1234)

data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# =============================================================================
# Load & prepare count matrix
# =============================================================================
counts <- read.delim(file.path(data_dir, "count_w_geneName.tsv"),
                     check.names = FALSE, row.names = 1)

samples <- colnames(counts)
meta <- tibble(sample = samples) |>
  mutate(
    Origin = if_else(str_detect(sample, "^Liver"), "Liver", "Adipo"),
    batch  = gsub("(.*)_(.*)_(.*)", "\\2", sample),
    Treat  = if_else(str_detect(sample, "ctrl"), "ctrl", "treat")
  ) |>
  mutate(across(c(Origin, batch, Treat), factor)) |>
  mutate(Origin = relevel(Origin, ref = "Adipo"),
         Treat  = relevel(Treat,  ref = "ctrl"))
rownames(meta) <- meta$sample

cnt  <- round(as.matrix(counts))
keep <- rowSums(cnt >= 10) >= ceiling(ncol(cnt) / 2)
cnt  <- cnt[keep, , drop = FALSE]

# =============================================================================
# DESeq2: per-organ differential expression (~ batch + Treat)
# =============================================================================
run_deseq <- function(origin_label) {
  idx <- meta$Origin == origin_label
  dds <- DESeqDataSetFromMatrix(cnt[, idx], droplevels(meta[idx, ]),
                                design = ~ batch + Treat)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, name = "Treat_treat_vs_ctrl")
  tibble(gene   = rownames(res),
         log2FC = as.numeric(res$log2FoldChange),
         pval   = as.numeric(res$pvalue),
         padj   = as.numeric(res$padj))
}

res_adipo <- run_deseq("Adipo")
res_liver <- run_deseq("Liver")

# =============================================================================
# GSEA: metabolic gene sets (Hallmark / KEGG / Reactome / GO:BP)
# =============================================================================
safe_msigdbr <- function(category, subcategory) {
  tryCatch(
    msigdbr(species = "Mus musculus", category = category,
            subcategory = subcategory),
    error = function(e)
      tryCatch(
        msigdbr(species = "Mus musculus", category = category,
                subcategory = paste0(subcategory, "_LEGACY")),
        error = function(e2) tibble(gs_name = character(), gene_symbol = character())
      )
  )
}

patterns <- paste0("(",  paste(c(
  "ADIPO","LIPID","FATTY","BETA.*OXID","TRIGLY","CHOLES","BILE",
  "PPAR","INSULIN","GLUCO","GLYCOL","GLUCONEO","OXIDATIVE","MITO",
  "THERMOGEN","BROWN","UCP","LIPOPROTEIN","FOXO","SREBP","MTOR","AMPK",
  "PEROXISOME","FAT.*METAB","CARNITINE","ADIPOGEN"
), collapse = "|"), ")")

gsets_df <- bind_rows(
  filter(msigdbr(species = "Mus musculus", category = "H"),          str_detect(gs_name, patterns)),
  filter(safe_msigdbr("C2", "CP:KEGG"),                              str_detect(gs_name, patterns)),
  filter(safe_msigdbr("C2", "CP:REACTOME"),                          str_detect(gs_name, patterns)),
  filter(msigdbr(species = "Mus musculus", category = "C5",
                 subcategory = "GO:BP"),                              str_detect(gs_name, patterns))
) |> distinct(gs_name, gene_symbol)
gsets <- split(gsets_df$gene_symbol, gsets_df$gs_name)

make_ranks <- function(res_tbl) {
  r <- setNames(res_tbl$log2FC, res_tbl$gene)
  r <- r[!is.na(r)]
  r <- tapply(r, names(r), function(z) z[which.max(abs(z))]) |> unlist()
  sort(r, decreasing = TRUE)
}

run_gsea <- function(res_tbl, out_tsv) {
  fg <- fgsea(pathways = gsets, stats = make_ranks(res_tbl),
              minSize = 5, maxSize = 500, nperm = 10000) |>
    arrange(padj, desc(NES)) |>
    mutate(leadingEdge = vapply(leadingEdge,
                                function(x) paste(x, collapse = ";"), ""))
  write.table(fg, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  fg
}

fg_adipo <- run_gsea(res_adipo, file.path(data_dir, "GSEA_fgsea_Adipo_filtered.tsv"))
fg_liver <- run_gsea(res_liver, file.path(data_dir, "GSEA_fgsea_Liver_filtered.tsv"))

extract_leading_edge <- function(fg, top_n = 5) {
  le <- bind_rows(arrange(fg, desc(NES)) |> slice_head(n = top_n),
                  arrange(fg, NES)       |> slice_head(n = top_n))$leadingEdge
  unique(unlist(strsplit(paste(le, collapse = ";"), ";")))
}

le_adipo <- extract_leading_edge(fg_adipo)
le_liver <- extract_leading_edge(fg_liver)

# =============================================================================
# functions
# =============================================================================
plot_volcano <- function(res_tbl, le_genes, organ_title,
                         lfc_thr = 1.0, alpha_p = 0.01,
                         xlim = NULL, out_png) {
  df <- res_tbl |>
    mutate(neg_log10p = -log10(pmax(pval, .Machine$double.xmin)),
           sig = case_when(
             !is.na(pval) & pval < alpha_p & log2FC >=  lfc_thr ~ "Up",
             !is.na(pval) & pval < alpha_p & log2FC <= -lfc_thr ~ "Down",
             TRUE ~ "not sig."
           ))
  lab_df <- filter(df, gene %in% le_genes, sig != "not sig.") |>
    arrange(pval, desc(abs(log2FC)))

  p <- ggplot(df, aes(x = log2FC, y = neg_log10p)) +
    geom_point(aes(color = sig), size = 1.4, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(values = c("Up" = "#EF4444", "Down" = "#3B82F6",
                                  "not sig." = "grey80"),
                       breaks = c("Up", "Down", "not sig.")) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(alpha_p),
               linetype = "dashed", linewidth = 0.4) +
    labs(title = organ_title,
         x = expression(Log[2] ~ "fold change"),
         y = expression(-Log[10] ~ "(P-value)"),
         color = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank())
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  if (nrow(lab_df) > 0)
    p <- p + geom_text_repel(data = lab_df, aes(label = gene),
                              size = 3, max.overlaps = Inf,
                              box.padding = 0.5, point.padding = 0.25,
                              min.segment.length = 0, segment.size = 0.3,
                              seed = 1234)
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
  invisible(p)
}

plot_gsea_bubble <- function(fg, organ_title, out_png,
                             top_n = 5, trunc_len = 65) {
  tb <- bind_rows(arrange(fg, desc(NES)) |> slice_head(n = top_n),
                  arrange(fg, NES)       |> slice_head(n = top_n)) |>
    mutate(pathway = str_trunc(pathway, trunc_len)) |>
    arrange(NES) |>
    mutate(pathway = factor(pathway, levels = pathway))
  p <- ggplot(tb, aes(x = NES, y = pathway,
                      size = -log10(pval), color = NES)) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(name = "-log10(p-value)", range = c(2.5, 15)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red",
                          midpoint = 0, name = "NES") +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL,
         title = organ_title) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(out_png, p, width = 10, height = 4.5, dpi = 300)
  invisible(p)
}

# =============================================================================
# Fig5c — Adipose volcano plot
# =============================================================================
plot_volcano(res_adipo, le_adipo, organ_title = "Adipose",
             lfc_thr = 1.0, alpha_p = 0.01, xlim = NULL,
             out_png = file.path(fig_dir, "Fig5.panelC.png"))

# =============================================================================
# Fig5d — Liver volcano plot
# =============================================================================
plot_volcano(res_liver, le_liver, organ_title = "Liver",
             lfc_thr = 1.0, alpha_p = 0.01, xlim = c(-6, 6),
             out_png = file.path(fig_dir, "Fig5.panelD.png"))

# =============================================================================
# Fig5e — Adipose GSEA bubble plot
# =============================================================================
plot_gsea_bubble(fg_adipo, organ_title = "GSEA — Adipose",
                 out_png = file.path(fig_dir, "Fig5.panelE.png"))

# =============================================================================
# Fig5f — Liver GSEA bubble plot
# =============================================================================
plot_gsea_bubble(fg_liver, organ_title = "GSEA — Liver",
                 out_png = file.path(fig_dir, "Fig5.panelF.png"))

# =============================================================================
# Fig5g — Venn diagram: shared upregulated GSEA terms (P < 0.05, NES > 0)
# =============================================================================
up_adipo <- fg_adipo |> filter(pval < 0.05, NES > 0) |> pull(pathway)
up_liver <- fg_liver |> filter(pval < 0.05, NES > 0) |> pull(pathway)

p_venn <- ggVennDiagram(list(Liver = up_liver, Adipose = up_adipo),
                        label = "both", label_alpha = 0,
                        set_color = c("#4DBBD5", "#E64B35"),
                        set_size  = 4) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "count") +
  scale_color_manual(values = c("#4DBBD5", "#E64B35"), guide = "none") +
  labs(subtitle = "Upregulated terms (GSEA: P < 0.05, NES > 0)") +
  theme(plot.subtitle = element_text(size = 9))

ggsave(file.path(fig_dir, "Fig5.panelG.png"), p_venn,
       width = 5, height = 4, dpi = 300)

# =============================================================================
# Fig5h — Shared upregulated pathway dot plot
# =============================================================================
shorten_pathway <- function(x) {
  x |>
    str_replace("^HALLMARK_", "H: ") |>
    str_replace("^KEGG_",     "K: ") |>
    str_replace("^REACTOME_", "R: ") |>
    str_replace("^GOBP_",     "GO: ") |>
    str_replace_all("_", " ")
}

common_paths <- intersect(up_adipo, up_liver)
dot_df <- bind_rows(
  fg_adipo |> filter(pathway %in% common_paths) |> mutate(Organ = "Adipose"),
  fg_liver |> filter(pathway %in% common_paths) |> mutate(Organ = "Liver")
) |> mutate(label = shorten_pathway(pathway))

path_order <- dot_df |>
  group_by(label) |> summarise(mean_NES = mean(NES), .groups = "drop") |>
  arrange(mean_NES) |> pull(label)
dot_df$label <- factor(dot_df$label, levels = path_order)

p_dot <- ggplot(dot_df, aes(x = NES, y = label,
                             color = Organ, size = -log10(pval))) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c("Adipose" = "#E64B35", "Liver" = "#4DBBD5")) +
  scale_size_continuous(name = "-log10(p-value)") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Fig5.panelH.png"), p_dot,
       width = 7, height = 5, dpi = 300)

sessionInfo()
