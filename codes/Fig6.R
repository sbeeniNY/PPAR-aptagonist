# Fig6.R — scRNA-seq: UMAP, stacked violin plots & GSEA (Adipose)
#
# Fig6l: UMAP colored by cell type
# Fig6m: Stacked violin plots (3 panels, genes split into groups)
# Fig6n: Adipose GSEA dot plot
# Fig6o: Enrichment plot (top NES pathway)
#
# Run this script from the repository root.
# Outputs are written to figures/

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(fgsea)
  library(msigdbr)
})

set.seed(1234)

data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# =============================================================================
# Load Seurat object
# =============================================================================
W <- readRDS(file.path(data_dir, "PPARG.scRNA.rds"))

# =============================================================================
# Fig6l — UMAP colored by cell type
# =============================================================================
pal <- c(
  "NK"          = "#85C7FF", "T"            = "#5FA9E6", "B"           = "#C8C2F0",
  "Neutrophil"  = "#6FB5D6", "Erythroid"    = "#9E9E9E",
  "Hepatocyte"  = "#E85C5C", "Cholangiocyte"= "#F5A1A1", "Kupper"      = "#E7C86E",
  "Adipocyte1"  = "#9FE07C", "Adipocyte2"   = "#79D6A0"
)

p_umap <- DimPlot(W, reduction = "umap.harmony", group.by = "celltype",
                  repel = TRUE, cols = pal) + NoLegend()
ggsave(file.path(fig_dir, "Fig6.panelL.png"), p_umap,
       width = 6, height = 6, dpi = 300)

# =============================================================================
# Fig6m — Stacked violin plots (Adipocyte & Hepatocyte subsets)
# =============================================================================
sub_W <- subset(W, subset = celltype %in%
                  c("Hepatocyte", "Cholangiocyte", "Adipocyte1", "Adipocyte2"))
Idents(sub_W) <- "celltype"

gene_list   <- c("PPARG","RXRA","SREBF2","ACACA","CD36",
                 "FABP4","PNPLA2","MGLL","ACADM","PDK4",
                 "PPARGC1A","INSR","AKT2","NAMPT","APOB","ABCA1")
stim_color  <- c("#d3d3d3", "#da4b35")

p_vln1 <- Stacked_VlnPlot(sub_W, features = gene_list[1:5],
                           split.by = "stim", colors_use = stim_color,
                           x_lab_rotate = TRUE)
ggsave(file.path(fig_dir, "Fig6.panelM.1.png"), p_vln1,
       width = 5, height = 10, dpi = 300)

p_vln2 <- Stacked_VlnPlot(sub_W, features = gene_list[6:10],
                           split.by = "stim", colors_use = stim_color,
                           x_lab_rotate = TRUE)
ggsave(file.path(fig_dir, "Fig6.panelM.2.png"), p_vln2,
       width = 5, height = 10, dpi = 300)

p_vln3 <- Stacked_VlnPlot(sub_W, features = gene_list[11:16],
                           split.by = "stim", colors_use = stim_color,
                           x_lab_rotate = TRUE)
ggsave(file.path(fig_dir, "Fig6.panelM.3.png"), p_vln3,
       width = 5, height = 10, dpi = 300)

# =============================================================================
# DEG: Adipocyte (Treated vs Control) via FindMarkers
# =============================================================================
Idents(W) <- "celltype.stim"
adipo_deg <- FindMarkers(
  W,
  ident.1 = c("Adipocyte1_Treated"),
  ident.2 = c("Adipocyte1_Control"),
  min.pct = 0, logfc.threshold = -Inf
)

# =============================================================================
# GSEA gene sets (Macaca mulatta: Hallmark / Reactome / GO:BP)
# =============================================================================
safe_msigdbr <- function(category, subcategory) {
  tryCatch(
    msigdbr(species = "Macaca mulatta", category = category,
            subcategory = subcategory),
    error = function(e)
      tryCatch(
        msigdbr(species = "Macaca mulatta", category = category,
                subcategory = paste0(subcategory, "_LEGACY")),
        error = function(e2) tibble(gs_name = character(), gene_symbol = character())
      )
  )
}

msig_H     <- msigdbr(species = "Macaca mulatta", category = "H")
msig_REACT <- safe_msigdbr("C2", "CP:REACTOME")
msig_BP    <- msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:BP")

patterns <- c(
  "ADIPO", "LIPID", "FATTY", "BETA.*OXID", "TRIGLY", "CHOLES", "BILE",
  "PPAR", "INSULIN", "GLUCO", "GLYCOL", "GLUCONEO", "OXIDATIVE", "MITO",
  "THERMOGEN", "BROWN", "UCP", "LIPOPROTEIN", "FOXO", "SREBP", "MTOR", "AMPK",
  "PEROXISOME", "FAT.*METAB", "CARNITINE", "ADIPOGEN"
)
rx <- paste0("(", paste(patterns, collapse = "|"), ")")

gsets_df <- bind_rows(
  filter(msig_H,     str_detect(gs_name, rx)),
  filter(msig_REACT, str_detect(gs_name, rx)),
  filter(msig_BP,    str_detect(gs_name, rx))
) |> distinct(gs_name, gene_symbol)
gsets <- split(gsets_df$gene_symbol, gsets_df$gs_name)

ranks_adipo <- {
  r <- setNames(adipo_deg$avg_log2FC, rownames(adipo_deg))
  r <- r[!is.na(r)]
  r <- tapply(r, names(r), function(z) z[which.max(abs(z))]) |> unlist()
  sort(r, decreasing = TRUE)
}

gsea_adipo <- fgsea(pathways = gsets, stats = ranks_adipo,
                    minSize = 5, maxSize = 500, nperm = 10000) |>
  arrange(padj, desc(NES))

# =============================================================================
# Fig6n — Adipose GSEA dot plot (top5 + bottom5 NES)
# =============================================================================
tb <- bind_rows(
  arrange(gsea_adipo, desc(NES)) |> slice_head(n = 5),
  arrange(gsea_adipo, NES)       |> slice_head(n = 5)
) |>
  mutate(pathway = str_trunc(pathway, 65)) |>
  arrange(NES) |>
  mutate(pathway = factor(pathway, levels = pathway))

p_gsea_dot <- ggplot(tb, aes(x = NES, y = pathway,
                              size = -log10(pval), color = NES)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "-log10(p-value)", range = c(2, 10)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "NES") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(fig_dir, "Fig6.panelN.png"), p_gsea_dot,
       width = 10, height = 4.5, dpi = 300)

# =============================================================================
# Fig6o — Enrichment plot (GOBP_POSITIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS)
# =============================================================================
enrich_id <- "GOBP_POSITIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS"

enrich_row <- gsea_adipo |> filter(pathway == enrich_id)
enrich_nes <- round(enrich_row$NES, 2)
enrich_fdr <- signif(enrich_row$padj, 3)

p_enrichment <- fgsea::plotEnrichment(gsets[[enrich_id]], ranks_adipo) +
  ggtitle(enrich_id) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("NES=%.2f  FDR=%.3g", enrich_nes, enrich_fdr),
           hjust = 1, vjust = 1.6, size = 3)

ggsave(file.path(fig_dir, "Fig6.panelO.png"), p_enrichment,
       width = 7, height = 4, dpi = 300)

sessionInfo()
