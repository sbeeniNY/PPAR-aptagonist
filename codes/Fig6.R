# Fig6.R — scRNA-seq: UMAP, stacked violin plots & GSEA (Adipose)
#
# Fig6l: UMAP colored by cell type
# Fig6m: Stacked violin plots (3 panels, genes split into groups)
# Fig6n-p: Adipose, Cholangiocyte, Hepatocyte GSEA dot plot
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

set.seed(1111)

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
# DEG: Adipocyte (Treated vs Control) via FindMarkers
# =============================================================================

# Idents(W) <- "celltype.stim"
# Adipo_response <- FindMarkers(W, ident.1 = c("Adipocyte1_Treated", "Adipocyte2_Treated"), ident.2 = c("Adipocyte1_Control", "Adipocyte2_Control"), min.pct =0,logfc.threshold=-Inf)
# Adipo1_response <- FindMarkers(W, ident.1 = c("Adipocyte1_Treated"), ident.2 = c("Adipocyte1_Control"), min.pct =0,logfc.threshold=-Inf)
# Adipo2_response <- FindMarkers(W, ident.1 = c("Adipocyte2_Treated"), ident.2 = c("Adipocyte2_Control"), min.pct =0,logfc.threshold=-Inf)
# Hepatocyte_response <- FindMarkers(W, ident.1 = c("Hepatocyte_Treated"), ident.2 = c("Hepatocyte_Control"), min.pct =0,logfc.threshold=-Inf)
# Cholangiocyte_response <- FindMarkers(W, ident.1 = c("Cholangiocyte_Treated"), ident.2 = c("Cholangiocyte_Control"), min.pct =0,logfc.threshold=-Inf)

# write.table(Adipo_response, "Adipo_response.txt", sep = "\t", quote = FALSE)
# write.table(Adipo1_response, "Adipo1_response.txt", sep = "\t", quote = FALSE)
# write.table(Adipo2_response, "Adipo2_response.txt", sep = "\t", quote = FALSE)
# write.table(Hepatocyte_response, "Hepatocyte_response.txt", sep = "\t", quote = FALSE)
# write.table(Cholangiocyte_response, "Cholangiocyte_response.txt", sep = "\t", quote = FALSE)


# =============================================================================
# Fig6m — gene level differential expression analysis heatmap
# =============================================================================
sub_W <- subset(W, subset = celltype %in%
                  c("Hepatocyte", "Cholangiocyte", "Adipocyte1", "Adipocyte2"))
Idents(sub_W) <- "celltype"

gene_list <- c("PPARG","CD36","FABP4","LEP","LIPE","MGLL","RBP7","ABCA1","SIK3","GSTM5","APOE","NR1H3","APOC2","APOA2","SOAT1")
# stim_color  <- c("#d3d3d3", "#da4b35")

# p_vln1 <- Stacked_VlnPlot(sub_W, features = gene_list[1:5],
#                            split.by = "stim", colors_use = stim_color,
#                            x_lab_rotate = TRUE)
# ggsave(file.path(fig_dir, "Fig6.panelM.1.png"), p_vln1,
#        width = 5, height = 10, dpi = 300)

# p_vln2 <- Stacked_VlnPlot(sub_W, features = gene_list[6:10],
#                            split.by = "stim", colors_use = stim_color,
#                            x_lab_rotate = TRUE)
# ggsave(file.path(fig_dir, "Fig6.panelM.2.png"), p_vln2,
#        width = 5, height = 10, dpi = 300)

# p_vln3 <- Stacked_VlnPlot(sub_W, features = gene_list[11:16],
#                            split.by = "stim", colors_use = stim_color,
#                            x_lab_rotate = TRUE)
# ggsave(file.path(fig_dir, "Fig6.panelM.3.png"), p_vln3,
#        width = 5, height = 10, dpi = 300)

dea_files <- c(
  Hepatocyte    = "Hepatocyte_response.txt",
  Cholangiocyte = "Cholangiocyte_response.txt",
  Adipocyte1    = "Adipo1_response.txt",
  Adipocyte2    = "Adipo2_response.txt"
)

dea_stats <- lapply(names(dea_files), function(ct) {
  path <- file.path(data_dir, dea_files[[ct]])
  df   <- read.table(path, sep = "\t", header = TRUE, row.names = 1)
  df %>%
    rownames_to_column("gene") %>%
    filter(gene %in% gene_list) %>%
    transmute(
      celltype     = ct,
      gene         = gene,
      avg_log2FC   = avg_log2FC,
      p_val        = p_val,
      p_val_adj    = p_val_adj,
      sig          = ifelse(p_val < 0.05, "*", "")
    )
}) %>% bind_rows()


celltype_order <- c("Hepatocyte","Cholangiocyte","Adipocyte1","Adipocyte2")

df_heat <- dea_stats %>%
  mutate(
    gene     = factor(gene,     levels = gene_list),
    celltype = factor(celltype, levels = rev(celltype_order))
  )

p_heat <- ggplot(df_heat, aes(x = gene, y = celltype)) +
  geom_tile(aes(fill = avg_log2FC), color = "grey85", linewidth = 0.3) +
  scale_fill_gradientn( 
    colors = c("#3B82F6", "white", "#EF4444"),
    values = rescale(c(-3, 0, 8)), 
    limits = c(-3, 8),
    oob = scales::squish,
    name = "avg log2FC\n(Treated-Control)",
    na.value = "grey90") +
  geom_text(data = df_heat %>% filter(sig != ""),
            aes(label = sig), size = 5, vjust = 0.8, fontface = "bold") +
  labs(caption = "* p_val < 0.05  (Seurat FindMarkers)") +
  theme_bw(base_size = 11) +
  theme(
    axis.title      = element_blank(),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y     = element_text(size = 10),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    plot.caption    = element_text(size = 8, color = "grey50"),
    legend.key.size = unit(0.45, "cm"),
    legend.title    = element_text(size = 9)
  )

ggsave("Fig6.panelM.png", p_heat, width = 9, height = 3)

# =============================================================================
# GSEA gene sets (Macaca mulatta: Hallmark / Reactome / GO:BP)
# =============================================================================
# safe_msigdbr <- function(category, subcategory) {
#   tryCatch(
#     msigdbr(species = "Macaca mulatta", category = category,
#             subcategory = subcategory),
#     error = function(e)
#       tryCatch(
#         msigdbr(species = "Macaca mulatta", category = category,
#                 subcategory = paste0(subcategory, "_LEGACY")),
#         error = function(e2) tibble(gs_name = character(), gene_symbol = character())
#       )
#   )
# }

# msig_H     <- msigdbr(species = "Macaca mulatta", category = "H")
# msig_REACT <- safe_msigdbr("C2", "CP:REACTOME")
# msig_BP    <- msigdbr(species = "Macaca mulatta", category = "C5", subcategory = "GO:BP")

# patterns <- c(
#   "ADIPO", "LIPID", "FATTY", "BETA.*OXID", "TRIGLY", "CHOLES", "BILE",
#   "PPAR", "INSULIN", "GLUCO", "GLYCOL", "GLUCONEO", "OXIDATIVE", "MITO",
#   "THERMOGEN", "BROWN", "UCP", "LIPOPROTEIN", "FOXO", "SREBP", "MTOR", "AMPK",
#   "PEROXISOME", "FAT.*METAB", "CARNITINE", "ADIPOGEN"
# )
# rx <- paste0("(", paste(patterns, collapse = "|"), ")")

# gsets_df <- bind_rows(
#   filter(msig_H,     str_detect(gs_name, rx)),
#   filter(msig_REACT, str_detect(gs_name, rx)),
#   filter(msig_BP,    str_detect(gs_name, rx))
# ) |> distinct(gs_name, gene_symbol)
# gsets <- split(gsets_df$gene_symbol, gsets_df$gs_name)

# ranks_adipo <- {
#   r <- setNames(adipo_deg$avg_log2FC, rownames(adipo_deg))
#   r <- r[!is.na(r)]
#   r <- tapply(r, names(r), function(z) z[which.max(abs(z))]) |> unlist()
#   sort(r, decreasing = TRUE)
# }

# gsea_adipo <- fgsea(pathways = gsets, stats = ranks_adipo,
#                     minSize = 5, maxSize = 500, nperm = 10000) |>
#   arrange(padj, desc(NES))

# =============================================================================
# Fig6n-p — Curated GSEA bubble plot (cell-type-specific GSEA)
# Pre-computed GSEA TSV files from upstream analysis 
# =============================================================================

gsea_adipo_all    <- read.table(file.path(data_dir,"GSEA_SC_Adipo_filtered.tsv"),        sep = "\t", header = TRUE)
gsea_cholangiocyte <- read.table(file.path(data_dir,"GSEA_SC_Cholangiocyte_filtered.tsv"), sep = "\t", header = TRUE)
gsea_hepatocyte    <- read.table(file.path(data_dir,"GSEA_SC_Hepatocyte_filtered.tsv"),    sep = "\t", header = TRUE)

entries_adipo <- c(
  "GOBP_POSITIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS",
  "GOBP_POSITIVE_REGULATION_OF_TRIGLYCERIDE_CATABOLIC_PROCESS",
  "GOBP_CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS",
  "GOBP_REGULATION_OF_VERY_LOW_DENSITY_LIPOPROTEIN_PARTICLE_REMODELING",
  "GOBP_CELLULAR_RESPONSE_TO_LIPID",
  "GOBP_REGULATION_OF_LIPOPROTEIN_METABOLIC_PROCESS",
  "GOBP_LIPID_DROPLET_FUSION")
entries_chol <- c(
  "GOBP_CELLULAR_RESPONSE_TO_LIPOPROTEIN_PARTICLE_STIMULUS",
  "GOBP_PLASMA_LIPOPROTEIN_PARTICLE_CLEARANCE",
  "GOBP_RESPONSE_TO_LIPOPROTEIN_PARTICLE",
  "REACTOME_NR1H3_NR1H2_REGULATE_GENE_EXPRESSION_LINKED_TO_CHOLESTEROL_TRANSPORT_AND_EFFLUX",
  "REACTOME_TRIGLYCERIDE_BIOSYNTHESIS",
  "GOBP_LIPOPROTEIN_BIOSYNTHETIC_PROCESS"
)
entries_hep <- c(
  "REACTOME_MITOCHONDRIAL_IRON_SULFUR_CLUSTER_BIOGENESIS",
  "GOBP_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_MITOCHONDRIAL_MEMBRANE",
  "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
  "REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
  "GOBP_REGULATION_OF_AUTOPHAGY_OF_MITOCHONDRION",
  "REACTOME_TRANSPORT_OF_BILE_SALTS_AND_ORGANIC_ACIDS_METAL_IONS_AND_AMINE_COMPOUNDS"
)

df_adipo <- gsea_adipo_all    %>% filter(pathway %in% entries_adipo)
df_chol  <- gsea_cholangiocyte %>% filter(pathway %in% entries_chol)
df_hep   <- gsea_hepatocyte    %>% filter(pathway %in% entries_hep)

# ── 공통 size 스케일: 3개 데이터의 -log10(pval) 전체 범위로 limits 통일 ──────
size_limit <- range(-log10(bind_rows(df_adipo, df_chol, df_hep)$pval),
                    na.rm = TRUE)

gsea_dot_theme <- theme_light() +
  theme(
    axis.title       = element_text(size = 12),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 10)
  )

make_dot <- function(df, title = NULL) {
  ggplot(df, aes(x = NES, y = reorder(pathway, NES),
                 size = -log10(pval), color = NES)) +
    geom_point(alpha = 0.7) +
    scale_size_continuous(
      name   = "-log10(p-value)",
      range  = c(2, 10),
      limits = size_limit          # 3개 플랏 동일 스케일
    ) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = "NES"
    ) +
    labs(x = "NES", y = NULL, title = title) +
    gsea_dot_theme
}

p_dot_adipo <- make_dot(df_adipo, "Adipocyte")
p_dot_chol  <- make_dot(df_chol,  "Cholangiocyte")
p_dot_hep   <- make_dot(df_hep,   "Hepatocyte")

ggsave(file.path(fig_dir, "Fig6.panelN.png"),        p_dot_adipo, width = 12, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "Fig6.panelO.png"), p_dot_chol,  width = 12, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "Fig6.panelP.png"),    p_dot_hep,   width = 12, height = 4, dpi = 300)


sessionInfo()
