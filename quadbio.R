library(simspec)
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(AnnotationHub)

counts_WT <- Read10X_h5("/multiome_project/WT/outs/filtered_feature_bc_matrix.h5")
counts_MT <- Read10X_h5("/multiome_project/MT/outs/filtered_feature_bc_matrix.h5")

ah <- AnnotationHub()
ensdbs <- query(ah, c("EnsDb.Hsapiens"))

ensdb_id <- ensdbs$ah_id[grep(paste0(" 98 EnsDb"), ensdbs$title)]
ensdb <- ensdbs[[ensdb_id]]

seqlevelsStyle(ensdb) <- "UCSC"
annotations <- GetGRangesFromEnsDb(ensdb = ensdb)
genome(annotations) <- "hg38"

seurat_WT <- CreateSeuratObject(counts = counts_WT$`Gene Expression`,
                                assay = "RNA",
                                project = "LM_WT")
seurat_WT[['ATAC']] <- CreateChromatinAssay(counts = counts_WT$`Peaks`,
                                            annotation = annotations,
                                            fragments = "/multiome_project/WT/outs/atac_fragments.tsv.gz",
                                            sep = c(":", "-"),
                                            genome = 'hg38')

seurat_MT <- CreateSeuratObject(counts = counts_MT$`Gene Expression`,
                                assay = "RNA",
                                project = "LM_MT")
seurat_MT[['ATAC']] <- CreateChromatinAssay(counts = counts_MT$`Peaks`,
                                            annotation = annotations,
                                            fragments = "/multiome_project/MT/outs/atac_fragments.tsv.gz",
                                            sep = c(":", "-"),
                                            genome = 'hg38')

seurat <- merge(seurat_WT, seurat_MT)

peaks <- reduce(unlist(as(c(seurat_WT@assays$ATAC@ranges,
                            seurat_MT@assays$ATAC@ranges),
                          "GRangesList")))
peakwidths <- width(peaks)
peaks <- peaks[peakwidths < 10000 & peakwidths > 20]

counts_atac_merged <- FeatureMatrix(seurat@assays$ATAC@fragments,
                                    features = peaks,
                                    cells = colnames(seurat))
seurat[['ATAC']] <- CreateChromatinAssay(counts_atac_merged,
                                         fragments = seurat@assays$ATAC@fragments,
                                         annotation = seurat@assays$ATAC@annotation,
                                         sep = c(":","-"),
                                         genome = "hg38")

library(BSgenome.Hsapiens.UCSC.hg38)
standard_chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
idx_standard_chroms <- which(as.character(seqnames(granges(seurat[['ATAC']]))) %in% standard_chroms)
seurat[["ATAC"]] <- subset(seurat[["ATAC"]],
                           features = rownames(seurat[["ATAC"]])[idx_standard_chroms])
seqlevels(seurat[['ATAC']]@ranges) <- intersect(seqlevels(granges(seurat[['ATAC']])),
                                                unique(seqnames(granges(seurat[['ATAC']]))))

seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
seurat <- NucleosomeSignal(seurat, assay = "ATAC")
seurat <- TSSEnrichment(seurat, assay = "ATAC")

VlnPlot(seurat,
        features = c("nFeature_RNA",
                     "percent.mt",
                     "nFeature_ATAC",
                     "TSS.enrichment",
                     "nucleosome_signal"),
        ncol = 5,
        pt.size = 0)

seurat <- subset(seurat,
                 subset = nFeature_RNA > 1000 &
                   nFeature_RNA < 20000 &
                   percent.mt < 30 &
                   nFeature_ATAC > 1000 &
                   nFeature_ATAC < 30000 &
                   TSS.enrichment > 1 &
                   nucleosome_signal < 2
)

DefaultAssay(seurat) <- "RNA"

seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20, reduction.name = "umap_rna", reduction.key = "UMAPRNA_")

p1 <- DimPlot(seurat, group.by = "orig.ident", reduction = "umap_rna") & NoAxes()
p2 <- FeaturePlot(seurat,
                  c("PECAM1","TWIST1","SNAI1","COL1A1"),
                  reduction = "umap_rna") & NoAxes() & NoLegend()
p1 | p2

library(simspec)
library(Pando)
seurat <- cluster_sim_spectrum(seurat,
                               label_tag = "orig.ident",
                               cluster_resolution = 0.6,
                               reduction.name = "css_rna",
                               reduction.key = "CSSRNA_")
seurat <- RunUMAP(seurat,
                  reduction = "css_rna",
                  dims = 1:ncol(Embeddings(seurat,"css_rna")),
                  reduction.name = "umap_css_rna",
                  reduction.key = "UMAPCSSRNA_")

p1 <- DimPlot(seurat, group.by = "orig.ident", reduction = "umap_css_rna") & NoAxes()
p2 <- FeaturePlot(seurat,
                  c("PECAM1","TWIST1","SNAI1","COL1A1"),
                  reduction = "umap_css_rna") & NoAxes() & NoLegend()
p1 | p2

seurat <- FindNeighbors(seurat,
                        reduction = "css_rna",
                        dims = 1:ncol(Embeddings(seurat,"css_rna"))) %>%
  FindClusters(resolution = 0.2)

DE_cl_rna <- presto::wilcoxauc(seurat, "RNA_snn_res.0.2")
top_markers <- DE_cl_rna %>%
  dplyr::filter(logFC > log(1.2) &
                  auc > 0.7 &
                  padj < 0.01 &
                  pct_in - pct_out > 30 &
                  pct_out < 30) %>%
  group_by(group) %>%
  top_n(1, wt = auc)

p1 <- DimPlot(seurat,
              group.by="RNA_snn_res.0.2",
              reduction="umap_css_rna", label=T) & NoAxes() & NoLegend()
p2 <- FeaturePlot(seurat,
                  features = unique(top_markers$feature),
                  reduction="umap_css_rna",
                  order = T,
                  ncol=3) & NoAxes() & NoLegend()
(p1 | p2) + patchwork::plot_layout(widths = c(2,3))


DefaultAssay(seurat) <- "ATAC"
seurat <- FindTopFeatures(seurat, min.cutoff = 50)

seurat <- RunTFIDF(seurat, method = 1)

seurat <- RunSVD(seurat, n = 50)



p1 <- ElbowPlot(seurat, ndims = 30, reduction="lsi")
p2 <- DepthCor(seurat, n = 30)
p1 | p2

seurat <- RunUMAP(seurat,
                  reduction = "lsi",
                  dims = 2:30,
                  reduction.name = "umap_atac",
                  reduction.key = "UMAPATAC_")
p1 <- DimPlot(seurat,
              group.by = "orig.ident",
              reduction = "umap_atac") & NoAxes()
p2 <- FeaturePlot(seurat,
                  c("PECAM1","TWIST1","SNAI1","COL1A1"),
                  reduction = "umap_atac") & NoAxes() & NoLegend()
p1 | p2

gene_act <- GeneActivity(seurat)
seurat[['RNA_inferred']] <- CreateAssayObject(gene_act) %>% NormalizeData()

DefaultAssay(seurat) <- "RNA_inferred"
beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
p3 <- FeaturePlot(seurat,
                  c("PECAM1","TWIST1","SNAI1","COL1A1"),
                  reduction = "umap_atac",
                  cols = beach_colscheme(30)) & NoAxes() & NoLegend()
p1 | p2 | p3

DefaultAssay(seurat) <- "ATAC"

library(simspec)
seurat <- cluster_sim_spectrum(seurat,
                               label_tag = "orig.ident",
                               use_dr = "lsi",
                               dims_use = 2:30,
                               cluster_resolution = 0.6,
                               reduction.name = "css_atac",
                               reduction.key = "CSSATAC_")
seurat <- RunUMAP(seurat,
                  reduction = "css_atac",
                  dims = 1:ncol(Embeddings(seurat,"css_atac")),
                  reduction.name = "umap_css_atac",
                  reduction.key = "UMAPCSSATAC_")

library(harmony)
seurat <- RunHarmony(seurat,
                     group.by.vars = "orig.ident",
                     reduction = "lsi",
                     dims.use = 2:30,
                     max.iter.harmony = 50,
                     reduction.save = "harmony_atac")
seurat <- RunUMAP(seurat,
                  reduction = "harmony_atac",
                  dims = 1:ncol(Embeddings(seurat,"harmony_atac")),
                  reduction.name = "umap_harmony_atac",
                  reduction.key = "UMAPHARMONYATAC_")

seurat <- FindMultiModalNeighbors(seurat,
                                  reduction.list = list("css_rna", "css_atac"),
                                  dims.list = list(1:ncol(Embeddings(seurat,"css_rna")),
                                                   1:ncol(Embeddings(seurat,"css_atac"))),
                                  modality.weight.name = c("RNA.weight","ATAC.weight"),
                                  verbose = TRUE)


seurat <- RunUMAP(seurat, nn.name = "weighted.nn", assay = "RNA")
seurat <- FindClusters(seurat, graph.name = "wsnn", resolution = 0.2)

p1 <- UMAPPlot(seurat, group.by = "orig.ident") & NoAxes()
p2 <- UMAPPlot(seurat, group.by = "wsnn_res.0.2", label=T) & NoAxes()
p3 <- FeaturePlot(seurat,
                  c("PECAM1","TWIST1","SNAI1","COL1A1"),
                  reduction = "umap") & NoAxes() & NoLegend()
p1 | p2 | p3

seurat$celltype <- setNames(rep(c("endothelial","mural","intermediate"), c(4,4,1)),
                            c(c(0,4,5,8),c(1,2,6,7),3))[as.character(seurat$wsnn_res.0.2)]

p1 <- UMAPPlot(seurat, group.by = "celltype", label=T) & NoAxes()
p2 <- FeaturePlot(seurat,
                  c("PDGFB","CLDN5","PDGFRB","COL5A1"),
                  order=T,
                  reduction = "umap") & NoAxes() & NoLegend()
p1 | p2

library(presto)

DefaultAssay(seurat) <- "RNA"
DE_ct <- wilcoxauc(seurat, "celltype", seurat_assay = "RNA")
top_markers_ct <- DE_ct %>%
  dplyr::filter(abs(logFC) > log(1.2) &
                  padj < 0.01 &
                  auc > 0.65 &
                  pct_in - pct_out > 30 &
                  pct_out < 20) %>%
  group_by(group) %>%
  top_n(10, wt = auc)

top_markers_ct

DefaultAssay(seurat) <- "ATAC"
DA_ct <- wilcoxauc(seurat, "celltype", seurat_assay = "ATAC")
top_peaks_ct <- DA_ct %>%
  dplyr::filter(abs(logFC) > log(1.1) &
                  padj < 0.01 &
                  auc > 0.55) %>%
  group_by(group) %>%
  top_n(100, wt = auc)

marker_peak_ct %>% top_n(5, wt=auc)

library(BSgenome.Hsapiens.UCSC.hg38)

seurat <- RegionStats(seurat,
                      genome = BSgenome.Hsapiens.UCSC.hg38)
seurat <- LinkPeaks(seurat,
                    peak.assay = "ATAC",
                    expression.assay = "RNA",
                    genes.use = top_markers_ct$feature)

DefaultAssay(seurat) <- "ATAC"
p1 <- CoveragePlot(seurat,
                   region = "RASGRP3",
                   features = "RASGRP3",
                   group.by = "celltype",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
p2 <- CoveragePlot(seurat,
                   region = "COL6A3",
                   features = "COL6A3",
                   group.by = "celltype",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
patchwork::wrap_plots(p1, p2, ncol = 1)


library(TFBSTools)
library(JASPAR2022)

pfm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
df_pfm <- data.frame(t(sapply(pfm, function(x)
  c(id=x@ID, name=x@name, symbol=ifelse(!is.null(x@tags$symbol),x@tags$symbol,NA)))))

seurat <- AddMotifs(seurat, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)

seurat <- RunChromVAR(seurat, genome = BSgenome.Hsapiens.UCSC.hg38)

library(chromVAR)
DefaultAssay(seurat) <- "chromvar"
DA_motifs_ct <- wilcoxauc(seurat, group_by = "celltype", seurat_assay = "chromvar") %>%
  mutate(symbol = setNames(ifelse(is.na(df_pfm$symbol), df_pfm$name, df_pfm$symbol),
                           df_pfm$id)[feature])

enriched_motifs_ct <- DA_motifs_ct %>%
  filter(padj < 0.01 & auc > 0.7) %>%
  group_by(group)
top_motifs_ct <- top_n(enriched_motifs_ct, 3, wt=auc)

tfs <- read.table("/multiome_project/Homo_sapiens_TF.txt", sep="\t", header=T)

tf_motifs_ct <- enriched_motifs_ct %>%
  filter(symbol %in% tfs$Symbol)
marker_tfs_ct <- DE_ct %>%
  filter(feature %in% tfs$Symbol &
           abs(logFC) > log(1.2) &
           padj < 0.01 &
           auc > 0.65 &
           pct_in - pct_out > 20) %>%
  inner_join(tf_motifs_ct,
             by = c("feature" = "symbol"),
             suffix = c("_tf","_motif")) %>%
  filter(group_tf == group_motif)

top_tfs_ct <- group_by(marker_tfs_ct, group_tf) %>%
  top_n(3, wt = auc_motif)

beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
DefaultAssay(seurat) <- "RNA"
p1 <- FeaturePlot(seurat,
                  top_tfs_ct$feature,
                  reduction = "umap",
                  order=T,
                  cols=beach_colscheme(30),
                  ncol=6) & NoAxes() & NoLegend()
DefaultAssay(seurat) <- "chromvar"
p2 <- FeaturePlot(seurat,
                  top_tfs_ct$feature_motif,
                  reduction = "umap",
                  order=T,
                  cols=bluered_colscheme(30),
                  ncol=6) & NoAxes() & NoLegend()
p1 / p2

library(Pando)

seurat <- initiate_grn(seurat,
                       regions=phastConsElements20Mammals.UCSC.hg38,
                       rna_assay = "RNA", peak_assay = "ATAC")

seurat <- find_motifs(seurat,
                      pfm = Pando::motifs,
                      motif_tfs = Pando::motif2tf,
                      genome = BSgenome.Hsapiens.UCSC.hg38)

library(doParallel)
registerDoParallel(20)
seurat <- infer_grn(seurat,
                    parallel = T,
                    tf_cor = 0.05,
                    method="glm",
                    family="gaussian",
                    scale=F,
                    verbose=T)
grn <- seurat@grn@networks$glm_network@coefs %>%
  filter(padj < 0.01)

grn

positive_regulons <- split(grn$target[grn$estimate>0], grn$tf[grn$estimate>0])
positive_regulons <- positive_regulons[lengths(positive_regulons) > 10]
negative_regulons <- split(grn$target[grn$estimate<0], grn$tf[grn$estimate<0])
negative_regulons <- negative_regulons[lengths(negative_regulons) > 10]

DefaultAssay(seurat) <- "RNA"
mod_act_pos <- AddModuleScore(seurat,
                              features = positive_regulons,
                              name = "regulon_")@meta.data
mod_act_pos <- mod_act_pos[,grep("^regulon_", colnames(mod_act_pos))] %>%
  setNames(paste0(names(positive_regulons),"(+)"))
mod_act_neg <- AddModuleScore(seurat,
                              features = negative_regulons,
                              name = "regulon_")@meta.data
mod_act_neg <- mod_act_neg[,grep("^regulon_", colnames(mod_act_neg))] %>%
  setNames(paste0(names(negative_regulons),"(-)"))

seurat[['regulon']] <- CreateAssayObject(data = t(cbind(mod_act_pos, mod_act_neg)))

DefaultAssay(seurat) <- "RNA"
p1 <- FeaturePlot(seurat,
                  top_tfs_ct$feature,
                  reduction = "umap",
                  cols = beach_colscheme(30),
                  order = T,
                  ncol = 6) & NoAxes() & NoLegend()
DefaultAssay(seurat) <- "regulon"
p2 <- FeaturePlot(seurat,
                  features = c(intersect(paste0(top_tfs_ct$feature,"(+)"), rownames(seurat)),
                               intersect(paste0(top_tfs_ct$feature,"(-)"), rownames(seurat))),
                  reduction = "umap",
                  cols = bluered_colscheme(30),
                  order = T,
                  ncol = 6) & NoAxes() & NoLegend()
(p1 / p2) + patchwork::plot_layout(height = c(1,2))



