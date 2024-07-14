library(shiny)
library(knitr)
library(tidyverse)
library(AnnotationDbi)
library(edgeR)
library(DESeq2)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86) 
library(bslib)
library(plotly)
library(ComplexHeatmap)

## set theme
theme_set(theme_bw())

## load gene count matrix
count_df <- read_tsv("https://raw.githubusercontent.com/joowkim/ngs-data-pipeline/master/star_read_cnt.tsv")

## get gene annotation
ensdb_genes <- genes(EnsDb.Hsapiens.v86) |> 
  as.data.frame() |>
  dplyr::select(gene_id, gene_name)

cnt_df <- left_join(count_df, ensdb_genes, by = "gene_id") |>
  dplyr::rename(gene_symbol = gene_name) |>
  relocate(gene_symbol, .after = gene_id)


## meta data
meta_df <- data.frame(
  stringsAsFactors = FALSE,
                              samp_name = c("HER2.1","HER2.2",
                                            "HER2.3","NBS1","NBS2","NBS3",
                                            "Non.TNBC1","Non.TNBC2","Non.TNBC3",
                                            "TNBC1","TNBC2","TNBC3"),
                                subtype = c("HER","HER","HER",
                                            "Normal","Normal","Normal","Non.TNBC",
                                            "Non.TNBC","Non.TNBC","TNBC",
                                            "TNBC","TNBC"),
                              samp_type = c("Tissue","Tissue",
                                            "Tissue","Organoid","Organoid",
                                            "Organoid","Tissue","Tissue","Tissue",
                                            "Tissue","Tissue","Tissue")
                     )


## total # reads
cnt_df_long <- cnt_df |> 
  dplyr::select(-gene_id, -gene_symbol) |>
  pivot_longer(cols = everything(), names_to = "sample_name", values_to = "cnt") |>
  group_by(sample_name) |> summarise(total_reads = sum(cnt))

plt_total_num_reads <- ggplot(cnt_df_long, aes(x=sample_name, y=total_reads, fill=sample_name)) + 
  geom_bar(stat = "identity") + 
  scale_fill_viridis_d() +
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  scale_y_continuous(name="total # reads", labels = scales::comma)

total_read_df <- cnt_df |> 
  dplyr::select(-gene_id, -gene_symbol) |> 
  summarise(across(where(is.numeric),sum)) |> 
  pivot_longer(everything(), names_to ="sample_name", values_to = "total_reads") |> 
  mutate(across(total_reads, ~ format (., big.mark=",")))

stopifnot(all(colnames(cnt_df |> dplyr::select(-gene_id, -gene_symbol)) == (meta_df$samp_name)))

cnt_all_df <- cnt_df |> 
  dplyr::select(-gene_symbol) |>
  data.frame() |> 
  column_to_rownames(var = "gene_id")

gene_names_df <- data.frame(row.names = count_df$gene_id)

# factor litter_num
dds <- DESeqDataSetFromMatrix(countData = cnt_all_df,
                              colData = meta_df,
                              design = ~ + subtype)

mcols(dds) <- gene_names_df

## Filtering to remove low counts
# https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.knit.md
y <- DGEList(counts=cnt_all_df, genes = cnt_df |> dplyr::select(gene_id, gene_symbol) |> column_to_rownames("gene_id"), samples = meta_df, group = meta_df$subtype) 
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
dds <- dds[keep,]

## define pca helper function - get pc scores and eigen values/vectors
# https://rdrr.io/bioc/DESeq2/src/R/plots.R
pca_helper <- function(object, ntop=500) {
  # object is either vst or rlog
  
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  return (pca)
}

## vst & regularized log
rlog_data<- rlog(dds, blind=FALSE)
# rlogdata<- rlog(dds, blind=FALSE)

# get pca_obj
pca_obj <- pca_helper(rlog_data, ntop = Inf)
# get pca scores
pca_scores <- pca_obj$x %>% data.frame()

# get eigen values
pca_eigen_val <- round((summary(pca_obj)$importance[2,] * 100),1)

# merge pca_socres and meta_df
merged_df <- pca_scores %>% rownames_to_column("samp_name") 
pca_merged_df <- left_join(merged_df, meta_df, by = 'samp_name')

## plot pca
## pca pc1/pc2
# see https://stackoverflow.com/questions/56547256/show-3-factors-ggplot-geom
pcaplot_pc1_pc2 <- ggplot(pca_merged_df, aes(x=PC1, y=PC2, color=subtype, shape=samp_type)) + 
  geom_point(size=5.5) +
  # coord_fixed() + 
  # scale_size_manual(values=c(7,5))+
  ggrepel::geom_text_repel(label=pca_merged_df$samp_name, box.padding = .5) +
  xlab(paste0("PC1: ", pca_eigen_val[1], " % var expl")) + 
  ylab(paste0("PC2: ", pca_eigen_val[2], " % var expl")) +
  ggtitle("PCA plot, all genes, rlog applied")  +
  theme_bw() +
  # theme(aspect.ratio=1) + 
  theme(text = element_text(size = 30))
  # scale_size_manual(values = c(8,5))

## scree plot
scree_plot <- pcaExplorer::pcascree(pca_obj, type = 'pev')
# print(scree_plot)

## spearman correlation - Michael Love prefers the ecluidean distance than correlation distance "the correlation between normalized, transformed sample counts is typically going to be very high unless you compare across different cell types" from Michael Love.
# https://support.bioconductor.org/p/100620/

# corrs <- cor(assay(vstdata), method = "spearman")
sampleDists <- dist(t(assay(rlog_data)))
sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vstdata$Genotype, vstdata$Treat, sep="-")
# colnames(sampleDistMatrix) <- paste(vstdata$Genotype, vstdata$Treat, sep="-")
# corr.dists <- as.dist(1 - corrs)

# ann <-data.frame(Treat = meta_df$Treat, Genotype = meta_df$Genotype)
# ann <- meta_df %>% dplyr::select(-samp_name, -Group)
ann <- data.frame(subtype = meta_df$subtype, samp_type = meta_df$samp_type)
# ann <-data.frame(Treat = meta_df$Treat) #, Treat = meta_df$Treat)
ann_colors <- list(
  subtype = c("HER" = "forestgreen", "Normal" = "gold", "Non.TNBC" = "cyan1", "TNBC" = 'yellow'),
  samp_type = c("Tissue" = 'cyan1', "Organoid" = "darkgoldenrod1")
  # Sex = c("Female" = "green", "Male" = "yellow")
  
)
heat_anno <- HeatmapAnnotation(df = ann,
                               col = ann_colors
)

# distance method = euclidean distances based on vst
hmap <- pheatmap(sampleDistMatrix,
                 # clustering_distance_columns = sampleDists,
                 # clustering_distance_rows = sampleDists,
                 # rect_gp = gpar(col = "white", lwd = 2), 
                 # # row_names_gp = gpar(fontsize = 5), 
                 # row_dend_side = 'left',
                 column_title = "Heatmap of the sample-to-sample based on Euclidean distances", 
                 heatmap_legend_param = list(
                   title = "Euclidean dist"),
                 top_annotation = heat_anno,
)

design <- model.matrix(~0 + subtype, data = y$samples)
colnames(design) <- colnames(design) %>% str_remove("subtype")
# print(design)

y <- estimateDisp(y, design, robust=TRUE)
# plotBCV(y, main=paste0("BCV plot"))
fit <- glmQLFit(y, design, robust=TRUE)
# plotQLDisp(fit, main=paste0("QLDisp plot"))
# paste0("common BCV: ", sqrt(y$common.dispersion))

## contrast
contrasts <- makeContrasts(
  HER_vs_Normal = HER-Normal,
  Non.TNBC_vs_Normal = Non.TNBC-Normal,
  TNBC_vs_Normal = TNBC-Normal,
  TNBC_vs_Non.TNBC = TNBC-Non.TNBC,
  TNBC_vs_HER = TNBC-HER,
  # KI effect within the Control group: 
  # Group2_vs_Group1 = Group2-Group1,
  # KI effect within the Drug group: 
  # Group4_vs_Group3 = Group4-Group3,
  # Drug treatment effect within the WT group: 
  # Group3_vs_Group1 = Group3-Group1,
  # Drug treatment effect within the KI group: 
  # Group4_vs_Group2 = Group4-Group2,
  # Overall differences between the Con and Drug group:
  # Group4_3_vs_Group1_2 = (Group4+Group3)/2 - (Group1+Group2)/2,
  # Overall differences between the WT and KI group:
  # Group2_4_vs_Group1_3 = (Group2+Group4)/2 - (Group1+Group3)/2,
  # BVsT=(groupLUNG_B+groupBRAIN_B)/2-(groupLUNG_T+groupBRAIN_T)/2,
  # LungVsBrain=(groupLUNG_B+groupLUNG_T)/2-(groupBRAIN_B+groupBRAIN_T)/2,
  # BVsT_Lung=groupLUNG_B-groupLUNG_T,
  # BVsT_Brain=groupBRAIN_B-groupBRAIN_T,
  # LungVsBrain_B=groupLUNG_B-groupBRAIN_B, 
  # LungVsBrain_T=groupLUNG_T-groupBRAIN_T, 
  levels=colnames(design))
rownames(contrasts) <- gsub("group", "", rownames(contrasts))


qlf <- lapply(rlang::set_names(colnames(contrasts), colnames(contrasts)), function(contrast){
  glmQLFTest(fit, contrast=contrasts[,contrast])
})

res <- lapply(qlf, function(contrast) topTags(contrast, n = Inf))

# lapply(res, function(contrast) table(as.data.frame(contrast)$FDR < 0.05))
# lapply(res, function(contrast) {
#   table(as.data.frame(contrast) %>% 
#           dplyr::mutate(signif=FDR < 0.05, dir=ifelse(logFC>0,"up","dwn")) %>% 
#           dplyr::select(signif, dir))
# })

out_objs <- list(dgelist=y, fit=fit, qlf=qlf, res=res)

volcano <- lapply(rlang::set_names(names(res),names(res)), function(contrast){
  toptable <- res[[contrast]]$table
  #toptable
  EnhancedVolcano::EnhancedVolcano(toptable=toptable, x="logFC", y="FDR", 
                                   lab=toptable$gene_symbol, title=contrast, pCutoff=0.05, FCcutoff = 1, ylab="-log10(FDR)",
                                   subtitle = "", legendDropLevels=FALSE, caption = paste0("total = ", nrow(toptable), " genes"),
                                   legendLabels = c("NS", expression(Log[2] ~ FC), "FDR", expression(FDR ~ and ~ log[2] ~ FC)),
                                   ylim = c(0, max(-log10(toptable$FDR), na.rm = TRUE) + 0.5),
                                   xlim = c(-max(abs(toptable$logFC))-0.5, max(abs(toptable$logFC))+0.5))
})

plt_volcano <- patchwork::wrap_plots(volcano) + patchwork::plot_annotation(title = 'Volcano plots') + patchwork::plot_layout(guides = "collect") & theme(legend.position = 'bottom', title = element_text(size=18, face='bold'))

HER_vs_Normal <- res$HER_vs_Normal |> data.frame() |> dplyr::arrange(FDR) |> head(30)

Non.TNBC_vs_Normal <- res$Non.TNBC_vs_Normal |> data.frame() |> dplyr::arrange(FDR) |> head(30)

TNBC_vs_Normal <- res$TNBC_vs_Normal |> data.frame() |> dplyr::arrange(FDR) |> head(30)

TNBC_vs_Non.TNBC  <- res$TNBC_vs_Non.TNBC  |> data.frame() |> dplyr::arrange(FDR) |> head(30)

TNBC_vs_HER <- res$TNBC_vs_HER |>  data.frame() |> dplyr::arrange(FDR) |> head(30)
## define UI

ui <- fluidPage(
  navbarPage("Breast Cancer Gene expression data analysis",
             tabPanel("Introduction",
                      mainPanel(
                        h5("Breast cancer bulk RNA-seq data from human cell lines were obtained from NCBI GEO data (PRJNA227137). 
                          The raw data was generated from an experiment with a one-factor design and four levels (Normal, HER2, Non-TNBC, and TNBC). 
                          There were three replicates for each group, resulting in 12 sequenced samples. The raw data was processed using bioinformatics tools, including STAR."
                        ), # h4 end
                        p(),
                        h3("Meta data information"),# h5 end
                        ## DT::DTOutput("total_num_reads_table"),
                        DT::DTOutput("meta_df"),
                        plotlyOutput("plt_total_num_reads")
                      ) # mainPanel end
             ), # tabPanel end
             tabPanel("Explanatory data analysis - PCA & Hclustering",
                      mainPanel(
                        h5("Here, we perform PCA and hierarchical clustering to identify the major source of variation. 
                           We assume that biological replicates from the same group or subtype will cluster together, indicating that the unexplained variance would be relatively small unless there is a batch effect or outliers present. 
                           We utilize all genes and apply a regularized-logarithm transformation to stabilize the variance across the mean. 
                           Our results show that PC1 clearly separates normal samples from breast cancer samples. 
                           Additionally, we observe that the sample type (organoids and tissue) is confounded with subtype."
                        ),
                        p(),
                        plotOutput("plt_pca"),
                        p(),
                        h3("Sample to sample distance"),
                        plotOutput("plt_dist_hmap")
                      ) # mainPanel end
             ), # tabPanel end
             tabPanel("DGE analysis",
               mainPanel(
                 h5(
                   "We perform differential expressed gene analysis using edgeR, https://pubmed.ncbi.nlm.nih.gov/19910308/. 
                   We can see that there are many differently expressed genes between breast cancer types and normal samples in volcano plots.
                   Differential expressed genes can be investigated further using GSEA or ORA analysis to look into associations between them and specific pathways.
                   "
                 ),
                 plotOutput("plt_volcano", height = "700px"),
                 h4("Top 30 sig DEG, HER vs Normal"),
                 DT::DTOutput("HER_vs_Normal"),
                 h4("Top 30 sig DEG, Non.TNBC_vs_Normal"),
                 DT::DTOutput("Non.TNBC_vs_Normal"),
                 h4("Top 30 sig DEG, TNBC_vs_Normal"),
                 DT::DTOutput("TNBC_vs_Normal"),
                 h4("Top 30 sig DEG, TNBC_vs_Non.TNBC"),
                 DT::DTOutput("TNBC_vs_Non.TNBC"),
                 h4("Top 30 sig DEG, TNBC_vs_HER"),
                 DT::DTOutput("TNBC_vs_HER")
               )
             )# tabPanel end
  ) # navbarPage end
) # fluidPage end
                        
                      
server <- function(input, output, session) {
  output$plt_total_num_reads <- renderPlotly(plt_total_num_reads)
  output$total_num_reads_table <- DT::renderDT(total_read_df, caption = htmltools::tags$caption( style = 'caption-side: top; text-align: center; color:black;  font-size:200% ;', 'Meta data') )
  output$meta_df <- DT::renderDT(meta_df)
  output$plt_pca <- renderPlot(pcaplot_pc1_pc2)
  output$plt_dist_hmap <- renderPlot(hmap)
  output$plt_volcano <- renderPlot(plt_volcano)
  output$HER_vs_Normal <- DT::renderDT(HER_vs_Normal)
  output$Non.TNBC_vs_Normal <- DT::renderDT(Non.TNBC_vs_Normal)
  output$TNBC_vs_Normal <- DT::renderDT(TNBC_vs_Normal)
  output$TNBC_vs_Non.TNBC  <- DT::renderDT(TNBC_vs_Non.TNBC)
  output$TNBC_vs_HER <- DT::renderDT(TNBC_vs_HER)
  
} # server end
shinyApp(ui, server)

# Run the app ----
shinyApp(ui = ui, server = server)