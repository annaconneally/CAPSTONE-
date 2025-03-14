
# ======
# Step 1: Load Required Libraries and Analyse Seurat Object
# ======
# Set the path to the .robj file
file_path_all_neuron_reads <- "/Users/annaconneally/Desktop/R_datasets/remote_memory_data/all_neuron_reads.Robj"
file_path_non_neuronal_cells <- "/Users/annaconneally/Desktop/R_datasets/remote_memory_data/TRAP2_non_neuronal_cells.Robj"
load(file_path_all_neuron_reads)
load(file_path_non_neuronal_cells)

# Load the file and libraries
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Inspect the Seurat Object
print(all_para_reads_filt)  # object structure
head(all_para_reads_filt@meta.data)  # View metadata

# Set Cell Type Identities & Isolate Microglia
Idents(all_para_reads_filt) <- all_para_reads_filt$celltype
table(Idents(all_para_reads_filt))  # Confirm cell types

#========
# STEP: 2 Isolating Microglia 
#=========
# This creates a seurat object called microglia, filtering from object all_para_reads_filt 
microglia <- subset(all_para_reads_filt, celltype == "3_Microglia_Il1a+" | celltype == "4_Microglia_Il1a-")
#removing cells that belong to a sample (orig.ident) with fewer than 3 cells.
microglia<- subset(microglia, cells = colnames(microglia)[table(microglia$orig.ident)[microglia$orig.ident] >= 3])
head(microglia)

# normalise - Accounts for differences in sequencing depth between cells.
microglia <- NormalizeData(microglia, normalization.method = "LogNormalize", scale.factor = 10000)
microglia <- FindVariableFeatures(microglia, selection.method = "vst", nfeatures = 2000)

# scale
all.genes <- rownames(microglia)
microglia<- ScaleData(microglia, features = all.genes)
microglia[["RNA"]]$scale.data

# top10hvg
microglia_top10_hvgenes <- head(VariableFeatures(microglia), 10)
plot1<- VariableFeaturePlot(microglia)
plot2 <- LabelPoints(plot = plot1, points = microglia_top10_hvgenes, repel = TRUE)
plot1 + plot2

# Visualse QC metrics of microglia
microglia[["percent.mt"]] <- PercentageFeatureSet(microglia, pattern = "^mt-")
VlnPlot(microglia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(microglia, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(microglia, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# PCAanalysis
microglia <- RunPCA(microglia, features = VariableFeatures(object = microglia))
print(microglia[["pca"]], dims = 1:05, nfeatures = 5) #test different PC's to choose most representative of biological variation
DimPlot(microglia, reduction = "pca")

#=======
# STEP 3 CREATE CLUSTERS & DEFINE CLUSTER MARKERS
#=======
# findclusters
microglia <- FindNeighbors(microglia, dims = 1:10)
microglia<- FindClusters(microglia, resolution = 0.5) #test different resolutions to choose optimal 
microglia <- RunUMAP(microglia, dims = 1:10)
DimPlot(microglia, reduction = "umap")
DimPlot(microglia, group.by = "orig.ident")

# UMAP colored by sample identity
DimPlot(microglia, group.by = "orig.ident") +
  labs(
    title = "UMAP Projection Colored by Sample Identity",
    x = "UMAP_1 (Primary Gene Expression Variation)",
    y = "UMAP_2 (Secondary Gene Expression Variation)",
    color = "Sample Identity") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10))
    
# Cell Count Analysis 
#Convert to a data frame for visualisation
cluster_counts <- table(microglia$seurat_clusters)
print(cluster_counts)
cluster_counts_df <- as.data.frame(cluster_counts)
colnames(cluster_counts_df) <- c("Cluster", "Cell_Count")
# Calculate summary statistics for total expression counts
summary_statistics <- summary(rowSums(microglia@assays$RNA@counts))
print(summary_statistics)
library(ggplot2)
# Plot the number of cells in each cluster
ggplot(cluster_counts_df, aes(x = Cluster, y = Cell_Count, fill = Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Count Distribution Across Microglial Clusters",
       x = "Cluster",
       y = "Number of Cells")
# Sum total transcript counts per cell
cell_counts <- colSums(microglia@assays$RNA@counts)
# Combine with cluster metadata
cell_metadata <- data.frame(Cluster = microglia$seurat_clusters, Total_Counts = cell_counts)
# Compute per-cluster summary statistics
cluster_summaries <- aggregate(Total_Counts ~ Cluster, data = cell_metadata, summary)
print(cluster_summaries)

# Find differentially expressed features (cluster biomarkers)
cluster0.markers<-FindMarkers(microglia, ident.1=0)
cluster1.markers<-FindMarkers(microglia, ident.1=1)
cluster2.markers<-FindMarkers(microglia, ident.1=2)
cluster3.markers<-FindMarkers(microglia, ident.1=3)
head(cluster0.markers, n = 5)

# Find markers for every cluster compared to all remaining cells
microglia.markers <- FindAllMarkers(microglia, min.pct = .1)
microglia.markers %>%
  group_by(cluster) %>%
  dplyr::filter(abs(avg_log2FC) > 0.25) 
markers<-read.csv("filtered_microglia_markers.csv",stringsAsFactors = FALSE)
# To retain only statistically significant genes (p_val_adj ≤ 0.05)  save them in a new file (downstream)
downstream <- microglia.markers %>% 
  dplyr::filter(p_val_adj <= 0.05)
# Save for downstream GO analysis
write.csv(downstream, "downstream_genes.csv", row.names = FALSE)

#==========
# STEP 4 : LOOPING IN T5 MICE AND GROUPING MICE BY IDENTITY AND CONDITION
#==========
# Add a new column to store experimental conditions
microglia$condition <- NA  #Initialise

# Assign conditions based on `orig.ident` or other identifiers
microglia$condition[microglia$orig.ident %in% c("FC", "T1")] <- "FC"
microglia$condition[microglia$orig.ident %in% c("Homecage", "T5")] <- "Homecage"
microglia$condition[microglia$orig.ident %in% c("Fearonly", "T4")] <- "Fearonly"
microglia$condition[microglia$orig.ident %in% c("Context", "T2", "T3")] <- "Context"
# UMAP grouped by condition
DimPlot(microglia, reduction = "umap", group.by = "condition")
# Split UMAP by condition
DimPlot(microglia, reduction = "umap", split.by = "condition")

# Add the 'seurat_clusters' column to the metadata
microglia$seurat_clusters <- Idents(microglia)
head(microglia@meta.data)
#=========
# STEP 5 : Propellor | Cell Abundace Analysis Across Condition
#=========
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(multcomp)  # Tukey's HSD test
library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare data for propeller analysis
meta_data <- microglia@meta.data
meta_data <- meta_data %>% dplyr::select(orig.ident, seurat_clusters, condition)
colnames(meta_data) <- c("sample_id", "cluster", "condition")
meta_data$cluster <- as.factor(meta_data$cluster)
meta_data$condition <- as.factor(meta_data$condition)

# Perform propeller analysis
propeller_results <- propeller(clusters = meta_data$cluster, 
                               sample = meta_data$sample_id, 
                               group = meta_data$condition)

# Convert propeller results to long format for visualisation
prop_data <- data.frame(
  Cluster = rownames(propeller_results), 
  Context = propeller_results$PropMean.Context,
  `Fear Recall` = propeller_results$PropMean.FC,  # Changed "FC" to "Fear Recall"
  Fearonly = propeller_results$PropMean.Fearonly,
  Homecage = propeller_results$PropMean.Homecage)

prop_long <- pivot_longer(prop_data, cols = -Cluster, names_to = "Condition", values_to = "Proportion")

# Normalise proportions so that they sum to 100% **within each condition**
prop_long <- prop_long %>%
  group_by(Condition) %>%
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()

# Convert to percentage for visualization
prop_long$Percentage <- round(prop_long$Proportion * 100, 1)

# Statistical Testing
anova_results <- aov(Proportion ~ Condition + Cluster, data = prop_long)
summary(anova_results)

# Tukey's HSD post-hoc test
tukey_results <- TukeyHSD(anova_results, "Cluster")
print(tukey_results)

# Kruskal-Wallis Test (Non-parametric)
kruskal_results <- kruskal.test(Proportion ~ Cluster, data = prop_long)
print(kruskal_results)

# Organise and print statistical summary table
stat_summary <- data.frame(
  Test = c("ANOVA - Condition", "ANOVA - Cluster", "Kruskal-Wallis - Cluster"),
  p_value = c(1.00000, 0.00167, kruskal_results$p.value),
  Significant = c("No", ifelse(0.00167 < 0.05, "Yes", "No"), ifelse(kruskal_results$p.value < 0.05, "Yes", "No")))
print(stat_summary)

# Step 3: Define muted cluster colors based on your request
muted_cluster_colors <- c("0" = "#66A5D9",  # Muted Blue
                          "1" = "#E78F62",  # Muted Coral
                          "2" = "#66C2A5",  # Muted Green
                          "3" = "#C390D4")  # Muted Purple

# Step 4: Generate graph including context
ggplot(prop_long, aes(x = Condition, y = Proportion, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_fill(vjust = 0.5), color = "black", size = 5) +
  scale_fill_manual(values = muted_cluster_colors, name = "Cluster") +
  labs(
    title = "Proportion of Clusters Across Experimental Conditions (Including Context)",
    x = "Condition", 
    y = "Proportion", 
    fill = "Cluster") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    axis.title.x = element_text(face = "bold", size = 14),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold and larger y-axis title
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12))

# Step 5: Generate Stacked Graph excluding context
prop_long_no_context <- prop_long %>% filter(Condition != "Context")

ggplot(prop_long_no_context, aes(x = Condition, y = Proportion, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = paste0(Percentage, "%")), position = position_fill(vjust = 0.5), color = "black", size = 5) +
  scale_fill_manual(values = muted_cluster_colors, name = "Cluster") +
  labs(
    title = "Proportion of Clusters Across Experimental Conditions (Excluding Context)",
    x = "Condition", 
    y = "Proportion", 
    fill = "Cluster") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
    axis.text.y = element_text(face = "bold"),  # Bold y-axis labels
    axis.title.x = element_text(face = "bold", size = 14),  # Bold and larger x-axis title
    axis.title.y = element_text(face = "bold", size = 14),  # Bold and larger y-axis title
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12))

# match UMAP to colours: # Load required library
library(Seurat)
library(ggplot2)

# Define muted cluster colors
muted_cluster_colors <- c("0" = "#66A5D9",  # Muted Blue
                          "1" = "#E78F62",  # Muted Coral
                          "2" = "#66C2A5",  # Muted Green
                          "3" = "#C390D4")  # Muted Purple

# Ensure clusters are stored as factors
microglia$seurat_clusters <- factor(microglia$seurat_clusters, levels = c("0", "1", "2", "3"))

# Generate UMAP with custom colors
DimPlot(microglia, reduction = "umap", cols = muted_cluster_colors) +
  labs(title = "UMAP Projection of Microglial Clusters",
       x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12))

# Generate UMAP colored by Sample Identity
DimPlot(microglia, group.by = "orig.ident") +
  labs(title = "UMAP Projection Colored by Sample Identity",
       x = "UMAP 1", y = "UMAP 2", color = "Sample Identity") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12))

# Step 10: Calculate cell counts per condition across clusters
cell_counts <- table(meta_data$cluster, meta_data$condition)

# Convert the cell counts table into a data frame
cell_counts <- as.data.frame(cell_counts)
colnames(cell_counts) <- c("Cluster", "Condition", "Freq")

# Now, 'Cluster' and 'Condition' columns available for plotting
ggplot(cell_data, aes(x = Cluster, y = Mean, color = Condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), 
                position = position_dodge(0.7), width = 0.2) + 
  geom_point(data = cell_counts, 
             aes(x = Cluster, y = Freq), 
             color = "black", position = position_jitter(width = 0.2)) +  # corrected jitter positioning
  labs(title = "Cell Abundance Across Clusters with SEM", x = "Cluster", y = "Cell Count") +
  theme_minimal() +
  scale_y_continuous(expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Visualisation 
# Group by mouse, condition, and cluster
mouse_cluster_counts <- microglia@meta.data %>%
  group_by(mouse, orig.ident, seurat_clusters) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# Rename columns for clarity
colnames(mouse_cluster_counts) <- c("Mouse", "Condition", "Cluster", "Cell_Count")

# Print summary
print(mouse_cluster_counts)

# Assign correct conditions to T-groups
mouse_cluster_counts <- mouse_cluster_counts %>%
  mutate(Condition = case_when(
    Condition %in% c("FC", "T1") ~ "FC",
    Condition %in% c("Homecage", "T5") ~ "Homecage",
    Condition %in% c("Fearonly", "T4") ~ "Fearonly",
    Condition %in% c("T2", "T3", "Context") ~ "Context",
    TRUE ~ Condition  # Keep other labels unchanged))
cell_summary <- mouse_cluster_counts %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean = mean(Cell_Count),
    SEM = sd(Cell_Count) / sqrt(n()),
    .groups = "drop")
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Ensure zero counts are included
mouse_cluster_counts <- complete(mouse_cluster_counts, Condition, Cluster, fill = list(Cell_Count = 0))

# Compute mean and SEM for each cluster-condition combination
cell_summary <- mouse_cluster_counts %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean = mean(Cell_Count),
    SEM = sd(Cell_Count) / sqrt(n()),
    .groups = "drop")

# Define custom condition colors
condition_colors <- c(
  "Fear Recall" = "#E78F62",  
  "Homecage"   = "#66C2A5",  
  "Fearonly"   = "#66A5D9",  
  "Context"    = "#C390D4" )

# Create the bar plot including the "Context" condition with jittered dots
ggplot(cell_summary, aes(x = as.factor(Cluster), y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "black") +  # Add bar borders
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.7), width = 0.3) +
  geom_point(data = mouse_cluster_counts, 
             aes(x = as.factor(Cluster), y = Cell_Count, fill = Condition),
             shape = 21,  # Circle with border
             size = 3,
             stroke = 1.2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +  
  scale_fill_manual(values = condition_colors, name = "Condition") +
  labs(title = "Microglial Cell Counts Across Conditions (with Context)",
       x = "Cluster", 
       y = "Mean Cell Count ± SEM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12))

# Filter out the "Context" condition from the dataset
mouse_cluster_counts_no_context <- mouse_cluster_counts %>% 
  filter(Condition != "Context")

# Recompute the summary statistics without the "Context" condition
cell_summary_no_context <- mouse_cluster_counts_no_context %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean = mean(Cell_Count),
    SEM = sd(Cell_Count) / sqrt(n()),
    .groups = "drop")

# Create the bar plot excluding the "Context" condition with jittered dots
ggplot(cell_summary_no_context, aes(x = as.factor(Cluster), y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7, color = "black") +  
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), position = position_dodge(width = 0.7), width = 0.3) +  
  geom_point(data = mouse_cluster_counts_no_context, 
             aes(x = as.factor(Cluster), y = Cell_Count, fill = Condition),
             shape = 21,
             size = 3,
             stroke = 1.2, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7)) +
  scale_fill_manual(values = condition_colors, name = "Condition") +
  labs(title = "Microglial Cell Counts Across Conditions (without Context)",
       x = "Cluster", 
       y = "Mean Cell Count ± SEM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 12))
  
#=======
# STEP 6 GENE ONTOLOGY ANALYSIS 
#=======
#Gene Ontology (GO) enrichment analysis & Gene Set Enrichment Analysis (GSEA) separately for each cluster using the (DEGs) from downstream file.
#GO Enrichment Analysis (enrichGO) - Identifies enriched Gene Ontology terms.
#Gene Set Enrichment Analysis (GSEA) (GSEA) - Performs ranking-based pathway analysis.
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "ggplot2", "dplyr"))

library(clusterProfiler)
library(org.Mm.eg.db)  # Use org.Hs.eg.db for human genes
library(enrichplot)
library(ggplot2)
library(dplyr)

downstream <- read.csv("downstream_genes.csv")  # Load DEGs file
cluster_list <- unique(downstream$cluster)  # Get unique cluster IDs

go_results_list <- list()
for (cl in cluster_list) {
  cluster_genes <- downstream %>%
    filter(cluster == cl) %>%
    pull(gene)  # Extract gene symbols for the cluster
  
  gene_ids <- bitr(cluster_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
  
  go_results <- enrichGO(
    gene         = gene_ids$ENTREZID,
    OrgDb        = org.Mm.eg.db,
    keyType      = "ENTREZID",
    ont          = "BP",  # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05)
  
  go_results_list[[paste0("Cluster_", cl)]] <- go_results
  
  # Save results for each cluster
  write.csv(as.data.frame(go_results), paste0("GO_Cluster_", cl, ".csv"), row.names=FALSE)
}

#Visualisation - replace for each cluster 
barplot(go_results_list[["Cluster_3"]], showCategory=10, title="Top GO Terms for Cluster 3")
dotplot(go_results_list[["Cluster_3"]], showCategory=10, title="GO Dotplot for Cluster 3")

#========
# STEP 7:  Gene Set Enrichment Analysis 
#========
# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(msigdbr)
library(dplyr)
library(enrichplot)
library(ggplot2)

# Initialise list to store GSEA results
gsea_results_list <- list()

# Ensure MSigDB gene sets are correctly formatted
msigdb_terms <- msigdbr(species = "Mus musculus", category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

# Convert gene symbols to Entrez IDs for MSigDB gene sets
msigdb_terms <- merge(msigdb_terms, bitr(msigdb_terms$gene_symbol, 
                                         fromType = "SYMBOL", 
                                         toType = "ENTREZID", 
                                         OrgDb = org.Mm.eg.db), 
                      by.x = "gene_symbol", 
                      by.y = "SYMBOL") %>%
  dplyr::select(gs_name, ENTREZID)

# Loop through each cluster
for (cl in unique(microglia.markers$cluster)) {
  
  cat(paste0("\nProcessing GSEA for Cluster ", cl, "...\n"))
  
  # Extract DEGs for the current cluster and rank by log fold change
  cluster_degs <- microglia.markers %>%
    filter(cluster == cl) %>%
    arrange(desc(avg_log2FC))  # Rank genes by log2FC

  # Extract gene names and log2FC values
  ranked_genes <- cluster_degs$avg_log2FC
  names(ranked_genes) <- cluster_degs$gene

  # Convert gene symbols to Entrez IDs
  gene_ids <- bitr(names(ranked_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

  # **Fix: Ensure Only Mapped Genes Are Kept**
  if (nrow(gene_ids) == 0) {
    cat(paste0("Skipping Cluster ", cl, ": No valid genes mapped.\n"))
    next  # Skip if no genes mapped
  }

  ranked_genes <- ranked_genes[gene_ids$SYMBOL]  # Keep only successfully mapped genes
  names(ranked_genes) <- gene_ids$ENTREZID  # Assign Entrez IDs as names

  # **Fix: Skip Cluster if No Valid Genes Remain**
  if (length(ranked_genes) == 0) {
    cat(paste0("Skipping Cluster ", cl, ": No valid genes remain after filtering.\n"))
    next  # Skip to the next cluster
  }

  # **Print Debugging Info 
  cat("Number of genes used in GSEA:", length(ranked_genes), "\n")
  
  # Run GSEA
  gsea_results <- tryCatch({
    GSEA(
      ranked_genes,
      TERM2GENE = msigdb_terms,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.1,
      scoreType = "pos",
      minGSSize = 5)}, error = function(e) {
    cat(paste0("Error in GSEA for Cluster ", cl, ": ", e$message, "\n"))
    return(NULL)})

  # Store results
  if (!is.null(gsea_results)) {
    gsea_results_list[[paste0("Cluster_", cl)]] <- gsea_results
  }

  # Save results if there are significant pathways
  if (!is.null(gsea_results) && nrow(gsea_results@result) > 0) {
    write.csv(as.data.frame(gsea_results@result), paste0("GSEA_Cluster_", cl, ".csv"), row.names = FALSE)
  }
}

# **Fix: Check If Any Clusters are Significant Enrichment Before Plotting**
if (length(gsea_results_list) == 0) {
  cat("\nNo significant GSEA results found for any cluster.\n")
} else {
  # plot results but code only consideres clusters with significant enrichment**
  for (cl in names(gsea_results_list)) {
      
      # **Check if the cluster contains enriched pathways**
      if (!is.null(gsea_results_list[[cl]]) && nrow(gsea_results_list[[cl]]@result) > 0) {
          
          plot_path <- paste0(cl, "_GSEA_dotplot.png")
          
          # Generate and Save Dot Plot
          dot_plot <- dotplot(gsea_results_list[[cl]], showCategory = 10) + 
                     ggtitle(paste0("Top Enriched Gene Sets (", cl, ")"))
          ggsave(plot_path, dot_plot, width = 8, height = 6, dpi = 300)
          
          print(dot_plot)  # Display dot plot inline

          # Generate and Save Enrichment Plot (Top 5 Pathways)
          gsea_top_terms <- head(gsea_results_list[[cl]]@result$ID, 5)
          enrichment_plot <- enrichplot::gseaplot2(gsea_results_list[[cl]], geneSetID = gsea_top_terms)
          ggsave(paste0(cl, "_GSEA_enrichment.png"), enrichment_plot, width = 8, height = 6, dpi = 300)
          
          print(enrichment_plot)  #Display enrichment plot inline

      } else {
          cat(paste0("\nSkipping Cluster ", cl, ": No enriched terms found.\n"))
      }
  }
}


# other graphs: library(Seurat)
library(ggplot2)
library(RColorBrewer)

# 1) Extract the first two principal components (PC1, PC2)
pca_coords <- Embeddings(microglia, "pca")[, 1:2]

# 2) Create a data frame with PC1, PC2, and Condition
plot_data <- data.frame(
  PC1 = pca_coords[, 1],
  PC2 = pca_coords[, 2],
  Condition = microglia$condition  # Ensure that 'condition' exists in your metadata
)

# 3) set colour palette using RColorBrewer's "Set1"
n_conditions <- length(unique(plot_data$Condition))
strong_palette <- brewer.pal(n = n_conditions, name = "Set1")

# 4) Plot the PCA 
ggplot(plot_data, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = strong_palette) +
  theme_minimal() +
  labs(
    title = "PCA: PC1 vs. PC2 by Condition (Strong Colours)",
    x = "Principal Component 1 (PC1)",
    y = "Principal Component 2 (PC2)")
# Analysing inter-mouse variability code for stats and visualisation 
# Install required packages (if not already installed)
install.packages(c("ggplot2", "dplyr", "tidyr", "scales"))

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)  # For percentage formatting
# Define muted cluster colours
muted_cluster_colors <- c("0" = "#66A5D9",  # Muted Blue
                          "1" = "#E78F62",  # Muted Coral
                          "2" = "#66C2A5",  # Muted Green
                          "3" = "#C390D4")  # Muted Purple

# Extract data for each condition
fc_mouse_cluster_counts <- table(microglia$mouse[microglia$condition == "FC"], microglia$seurat_clusters)
homecage_mouse_cluster_counts <- table(microglia$mouse[microglia$condition == "Homecage"], microglia$seurat_clusters)
fearonly_mouse_cluster_counts <- table(microglia$mouse[microglia$condition == "Fearonly"], microglia$seurat_clusters)

# Convert tables to data frames
fc_mouse_cluster_long <- as.data.frame(fc_mouse_cluster_counts)
homecage_mouse_cluster_long <- as.data.frame(homecage_mouse_cluster_counts)
fearonly_mouse_cluster_long <- as.data.frame(fearonly_mouse_cluster_counts)

# Rename columns
colnames(fc_mouse_cluster_long) <- c("Mouse", "Cluster", "Cell_Count")
colnames(homecage_mouse_cluster_long) <- c("Mouse", "Cluster", "Cell_Count")
colnames(fearonly_mouse_cluster_long) <- c("Mouse", "Cluster", "Cell_Count")
# Function to compute percentages per mouse
compute_percentages <- function(data) {
    data %>%
        group_by(Mouse) %>%
        mutate(Percentage = round((Cell_Count / sum(Cell_Count)) * 100, 1))
}

# Compute percentages for each condition
fc_mouse_cluster_long <- compute_percentages(fc_mouse_cluster_long)
homecage_mouse_cluster_long <- compute_percentages(homecage_mouse_cluster_long)
fearonly_mouse_cluster_long <- compute_percentages(fearonly_mouse_cluster_long)
# Function to generate stacked bar plots with percentage labels
generate_cluster_plot <- function(data, title) {
    ggplot(data, aes(x = Mouse, y = Cell_Count, fill = Cluster)) +
        geom_bar(stat = "identity", position = "fill", color = "black") +  # Stacked bar chart
        geom_text(aes(label = paste0(Percentage, "%")), 
                  position = position_fill(vjust = 0.5),  # Centers text inside bars
                  color = "white",  # White text for readability
                  size = 5) +  # Adjust font size
        scale_fill_manual(values = muted_cluster_colors, name = "Cluster") +  # Apply muted colors
        theme_minimal() +
        labs(
            title = title,
            x = "Mouse",
            y = "Proportion of Cells",
            fill = "Cluster"
        ) +
        scale_y_continuous(labels = scales::percent_format()) +  # Show y-axis as percentages
        theme(
            axis.text.x = element_text(face = "bold", size = 12),
            axis.text.y = element_text(face = "bold", size = 12),
            axis.title.x = element_text(face = "bold", size = 14),
            axis.title.y = element_text(face = "bold", size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(face = "bold", size = 12)
        )
}
fc_plot <- generate_cluster_plot(fc_mouse_cluster_long, "Proportion of Clusters Across FC Mice")
print(fc_plot)
homecage_plot <- generate_cluster_plot(homecage_mouse_cluster_long, "Proportion of Clusters Across Homecage Mice")
print(homecage_plot)
fearonly_plot <- generate_cluster_plot(fearonly_mouse_cluster_long, "Proportion of Clusters Across Fearonly Mice")
print(fearonly_plot)




