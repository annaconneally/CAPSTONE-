# CAPSTONE-
Final year undergraduate capstone project 2024/2025

# Step 1: Load Required Libraries & Files

#Set the path to the .robj file
file_path_all_neuron_reads <- "/Users/annaconneally/Desktop/R_datasets/remote_memory_data/all_neuron_reads.Robj"
file_path_non_neuronal_cells <- "/Users/annaconneally/Desktop/R_datasets/remote_memory_data/TRAP2_non_neuronal_cells.Robj"
load(file_path_all_neuron_reads)
load(file_path_non_neuronal_cells)

#Load the file and libraries
library(dplyr)
library(Seurat)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

# =====================================
#  Step 2: Inspect the Seurat Object
# =====================================
print(all_para_reads_filt)  # Check object structure
head(all_para_reads_filt@meta.data)  # View metadata
# Visualise RNA Counts Across Cell Types
ggplot(all_para_reads_filt@meta.data, aes(x = celltype, y = nCount_RNA)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "RNA Counts Across Non-Neuronal Cell Types",
       x = "Cell Type", 
       y = "Total RNA Count (nCount_RNA)")

#=====================================
# STEP: 3 Isolating Microglia 
#=====================================
#This creates a seurat object called microglia, filtering from object all_para_reads_filt 
microglia <- subset(all_para_reads_filt, celltype == "3_Microglia_Il1a+" | celltype == "4_Microglia_Il1a-")
#removing cells that belong to a sample (orig.ident) with fewer than 3 cells.
microglia<- subset(microglia, cells = colnames(microglia)[table(microglia$orig.ident)[microglia$orig.ident] >= 3])
head(microglia)

#normalise - Accounts for differences in sequencing depth between cells.
microglia <- NormalizeData(microglia, normalization.method = "LogNormalize", scale.factor = 10000)
microglia <- FindVariableFeatures(microglia, selection.method = "vst", nfeatures = 2000)

#scale
all.genes <- rownames(microglia)
microglia<- ScaleData(microglia, features = all.genes)
microglia[["RNA"]]$scale.data

#top10hvg
microglia_top10_hvgenes <- head(VariableFeatures(microglia), 10)
plot1<- VariableFeaturePlot(microglia)
plot2 <- LabelPoints(plot = plot1, points = microglia_top10_hvgenes, repel = TRUE)
plot1 + plot2

#Visualse QC metrics of microglia
microglia[["percent.mt"]] <- PercentageFeatureSet(microglia, pattern = "^mt-")
VlnPlot(microglia, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(microglia, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(microglia, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#PCAanalysis
microglia <- RunPCA(microglia, features = VariableFeatures(object = microglia))
print(microglia[["pca"]], dims = 1:05, nfeatures = 5)
DimPlot(microglia, reduction = "pca")

#================================================
# STEP 3 CREATE CLUSTERS & DEFINE CLUSTER MARKERS
#================================================
#findclusters
microglia <- FindNeighbors(microglia, dims = 1:10)
microglia<- FindClusters(microglia, resolution = 0.2)
microglia <- RunUMAP(microglia, dims = 1:10)
DimPlot(microglia, reduction = "umap")
DimPlot(microglia, group.by = "orig.ident")

#Finding differentially expressed features (cluster biomarkers)
cluster0.markers<-FindMarkers(microglia, ident.1=0)
cluster1.markers<-FindMarkers(microglia, ident.1=1)
cluster2.markers<-FindMarkers(microglia, ident.1=2)
cluster3.markers<-FindMarkers(microglia, ident.1=3)
head(cluster0.markers, n = 5)

#find markers for every cluster compared to all remaining cells
microglia.markers <- FindAllMarkers(microglia, min.pct = .1)
microglia.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
markers<-read.csv("filtered_microglia_markers.csv",stringsAsFactors = FALSE)     

#========================================================================
# STEP 4 : LOOPING IN T5 MICE AND GROUPING MICE BY IDENTITY AND CONDITION
#=======================================================================
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

#==============================
# STEP 5 GENE ONTOLOGY ANALYSIS 
#==============================
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)   # Use org.Hs.eg.db for human data if needed
library(ggplot2)
if (!requireNamespace("enrichplot", quietly = TRUE))
  BiocManager::install("enrichplot")
library(enrichplot)

# Step 1: Load and Filter Marker Data
# Read your CSV file of significant markers (e.g., filtered to have p_val_adj < 0.05 and avg_log2FC > 1)
markers <- read.csv("filtered_microglia_markers.csv", stringsAsFactors = FALSE) #stringsAsFactors = FALSE ensures that text columns (e.g., gene names) remain as character strings rather than being automatically converted into factors.
cat("Loaded", nrow(markers), "markers\n") 
# Print the number of markers loaded
head(markers)

# Step 2: Perform Enrichment Analysis per Cluster
clusters <- unique(significant_markers$cluster)
# Create a list to store enrichment results for each cluster
enrichment_results <- list()

for (cl in clusters) {
  # Extract unique gene symbols for the current cluster
  genes <- significant_markers %>% 
    filter(cluster == cl) %>% 
    pull(gene) %>% 
    unique()
  
  # Perform GO enrichment analysis (over-representation analysis) for Biological Process (BP)
  ego <- enrichGO(
    gene         = genes,
    OrgDb        = org.Mm.eg.db,  # Change this if using another organism
    keyType      = "SYMBOL",
    ont          = "BP",          # Choose "BP", "MF", "CC", or "ALL" as needed
    pAdjustMethod= "BH",
    qvalueCutoff = 0.05,
    readable     = TRUE
  )
  
  # Store the enrichment result in the list
  enrichment_results[[paste0("Cluster_", cl)]] <- ego
  cat("Completed enrichment analysis for Cluster", cl, "\n")
}

# Step 3: Visualize the Enrichment Results for Each Cluster with Colour Gradients
# For each cluster’s enrichment result, create a dot plot where:
# The x-axis shows the GeneRatio,
# The y-axis shows the GO Term (Description),
# - The dot size corresponds to the Count (number of genes),
# - The dot color reflects the significance (-log10(p.adjust)) with a color gradient.
plots <- list()
for (cl in names(enrichment_results)) {ego_obj <- enrichment_results[[cl]]
  
  # Convert the enrichment result to a data frame
  df_ego <- as.data.frame(ego_obj)
  
  # Skip if there are no enriched terms for this cluster
  if(nrow(df_ego) == 0) {cat("No enrichment results for", cl, "\n")next}
  
  # Select the top 10 GO terms by adjusted p-value (lowest first)
  df_top10 <- df_ego %>%
    arrange(p.adjust) %>%
    head(10) %>%
    mutate(logP = -log10(p.adjust))
  
  # Create a dot plot using ggplot2
  p <- ggplot(df_top10, aes(x = GeneRatio, y = reorder(Description, logP))) +
    geom_point(aes(size = Count, color = logP)) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = paste("GO Enrichment for", cl),
         x = "Gene Ratio",
         y = "GO Term",
         color = "-log10(p.adjust)") +
    theme_minimal()
  
  # Store and print the plot for the cluster
  plots[[cl]] <- p
  print(p)}

# Save All Enrichment Results
saveRDS(enrichment_results, file = "enrichment_results_by_cluster.rds")
cat("Enrichment results saved to enrichment_results_by_cluster.rds\n")

#=========================================
#STEP 6 : Differential Abundance analysis 
#=========================================
Load required libraries
library(speckle)   # Propeller
library(Seurat)    # Seurat functions
library(limma)     # contrasts and model matrix
library(ggplot2)   # for plotting
library(dplyr)     # data manipulation
library(tidyr)     # data reshaping
library(scales)    # for formatting

# Load your Seurat objects (adjust file paths as needed)
file_path_all_neuron_reads <- "/Users/annaconneally/Desktop/R_datasets/remote_memory_data/all_neuron_reads.Robj"
file_path_non_neuronal_cells <- "/Users/annaconneally/Desktop/R_datasets/remote_memory_data/TRAP2_non_neuronal_cells.Robj"
load(file_path_all_neuron_reads)
load(file_path_non_neuronal_cells)

# Check metadata structure of your loaded object
print(all_para_reads_filt)
head(all_para_reads_filt@meta.data)

# Step 1: Subset microglia cell types
microglia <- subset(all_para_reads_filt, celltype == "3_Microglia_Il1a+" | celltype == "4_Microglia_Il1a-")

# Filter out samples with fewer than 3 cells
microglia <- subset(microglia, cells = colnames(microglia)[table(microglia$orig.ident)[microglia$orig.ident] >= 3])
head(microglia)

# Step 2: Normalization and finding variable features
microglia <- NormalizeData(microglia, normalization.method = "LogNormalize", scale.factor = 10000)
microglia <- FindVariableFeatures(microglia, selection.method = "vst", nfeatures = 2000)

# Step 3: Scale the data
all.genes <- rownames(microglia)
microglia <- ScaleData(microglia, features = all.genes)

# Step 4: Perform PCA and visualize
microglia <- RunPCA(microglia, features = VariableFeatures(object = microglia))

# Step 5: Find nearest neighbors and clusters
microglia <- FindNeighbors(microglia, dims = 1:5)
microglia <- FindClusters(microglia, resolution = 0.2)

# Step 6: Add experimental conditions to metadata
microglia$condition <- NA  # Initialise
microglia$condition[microglia$orig.ident %in% c("FC", "T1")] <- "FC"
microglia$condition[microglia$orig.ident %in% c("Homecage", "T5")] <- "Homecage"
microglia$condition[microglia$orig.ident %in% c("Fearonly", "T4")] <- "Fearonly"
microglia$condition[microglia$orig.ident %in% c("Context", "T2", "T3")] <- "Context"

# Step 7: Prepare data for propeller analysis
meta_data <- microglia@meta.data
meta_data <- meta_data %>% dplyr::select(orig.ident, seurat_clusters, condition)
colnames(meta_data) <- c("sample_id", "cluster", "condition")
meta_data$cluster <- as.factor(meta_data$cluster)
meta_data$condition <- as.factor(meta_data$condition)

# Step 8: Perform propeller analysis
propeller_results <- propeller(clusters = meta_data$cluster, 
                               sample = meta_data$sample_id, 
                               group = meta_data$condition)

# Step 9: Visualize Proportions Across Clusters

# Prepare the data
prop_data <- data.frame(
  Cluster = rownames(propeller_results),
  Context = propeller_results$PropMean.Context,
  FC = propeller_results$PropMean.FC,
  Fearonly = propeller_results$PropMean.Fearonly,
  Homecage = propeller_results$PropMean.Homecage)

# Convert data to long format for ggplot
prop_long <- pivot_longer(prop_data, cols = -Cluster, names_to = "Condition", values_to = "Proportion")

# Ensure 'Condition' is a factor with the correct order
prop_long$Condition <- factor(prop_long$Condition, levels = c("FC", "Homecage", "Fearonly", "Context"))

# To ensure each cluster's proportions sum to 100%
prop_long <- prop_long %>%
  group_by(Cluster) %>%
  mutate(Proportion = Proportion / sum(Proportion)) %>%
  ungroup()  # Normalize the proportions within each cluster

# Add a percentage label for visualization
prop_long$Percentage <- round(prop_long$Proportion * 100, 1)

# Create the bar plot to show proportions across conditions within clusters
ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +  # Normalized to 100% per cluster
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_fill(vjust = 0.5), 
            color = "black", size = 5) +  # Add percentage labels inside bars
  theme_minimal() +
  labs(title = "Proportion of Experimental Conditions Across Clusters",
       x = "Cluster", 
       y = "Proportion", 
       fill = "Condition") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Step 10: Calculate cell counts per condition across clusters
cell_counts <- table(meta_data$cluster, meta_data$condition)

# Convert the cell counts table into a data frame
cell_counts <- as.data.frame(cell_counts)
colnames(cell_counts) <- c("Cluster", "Condition", "Freq")

# Now, the 'Cluster' and 'Condition' columns will be available for plotting
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



propeller_results <- propeller(clusters = meta_data$cluster, 
                               sample = meta_data$sample_id, 
                               group = meta_data$condition)
#Visualisation 
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
    TRUE ~ Condition  # Keep other labels unchanged
  ))
cell_summary <- mouse_cluster_counts %>%
  group_by(Cluster, Condition) %>%
  summarise(
    Mean = mean(Cell_Count),
    SEM = sd(Cell_Count) / sqrt(n()),
    .groups = "drop"
  )
ggplot(cell_summary, aes(x = Cluster, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), 
                position = position_dodge(0.7), width = 0.3, color = "black") +  
  labs(title = "Cell Counts Across Conditions with SEM",
       x = "Cluster", y = "Mean Cell Count ± SEM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Ensure correct condition mapping
cell_summary$Condition <- factor(cell_summary$Condition, levels = c("FC", "Homecage", "Fearonly", "Context"))
prop_long$Condition <- factor(prop_long$Condition, levels = c("FC", "Homecage", "Fearonly", "Context"))
# Define consistent colors
condition_colors <- c("FC" = "#D55E00", "Homecage" = "#009E73", "Fearonly" = "#56B4E9", "Context" = "#CC79A7")
ggplot(cell_summary, aes(x = Cluster, y = Mean, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), 
                position = position_dodge(0.7), width = 0.3, color = "black") +  
  scale_fill_manual(values = condition_colors) +  # Apply consistent colors
  labs(title = "Cell Counts Across Conditions with SEM",
       x = "Cluster", y = "Mean Cell Count ± SEM") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = paste0(Percentage, "%")), 
            position = position_fill(vjust = 0.5), 
            color = "black", size = 5) +
  scale_fill_manual(values = condition_colors) +  # Apply consistent colors
  labs(title = "Proportion of Experimental Conditions Across Clusters",
       x = "Cluster", y = "Proportion", fill = "Condition") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# ============================================================================================================================
# A: Side analysis of the Ila+ and Ila- labelled microglia: Identify Differentially Expressed Genes (DEGs) Between Il1a+ and Il1a-
# =============================================================================================================================
microglia_DEG_ila <- FindMarkers(all_para_reads_filt, ident.1 = "3_Microglia_Il1a+", ident.2 = "4_Microglia_Il1a-", min.pct = 0.25)
head(microglia_DEG_ila)  # View top DEGs

# ================================================================
# Step 7: GO Enrichment for All Microglia (Il1a+ & Il1a- Together)
# ================================================================
# Extract highly expressed genes - maybe use hvg for this
microglia_avg_exp <- AverageExpression(microglia_cells, assay = "RNA", return.seurat = FALSE)
# Convert Symbols to Entrez IDs (Top 500 Genes)
top_genes <- rownames(microglia_avg_exp$RNA[order(-microglia_avg_exp$RNA[, 1]), ])[1:500]
gene_entrez_all <- na.omit(mapIds(org.Mm.eg.db, keys = top_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
# Perform GO Enrichment for All Microglia
go_results_all <- enrichGO(
  gene = gene_entrez_all,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)
# Visualise GO Enrichment for All Microglia
dotplot(go_results_all, showCategory = 10) +
  ggtitle("GO Enrichment for All Microglia (Il1a+ and Il1a-)")

# =====================================
# Step 8: Separate GO Enrichment for Il1a+ and Il1a-
# =====================================
# Extract Upregulated Genes for Each Group
Il1a_plus_DEGs <- rownames(microglia_markers[microglia_markers$avg_log2FC > 0 & microglia_markers$p_val_adj < 0.05, ])
Il1a_minus_DEGs <- rownames(microglia_markers[microglia_markers$avg_log2FC < 0 & microglia_markers$p_val_adj < 0.05, ])
# Convert to Entrez IDs
Il1a_plus_entrez <- na.omit(mapIds(org.Mm.eg.db, keys = Il1a_plus_DEGs, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
Il1a_minus_entrez <- na.omit(mapIds(org.Mm.eg.db, keys = Il1a_minus_DEGs, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first"))
# Perform GO Enrichment for Il1a+
go_results_Il1a_plus <- enrichGO(
  gene = Il1a_plus_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)
# Perform GO Enrichment for Il1a-
go_results_Il1a_minus <- enrichGO(
  gene = Il1a_minus_entrez,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

# ================================================
# Step 9: Compare GO Terms Between Il1a+ and Il1a-
# ================================================
# Convert GO results to dataframes
go_plus_df <- as.data.frame(go_results_Il1a_plus)
go_minus_df <- as.data.frame(go_results_Il1a_minus)

# Extract GO Term Descriptions
go_plus_terms <- go_plus_df$Description
go_minus_terms <- go_minus_df$Description

# Identify Shared & Unique GO Terms
shared_terms <- intersect(go_plus_terms, go_minus_terms)
unique_plus <- setdiff(go_plus_terms, go_minus_terms)
unique_minus <- setdiff(go_minus_terms, go_plus_terms)

# Print Results
cat("Shared GO Terms:\n", shared_terms, "\n")
cat("\nGO Terms Unique to Il1a+:\n", unique_plus, "\n")
cat("\nGO Terms Unique to Il1a-:\n", unique_minus, "\n")

# ===================================================
# Step 10: Visualise GO Enrichment for Il1a+ vs Il1a-
# ====================================================
# Il1a+ GO Enrichment Plot
dotplot(go_results_Il1a_plus, showCategory = 10) + ggtitle("GO Enrichment for Il1a+ Microglia")

# Il1a- GO Enrichment Plot
dotplot(go_results_Il1a_minus, showCategory = 10) + ggtitle("GO Enrichment for Il1a- Microglia")

# =====================================
# Step 11: Save Results
# =====================================
write.csv(as.data.frame(go_results_all), "GO_All_Microglia.csv")
write.csv(as.data.frame(go_results_Il1a_plus), "GO_Il1a_Plus.csv")
write.csv(as.data.frame(go_results_Il1a_minus), "GO_Il1a_Minus.csv")
saveRDS(all_para_reads_filt, "processed_microglia_seurat.rds")
# Show First 10 GO Terms and Associated Genes
go_results_df <- as.data.frame(go_results_all)
head(go_results_df[, c("Description", "geneID")], 10)
print(all_para_reads_filt)

