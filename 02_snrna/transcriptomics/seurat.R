#Basic script for Seurat Analysis Workflow 
#R 4.4.2

library(Seurat)
library (dplyr)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Run downstream analyses
integrated_object <- RunPCA(integrated_object)
ElbowPlot(integrated_object)
integrated_object <- FindNeighbors(object = integrated_object, dims = 1:20)
integrated_object <- FindClusters(object = integrated_object, resolution = 0.8)
integrated_object <- RunUMAP(object = integrated_object, reduction = "pca", dims = 1:20)
DimPlot(integrated_object, reduction = "umap")

#Differential gene expression testing
integrated_object = PrepSCTFindMarkers(integrated_object)
DE = FindMarkers(integrated_object, ident.1 = "ident.1", ident.2 = "ident.2", test.use="wilcox", min.pct = 0.1, log2FC_threshold = 0.50)

#Visualization
EnhancedVolcano(DE, lab = rownames(DE), x = 'avg_log2FC', y = 'p_val_adj', title = 'title', pCutoff = 0.05, FCcutoff = 0.50, pointSize = 2.0, labSize = 4.5, 
selectLab = c("Gene1", "Gene2",...), boxedLabels = TRUE, drawConnectors = TRUE, col=c('gray34', 'gray34', 'blue', 'red2'), colAlpha = 0.6)

#GO analysis with clusterProfiler

library(clusterProfiler)
library(org.Mm.eg.db)

#GO Biological Process Overrepresentation plots

# Extract upregulated (logFC > 0) and downregulated (logFC < 0) genes
upregulated_genes <- rownames(DE %>% filter(avg_log2FC > 0))
downregulated_genes <- rownames(DE %>% filter(avg_log2FC < 0))

library(clusterProfiler)
library(org.Mm.eg.db)
# Perform GO enrichment for upregulated genes
ego_up <- enrichGO(gene = upregulated_genes,
                   OrgDb = org.Mm.eg.db,  # Change to org.Mm.eg.db for mouse
                   keyType = "SYMBOL",
                   ont = "BP",  # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# Perform GO enrichment for downregulated genes
ego_down <- enrichGO(gene = downregulated_genes,
                     OrgDb = org.Mm.eg.db,
                     keyType = "SYMBOL",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)

# Convert results to data frames
df_up <- as.data.frame(ego_up)
df_down <- as.data.frame(ego_down)

# Remove NA values
df_up <- na.omit(df_up)
df_down <- na.omit(df_down)

# Extract top GO terms for upregulated genes
top_up <- as.data.frame(ego_up)[1:N, c("ID", "Description", "p.adjust")]

# Extract top 5 GO terms for downregulated genes
top_down <- as.data.frame(ego_down)[1:N, c("ID", "Description", "p.adjust")]


# Ensure ego_up and ego_down are in data frame format
top_up <- as.data.frame(ego_up)[1:x, c("Description", "p.adjust")]
top_down <- as.data.frame(ego_down)[1:N, c("Description", "p.adjust")]


# Add regulation type
top_up$Regulation <- "Upregulated"
top_down$Regulation <- "Downregulated"

# Convert p.adjust to -log10(p.adjust) for better visualization
top_up$logP <- -log10(top_up$p.adjust)
top_down$logP <- log10(top_down$p.adjust)  # Negative for left-side bars

# Combine into one data frame
top_go <- rbind(top_up, top_down)

# Remove rows with NA values in p.adjust
top_go <- na.omit(top_go)

top_go <- top_go[!duplicated(top_go$Description), ]

# Reorder Description: Down on top, Up on bottom
top_go <- top_go %>%
  arrange(Regulation, logP) %>%
  mutate(Description = factor(Description, levels = rev(Description)))

# Convert Description to factor for ordered plotting
top_go$Description <- factor(top_go$Description, levels = rev(unique(top_go$Description)))

ggplot(top_go, aes(x = Description, y = logP, fill = Regulation)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for horizontal bars
  theme_minimal() +
  labs(title = "title",
       x = "GO Term", y = "-log10 Adjusted P-value") +
  scale_fill_manual(values = c("Upregulated" = "salmon", "Downregulated" = "lightblue3")) +
  theme(text = element_text(size = 12), axis.text.y = element_text(size = 12)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")  # Add center line

#Heatmaps and DotPlots

avg_exp <- AggregateExpression(object, 
                             assay = "SCT", 
                             group.by = "grouping variable", 
                             slot = "data")  # use log-normalized data

avg_mat <- avg_exp$SCT


library(pheatmap)
pheatmap(avg_mat, 
         scale = "row",         # standardize gene rows (Z-score)
         cluster_cols = TRUE or FALSE, 
         cluster_rows = TRUE,
         show_rownames = TRUE or FALSE,
         main = "title",
         fontsize_row = 10, 
         fontsize_col = 10)

DotPlot(object, features = "selected features", dot.scale = 12, cluster.idents = FALSE) +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  theme(axis.text.y.left = element_text(size = 10)) +
  ggtitle("title") +
  scale_color_gradient(low = "blue", high = "red2")

#Cell-Cell Communication Analysis 

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

# Use appropriate database for your species
CellChatDB <- CellChatDB.mouse  # or CellChatDB.human

# Subset Seurat object for CS-injury vs CTX-control
injury <- subset(object, subset = condition == "injury_x")
control <- subset(object, subset = condition == "control_x")

# Run CellChat for CS-injury
cellchat_injury <- createCellChat(object = injury, group.by = "cluster_annotation")
cellchat_injury@DB <- CellChatDB.mouse
cellchat_injury <- subsetData(cellchat_injury)
cellchat_injury <- identifyOverExpressedGenes(cellchat_injury)
cellchat_injury <- identifyOverExpressedInteractions(cellchat_injury)
cellchat_injury <- computeCommunProb(cellchat_injury)
cellchat_injury <- filterCommunication(cellchat_injury, min.cells = 10)
cellchat_injury <- computeCommunProbPathway(cellchat_injury)
cellchat_injury <- aggregateNet(cellchat_injury)

# Run CellChat for CTX-control
cellchat_control <- createCellChat(object = control, group.by = "cluster_annotation")
cellchat_control@DB <- CellChatDB.mouse
cellchat_control <- subsetData(cellchat_control)
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
cellchat_control <- computeCommunProb(cellchat_control) 
cellchat_control <- filterCommunication(cellchat_control, min.cells = 10)
cellchat_control <- computeCommunProbPathway(cellchat_control)
cellchat_control <- aggregateNet(cellchat_control)

# Merge CellChat objects
cellchat_merged <- mergeCellChat(list(cellchat_injury, cellchat_control), add.names = c("X" "x"))


# Visualize
compareInteractions(cellchat_merged, show.legend = F, group = c(1, 2)) #Bar plot showing total number or strength of interactions across groups.
compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")


netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE) #Circle plot showing gain/loss in communication between conditions.
netVisual_diffInteraction(cellchat_merged,
                          weight.scale = T,
                          measure = "weight")  # or "count"


rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE) #Ranks strongest sender/receiver roles in the network.
rankNet(cellchat_merged, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)

# Extract the ranked pathways
ranked_pathways <- rankNet(cellchat_merged, mode = "comparison")

head(ranked_pathways)

# Subset for a top pathway (e.g., "XX")
df_XX <- subsetCommunication(cellchat_injury, signaling = "XX")

library(dplyr)

# Rank sender (source) cell types
sender_rank <- df_XX %>%
  group_by(source) %>%
  summarize(total_prob = sum(prob)) %>%
  arrange(desc(total_prob))

# Rank receiver (target) cell types
receiver_rank <- df_XX %>%
  group_by(target) %>%
  summarize(total_prob = sum(prob)) %>%
  arrange(desc(total_prob))

# Combine and tag the data
sender_rank <- df_XX %>%
  group_by(cell_type = source) %>%
  summarize(total_prob = sum(prob)) %>%
  mutate(role = "Sender")

receiver_rank <- df_XX %>%
  group_by(cell_type = target) %>%
  summarize(total_prob = sum(prob)) %>%
  mutate(total_prob = -total_prob, role = "Receiver")  # Make negative for mirroring

# Combine sender and receiver into one dataframe
combined_df <- bind_rows(sender_rank, receiver_rank)

# Order cell types by sender strength to keep consistent bars
combined_df$cell_type <- factor(combined_df$cell_type, 
                                levels = unique(c(
                                  sender_rank$cell_type[order(sender_rank$total_prob)],
                                  receiver_rank$cell_type[order(receiver_rank$total_prob)])
                                ))

# Plot
ggplot(combined_df, aes(x = cell_type, y = total_prob, fill = role)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_y_continuous(labels = abs, name = "Total Communication Probability") +
  scale_fill_manual(values = c("Sender" = "orange", "Receiver" = "lightblue2")) +
  labs(title = "YYY") +
  theme_minimal()

# cell pair interactions for a specific signaling pathway
df_XX <- subsetCommunication(cellchat_injury, signaling = "XX")

# Find the top sender-receiver pairs
top_XX_pairs <- df_XX[order(-df_XX$prob), ]

# View top interactions
head(top_XX_pairs, 5)

# Subset top interactions
top_XX_pairs <- subsetCommunication(cellchat_injury, signaling = "XX") %>%
  arrange(desc(prob)) %>%
  slice(1:5)

# Create a new column with the L-R pair (ligand + receptor)
top_XX_pairs$LR_pair <- paste(top_XX_pairs$ligand, top_XX_pairs$receptor, sep = "_")

# Plot as a dot plot with annotations
ggplot(top_XX_pairs, aes(x = target, y = source, size = prob, color = prob)) +
  geom_point() +
  geom_text(aes(label = paste0("LR: ", LR_pair, "\n", "Prob: ", round(prob, 2))),
            position = position_jitter(width = 0.1, height = 0.4), size = 3, hjust = 1.5, vjust = 0.5) +
  scale_color_gradient(low = "blue", high = "red2") +
  labs(title = "Top XX Interactions",
       x = "Receiver Cell Type", y = "Sender Cell Type", size = "Probability", color = "Probability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
