library(uwot)
library(dbscan)
library(scales)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(caret)
library(tidyverse)
library(magrittr)

setwd("/Users/jalil/phd/PhD/Psilantro/Sudmex_alcohol_stress/")

Microglia_classes <- read_csv("Data/Inmunofluorescence/Morphologies_UMAP.csv")
metrics_umap <- read_table("Data/Inmunofluorescence/metrics4UMAP.txt")
Morph_metrics_CA_features <- read_csv("Data/Inmunofluorescence/Morph_metrics_CA_features.csv") %>% drop_na()
Morph_metrics_CA_features <- Morph_metrics_CA_features %>% filter(is.finite(Circularity))

n.neighbours <- 15
n.components <- 2

set.seed(420)
# Find projection using annotated cells
u1 <- umap(Microglia_classes %>% dplyr::select(metrics_umap$x), y = Microglia_classes %>% dplyr::select(metrics_umap$x), 
           scale = TRUE, ret_model = TRUE, n_neighbors = n.neighbours, n_components = n.components)
u1$layout <- Microglia_classes %>% dplyr::select(UMAP.V1,UMAP.V2) %>% as.matrix()

# Project the whole data set in the new space
u2 <- umap_transform(Morph_metrics_CA_features %>% dplyr::select(metrics_umap$x), u1)

df <- as.data.frame(u2) %>% mutate(Morphology = "Not annotated")

# Extract clusters in UMAP space
dbr <- hdbscan(df[,1:n.components], minPts = n.neighbours)
clusters <- factor(dbr$cluster)

# Confusion matrix
# annotated.idx <- which(df$Morphology != 'Not annotated')
# predictions <- clusters[annotated.idx]
# reference <- df$Morphology[annotated.idx]
# unassigned.idx <- which(predictions == 'Unassigned')
# predictions <- as.factor(predictions[-unassigned.idx])
# predictions <- droplevels(predictions)
# reference <- as.factor(reference[-unassigned.idx])
# confusion.matrix <- caret::confusionMatrix(predictions, reference, mode = "everything")

# Read back original data to add predictions
CA_phenotypes <- df %>% add_column(Predictions.from.UMAP = clusters,
                                    UMAP.V1 = df$V1, UMAP.V2 = df$V2) %>% 
      as_tibble() %>% select(-c(V1,V2))

classes <- c("Round", "Inflamed_Ameboid", "Fried_egg", "Inflamed_Fried_egg", "Hypertrophic", "Inflamed_Hypertrophic", "Bipolar")

Morph_metrics_CA_phenotype <- Morph_metrics_CA_features %>% add_column(CA_phenotypes %>% select(-Morphology))

CA_confusion_matrix <- table(Morph_metrics_CA_phenotype$Group,Morph_metrics_CA_phenotype$Predictions.from.UMAP) %>% as_tibble(.name_repair = "minimal") %>% set_colnames(c("Group","cluster","n"))

#CA_confusion_matrix %>% mutate(percentage = n/nrow(CA_phenotypes) * 100) %>% write_csv("Data/Inmunofluorescence/clusters_CA.csv")

# Plot CA matrix ---------------------------------------------------------

classes <- c("Round", "Inflamed_Ameboid", "Fried_egg", "Inflamed_Fried_egg", "Hypertrophic", "Inflamed_Hypertrophic", "Bipolar", "Not annotated")

CA_confusion_matrix %>% 
  mutate(morphology = case_when(cluster == 0 ~ classes[1],
                                cluster == 1 ~ classes[2],
                                cluster == 2 ~ classes[3],
                                cluster == 3 ~ classes[4],
                                cluster == 4 ~ classes[5],
                                cluster == 5 ~ classes[6],
                                cluster == 6 ~ classes[7],
                                cluster == 7 ~ classes[8])) %>%
  ggplot(aes(x = Group, y = morphology, fill = n)) +
  geom_tile() +
  #geom_text(aes(label = n), vjust = 1) +
  scale_fill_gradientn(colours = brewer.pal(9, "Reds")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Group", y = "Microglia subtype", fill = "Number of cells") 

