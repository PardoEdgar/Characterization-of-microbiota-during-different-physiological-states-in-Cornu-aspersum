library(tidyverse)
library(vegan)
library(ggrepel)
library(igraph)
library(viridis) 
path <- ("/Data_Bacteria_Reads.csv") 
Data_bacteria_reads <- read_csv(path)
Phylum_abundances <- Data_bacteria_reads %>%
  group_by(Phylum) %>%
  summarise(
    Total_abundance = sum(Max, na.rm = TRUE),
    Active = sum(Active, na.rm = TRUE),
    Post_aestivation = sum(Post_aestivation, na.rm = TRUE),
    Aestivation = sum(Aestivation, na.rm = TRUE)
  ) %>%
  filter(!is.na(Phylum)) %>%  
  arrange(desc(Active))  
Family_abundances <- Data_bacteria_reads %>%
  group_by(Family) %>%
  summarise(
    Total_abundance = sum(Max, na.rm = TRUE),
    Active = sum(Active, na.rm = TRUE),
    Post_aestivation = sum(Post_aestivation, na.rm = TRUE),
    Aestivation = sum(Aestivation, na.rm = TRUE)
  ) %>%
  filter(!is.na(Family)) %>%  # Eliminar filas donde Family es NA
  arrange(desc(Active))  # Ordenar por la Abundance activa
Top_families <- Family_abundances %>%
  arrange(desc(Total_abundance)) %>%
  slice_head(n = 15)
Top_long_familes <- pivot_longer(Top_families, 
             cols = c("Active", "Post_aestivation", "Aestivation"),
             names_to = "Group", 
             values_to = "Abundance")
Top_long_familes <- Top_long_familes %>%
  group_by(Group) %>%
  mutate(Relative_abundance = Abundance / sum(Abundance))
Top_phylum <- Phylum_abundances %>%
  arrange(desc(Total_abundance)) %>%
  slice_head(n = 10)
Top_long_phylum <- pivot_longer(Top_phylum, 
             cols = c("Active", "Post_aestivation", "Aestivation"),
             names_to = "Group", 
             values_to = "Abundance")
Top_long_phylum <- Top_long_phylum %>%
  group_by(Group) %>%
  mutate(Relative_abundance = Abundance / sum(Abundance))
Abundance_long_families <- pivot_longer(Family_abundances, 
                                        cols = c("Active", "Post_aestivation", "Aestivation"),
                                        names_to = "Group", 
                                        values_to = "Abundance")
Abundance_long_families <- Abundance_long_families %>%
  group_by(Group) %>%
  mutate(Relative_abundance = Abundance / sum(Abundance))
Abundance_long_families_top <- pivot_longer(Top_families, 
                                        cols = c("Active", "Post_aestivation", "Aestivation"),
                                        names_to = "Group", 
                                        values_to = "Abundance")
Abundance_long_families_top <- Abundance_long_families_top %>%
  group_by(Group) %>%
  mutate(Relative_abundance = Abundance / sum(Abundance))
Phylum_abundances_long <- pivot_longer(Phylum_abundances, 
                                        cols = c("Active", "Post_aestivation", "Aestivation"),
                                        names_to = "Group", 
                                        values_to = "Abundance")
Phylum_abundances_long <- Phylum_abundances_long %>%
  group_by(Group) %>%
  mutate(Relative_abundance = Abundance / sum(Abundance))
windows()
ggplot(Phylum_abundances_long, aes(x = Group, y = Abundance, fill= Group)) +
  geom_bar(stat = "identity") +
  labs(x = "Phylum", y = "Total abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15)) +
  theme_classic()
ggplot(Top_families, aes(x = reorder(Family, -Total_abundance), y = Total_abundance, fill = Family)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundance Total por Family", x = "Family", y = "Abundance Total") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_classic()
ggplot(Top_phylum, aes(x = reorder(Phylum, -Total_abundance), y = Total_abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(title = "Abundance Total by phylum", x = "Phylum", y = "Abundance Total") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_classic()
ggplot(data = Phylum_abundances_long, aes(x = Group, y = Phylum, fill = Relative_abundance)) +
  geom_raster() + 
  scale_fill_viridis(option = "cividis",name = "Abundance Relativa") +
  labs(title = "Abundance Relativa por Phylum y Tratamiento",
       x = "Tratamiento",
       y = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))
ggplot(data = Abundance_long_families_top, aes(x = Group, y = Family, fill = Relative_abundance)) +
  geom_raster() + 
  scale_fill_viridis(option = "cividis",name = "Abundance Relativa") +
  labs(title = "Abundance Relativa por Family y Tratamiento",
       x = "Tratamiento",
       y = "Family") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), axis.text.y = element_text(size = 15))
windows()
ggplot(data = Top_long_familes, aes(x = Family, y = Relative_abundance, fill=Group)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent) + 
  labs(title = "Abundance Relativa por Family y Tratamiento",
       x = "Tratamiento",
       y = "Abundance Relativa (%)",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13))
ggplot(data = Top_long_phylum, aes(x = Phylum, y = Relative_abundance, fill=Group)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent) + 
  labs(title = "Abundance Relativa por Family y Tratamiento",
       x = "Tratamiento",
       y = "Abundance Relativa (%)",
       fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 13))

#Alpha diversity
Abundance_matrix <- as.matrix(Family_abundances[,3:5])
Abundance_matrix <- t(as.matrix(Abundance_matrix))
shannon_index <- vegan::diversity(Abundance_matrix, index = "shannon")
shannon_index
simpson_index <- vegan::diversity(Abundance_matrix, index = "simpson")
simpson_index
Richness <- specnumber(Abundance_matrix)
Richness

#Beta diversity
Abundance_matriz <- as.matrix(Family_abundances[,3:5])
Abundance_matriz
colnames(Abundance_matriz) <- c("Active", "Post_aestivation", "Aestivation")
# PCA
pca_result <- prcomp(Abundance_matriz, scale = TRUE)
summary(pca_result)
var_exp <- summary(pca_result)$importance[2, ] * 100
x_label <- paste0("PC1 (", round(var_exp[1], 1), "% of Variance Explained)")
y_label <- paste0("PC2 (", round(var_exp[2], 1), "% of Variance Explained)")

biplot(pca_result, scale = 0, xlab = x_label, ylab = y_label)
pca_df <- as.data.frame(pca_result$x)
pca_df$Group <- Family_abundances$Family

ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 1) +
  geom_text_repel(aes(label = Family_abundances$Family), size = 3, max.overlaps = 100) +
  theme_minimal() +
  labs(x = x_label, y = y_label)

#Microbial Network
Abundance_resumida <- Abundance_long_families %>%
  group_by(Group, Family) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE)) %>%
  ungroup()

Abundance_long_families_t <- Abundance_resumida %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = list(Total_abundance = 0))
Abundance_long_families_t_matrix <- as.matrix(Abundance_long_families_t)
cor_matrix <- cor(Abundance_long_families_t[, -1], method = "spearman")
cor_matrix_filtered <- cor_matrix
cor_matrix_filtered[abs(cor_matrix) < 0.7] <- 0
cor_matrix[is.nan(cor_matrix)] <- 0
cor_matrix[is.na(cor_matrix)] <- 0
cor_matrix[cor_matrix < 0.7] <- 0
graph <- graph.adjacency(cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
plot(graph, vertex.size = 10, vertex.label.cex = 0.8, edge.width = E(graph)$weight*2)
V(graph)$color <- ifelse(V(graph)$name %in% c("Active"), "blue",
                         ifelse(V(graph)$name %in% c("Aestivation"), "red", "green"))
plot(graph, vertex.size = 10, vertex.label.cex = 0.8, edge.width = E(graph)$weight * 2, vertex.color = V(graph)$color)
Groups_Familys <- Abundance_long_families %>%
  group_by(Family) %>%
  summarise(Group_predominante = Group[which.max(Abundance)]) %>%
  ungroup()
colores_Groups <- c("Active" = "orange", "Aestivation" = "red", "Post_aestivation" = "darkgreen")
Groups_Familys$color <- colores_Groups[Groups_Familys$Group_predominante]
V(graph)$name <- trimws(tolower(V(graph)$name))
Groups_Familys$Family <- trimws(tolower(Groups_Familys$Family))
windows()
V(graph)$color <- Groups_Familys$color[match(V(graph)$name, Groups_Familys$Family)]

plot(graph, vertex.size = 8, vertex.label.cex =0.5, edge.width = E(graph)$weight * 2, vertex.color = V(graph)$color, vertex.label.color = "black")
tkplot(graph)
