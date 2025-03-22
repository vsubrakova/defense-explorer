library(WGCNA)
library(gplots)
library(patchwork)
library(tidyverse)

data <- read.csv("scripts/jaccard_genome_cluster.csv", sep=",")
rownames(data) <- data$X
data$X <- NULL


##########################################
###  Choose soft threshold parameters ####
##########################################


# Call the network topology analysis function
sft <- pickSoftThreshold(data,
                         powerVector = powers,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot2::ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot2::ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
a1+a2
print(a1)
print(a2)
#grid.arrange(a1, a2, nrow = 2)

power=6 # we choose power=6


threshold <- pickHardThreshold(data, RsquaredCut = 0.85)

threshold.data <- threshold$fitIndices

# visualization to pick power

a1 <- ggplot2::ggplot(threshold.data, aes(Cut, SFT.R.sq, label = Cut)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Hard threshold cut', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot2::ggplot(threshold.data, aes(Cut, mean.k., label = Cut)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Hard threshold cut', y = 'Mean Connectivity') +
  theme_classic()

a1 + a2


threshold <- 0.7 # we choose threshold = 0.7


##########################################
###  Transform similarity matrix to adjacency matrix ####
##########################################

data_power <- data ^ power
adjacency_matrix <- ifelse(data_power<threshold, 0,1)


##########################################
###  Perform hierearchical clustering ####
##########################################


cluster_tree = hclust(as.dist(1-adjacency_matrix), method = "average")

plot(cluster_tree, xlab="", sub="", main = "Protein clustering on adjacency matrix",
     labels = FALSE, hang = 0.04);
dev.off()

dynamicMods = cutreeDynamic(dendro = cluster_tree, distM = as.matrix(1 - adjacency_matrix),
                            pamRespectsDendro = TRUE,
                            minClusterSize = 30)
table(dynamicMods)
length(table(dynamicMods))

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(cluster_tree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()


TOMplot(adjacency_matrix, cluster_tree, dynamicColors, main = "clusters")

# extract protein names in clusters

data_colors <- data.frame(rownames(data), dynamicColors, stringsAsFactors = FALSE)

cluster_turquoise <- data_colors %>% filter(dynamicColors == "turquoise") %>% select(1) %>% as.vector() %>% unlist() %>% unname()

cluster_blue <- data_colors %>% filter(dynamicColors == "blue") %>% select(1) %>% as.vector() %>% unlist() %>% unname()

