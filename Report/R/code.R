Df <- read.table("~/Downloads/hgdp.txt")
colnames(Df)[1:3] <- c("Individual_ID","location","continent") genotypes <- Df %>% select(starts_with("V"))
convert_genotypes <- function(genotype_data) {
  freq_table <- table( c(substr(genotype_data,1,1),substr(genotype_data,2,2)))
  most_freq <- names(freq_table)[which.max(freq_table)]
  ref_char <- most_freq
  
  encoded_vector<- ifelse(substr(genotype_data,1,1) == ref_char, 1,0) + ifelse(substr(genotype_data,2,2) == ref_char,1,0)
}

encoded_df <- apply(genotypes, 2, convert_genotypes ) # Perform PCA
pca <- prcomp(encoded_df, scale = TRUE)
#Plot the scree plot to visualize the amount of variation explained by each component.
fviz_eig(pca, xlab = "Principal Component" ,addlabels = TRUE, ylim = c(0, 10))

#Choose the K
pca_scores <- pca$x[, 1:2]
set.seed(123)
wcss = matrix(cbind(c(1:15),rep(0,15)),ncol = 2) for (i in 1:15){
  wcss[i,2] <- kmeans(pca_scores, i, 10)$tot.withinss
}
df_wcss <- as.data.frame(wcss)
ggplot(df_wcss, aes(x = V1, y = V2)) + geom_line() + geom_point() +
  labs(x = "K", y = "WCSS") +
  ggtitle("WCSS vs. K")

# Plot the first two principal components
plot(pca$x[,1], pca$x[,2], pch = 20, col = "blue", main = "PCA plot")
# Plot the PCA plot with 3-means cluster assignments
# Perform K-means clustering on the PCA scores
kmeans_res <- kmeans(pca$x[,1:2], centers = 3)
plot(pca$x[,1], pca$x[,2], pch = 20, col = kmeans_res$cluster, main = "PCA plot with 3-means clusters")

color <- as.numeric(factor(Df$continent))
colors <- c("red", "blue", "green", "orange", "purple", "black","brown")
color_vector <- colors[color]
color_table <- matrix(c(names(table(Df$continent)),colors),ncol = 2)


# Plot the first two principal components with continental information
plot(pca$x[,1], pca$x[,2], pch = 20, col = color_vector, main = "PCA plot with continental information")

# Perform 7-means clustering on the PCA scores
kmeans_res2 <- kmeans(pca$x[,1:2], centers = 7)
# Plot the PCA plot with 7-means cluster assignments
plot(pca$x[,1], pca$x[,2], pch = 20, col = kmeans_res2$cluster, main = "PCA plot with 7-means clusters")


# Perform 6-means clustering on the PCA scores
kmeans_res3 <- kmeans(pca$x[,1:2], centers = 6)

# Plot the PCA plot with 6-means cluster assignments
plot(pca$x[,1], pca$x[,2], pch = 20, col = kmeans_res3$cluster, main = "PCA plot with 6-means clusters")


set.seed(1234)
tsne_result <- Rtsne(genotypes, dims=2,verbose= FALSE)
plot(tsne_result$Y, col=color_vector, pch=19, main="t-SNE Result")

