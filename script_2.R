#Importation des données
# --------------------------------------------------------------
expData <- read.table("cell-cycle_SCERE_DUO.txt", row.names = 1, 
                      sep = "\t", header = T)


# Représentations graphiques des données (boxplot - GeneProfiles2)
# --------------------------------------------------------------

#-- Exemple 1 (Boxplot) :
boxplot(expData)
#-- Exemple 2 (Profils d'expression des gènes - fonction de Gaëlle) :
plotGenes2 <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    yMax = max(expData)
    
  }
  
  # Representation of the first expression profile
  plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
       ylim = c(floor(yMin), ceiling(yMax)),
       xlab = "Time point", ylab = "Gene expression level",
       main = title)
  
  # Add expression profile for other genes
  for(i in 2:nrow(expData)){
    
    lines(1:ncol(expData), expData[i,], col = "grey")
    
    # end of for()  
  }
  
  
  # Average expression profile
  if(meanProfile == TRUE){
    expMean = apply(expData, 2, mean)
    lines(1:ncol(expData), expMean, col = "red", 
          lwd = 1.5, lty = "dashed")
  }
  
  # end of function plotGenes()  
} 

plotGenes2(expData) 

# --------------------------------------------------------------
# Clustering avec la méthode k-means
# --------------------------------------------------------------

N <- 5
method = "kmeans"
distance = "euclidean"

if((method == "kmeans")&(distance == "euclidean")){
  
  for(i in 1:N){
    resKmeans <- kmeans(expData, centers = N)
    cluster_km_eu <- expData[which(resKmeans$cluster == i),]
    plotGenes2(cluster_km_eu)
    heatmap(as.matrix(cluster_km_eu))
  }}else if((method == "kmeans")&(distance == "correlation")){
    
    cor(expData)
    cor(t(expData))
    matDist <- as.dist(1 - cor(t(expData)))
    
    resKmeans_Cor <- kmeans(matDist, centers = N)
    cluster_Km_Cor <- expData[which(resKmeans_Cor$cluster == i),]
    plotGenes2(cluster_Km_Cor)
    heatmap(as.matrix(cluster_Km_Cor))
  }



#nommer les heatmaps et les plots 
assign(paste("cluster", i, "Km"), cluster, envir = parent.frame(0))


# --------------------------------------------------------------
# Changer l'algorithm: Classification HCL 
# --------------------------------------------------------------
N <- 5

#distance euclidean 
dist <- dist(expData)

#distance de correlation



#Mon code ne marche pas
if((method == "HCL")&(distance == "euclidean")){
  
  resHCL <- hclust(dist)
  for(i in 1:N){
  cluster_dist_HCL <- expData[which(cutree(resHCL, k = N) == i),]
  plotGenes2(cluster_dist_HCL)
  heatmap(as.matrix(cluster_dist_HCL))}
}else if((method == "HCL")&(distance == "correlation")){
  cor(expData)
  cor(t(expData))
  matDist <- as.dist(1 - cor(t(expData)))
  
  for(i in 1:N){
    
    cluster_Cor_HCL <- expData[which(cutree(resHCL, k = N) == i),]
    plotGenes2(cluster_Cor_HCL)
    heatmap(as.matrix(cluster_Cor_HCL))
  }
  
}

