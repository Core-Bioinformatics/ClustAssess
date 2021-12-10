## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(1001)

## ----setup, fig.width=7, fig.height=4.5---------------------------------------
library(ClustAssess)
library(e1071)
library(dbscan)
library(ggplot2)
suppressPackageStartupMessages(library(dendextend))

# we will use the Iris dataset for this vignette
data = as.matrix(iris[,1:4])
df.iris = as.data.frame(prcomp(data)$x)
df.iris$species = iris$Species
ggplot(df.iris, aes(x=PC1, y=PC2, color=species)) + 
  geom_point() + 
  labs(title='Iris PCA')

## ----clustering, fig.width=7, fig.height=4.5----------------------------------
# a flat, disjoint clustering with DBscan
db.res = dbscan(data, eps=1)$cluster
df.iris$db.cluster = as.factor(db.res)
ggplot(df.iris, aes(x=PC1, y=PC2, color=db.cluster)) + 
  geom_point() + 
  labs(title='DBScan clustering')

# an overlapping clustering with soft k-means
cmeans.res = cmeans(data, centers=6)$membership
# get the strongest cluster assignment for each observation to plot
# but we will still use the soft cluster assignments to calculate ECS
df.iris$cmeans.cluster = as.factor(apply(cmeans.res, 1, which.max))
ggplot(df.iris, aes(x=PC1, y=PC2, color=cmeans.cluster)) + 
  geom_point() + 
  labs(title='c-means clustering')

# a hierarchical clustering
distances = dist(data, method='euclidean')
hc.res = hclust(distances, method='complete')
# plot the resulting dendrogram
node.colors = c('blue', 'red', 'green')
hc.res %>% as.dendrogram() %>% 
  set("leaves_pch", 19) %>% 
  set("leaves_cex", 0.5) %>% 
  set("leaves_col", node.colors[df.iris$species]) %>% 
  plot(main='Complete linkage hierarchical clustering', leaflab='none')

## ----ecs, fig.width=7, fig.height=4.5-----------------------------------------
# which observations are clustered more similarly?
# first compare flat disjoint and flat soft clusterings
df.iris$dbscan.cmeans.ecs = element_sim_elscore(db.res, 
                                                cmeans.res)
ggplot(df.iris, aes(x=PC1, y=PC2, color=dbscan.cmeans.ecs)) + 
  geom_point() + 
  labs(title='DBScan vs c-means similarity')
mean(df.iris$dbscan.cmeans.ecs)

# next compare flat disjoint and hierarchical disjoint clusterings
df.iris$dbscan.hc.ecs = element_sim_elscore(db.res, 
                                    hc.res)
ggplot(df.iris, aes(x=PC1, y=PC2, color=dbscan.hc.ecs)) + 
  geom_point() + 
  labs(title='DBScan vs hclust similarity')
mean(df.iris$dbscan.hc.ecs)

# finally compare flat soft and hierarchical disjoint clusterings
df.iris$cmeans.hc.ecs = element_sim_elscore(cmeans.res, 
                                            hc.res)
ggplot(df.iris, aes(x=PC1, y=PC2, color=cmeans.hc.ecs)) + 
  geom_point() + 
  labs(title='c-means vs hclust similarity')
mean(df.iris$cmeans.hc.ecs)

## ----sessinf------------------------------------------------------------------
sessionInfo()

