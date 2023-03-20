#call library
library(igraph)
library(RColorBrewer)
library(visNetwork)
library(ggplot2)
library(tidyverse)
library(ggpubr)
setwd("C:/Users/tinot/uu/human_network_lab_2")
#datainput
highschool_edge<-read.csv("Highschool_network_edge.csv",header=FALSE)
highschool_att<-read.csv("Highschool_network_att.csv",header = TRUE)
facebook_edge<-read.csv("Facebook_network_edge.csv",header=FALSE)
facebook_att<-read.csv("Facebook_network_att.csv",header = TRUE)

#build high school network
highschool_nodes<-data.frame(name=as.character(highschool_att$NodeID),
 gender=as.character(highschool_att$Gender),
 hall=as.character(highschool_att$Hall))
highschool_edges<-data.frame(from=c(as.character(highschool_edge[,1])),
 to=c(as.character(highschool_edge[,2])))
Highschool<-graph_from_data_frame(highschool_edges,directed = FALSE,vertices =
highschool_nodes)
co <- components(Highschool)
Highschool <- induced.subgraph(Highschool, which(co$membership == which.max(co$csize)))

#use only the largest component for analysis
summary(Highschool)

#build facebook network
facebook_nodes<-data.frame(name=as.character(facebook_att$NodeID))
facebook_edges<-data.frame(from=c(as.character(facebook_edge[,1])),
 to=c(as.character(facebook_edge[,2])))
Facebook<-graph_from_data_frame(facebook_edges,directed = FALSE,vertices = facebook_nodes)
summary(Facebook)

# Question 1
#function to calculate centrality metrics
degree_n = (degree(Highschool, mode = "all"))
closeness_n = (closeness(Highschool, normalized = TRUE))
betweenness_n = betweenness(Highschool, directed = FALSE, normalized = TRUE)
eigen_n = eigen_centrality(Highschool)

which.max(degree_n)
which.max(closeness_n)
which.max(betweenness_n)
which.max(eigen_n$vector)
#function to visualize the network (with interaction)
set.seed(100)
Highschool_interactive_layout<-visNetwork(data.frame(id=V(Highschool)$name),
                                          highschool_edges, main = "Highschool",submain="Can zoom in/out to check the IDs and
ties") %>%
  visIgraphLayout(layout = "layout_nicely",smooth = FALSE) %>%
  visNodes(shape="circle",label = TRUE) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)
Highschool_interactive_layout

# Question 2
# Correlation High school network
ggplot() + geom_point(aes(degree_n, betweenness_n)) + geom_smooth(aes(degree_n, betweenness_n), method="lm") + stat_cor(aes(degree_n, betweenness_n), method = "pearson")
ggplot() + geom_point(aes(degree_n, closeness_n)) + geom_smooth(aes(degree_n, closeness_n), method="lm") + stat_cor(aes(degree_n, closeness_n), method = "pearson")
ggplot() + geom_point(aes(degree_n, eigen_n$vector)) + geom_smooth(aes(degree_n, eigen_n$vector), method="lm") + stat_cor(aes(degree_n, eigen_n$vector), method = "pearson")

# Correlation Facebook network
degree_n_fb = (degree(Facebook, mode = "all"))
closeness_n_fb = (closeness(Facebook, normalized = TRUE))
betweenness_n_fb = betweenness(Facebook, directed = FALSE, normalized = TRUE)
eigen_n_fb = eigen_centrality(Facebook)
ggplot() + geom_point(aes(degree_n_fb, betweenness_n_fb)) + geom_smooth(aes(degree_n_fb, betweenness_n_fb), method="lm") + stat_cor(aes(degree_n_fb, betweenness_n_fb), method = "pearson")
ggplot() + geom_point(aes(degree_n_fb, closeness_n_fb)) + geom_smooth(aes(degree_n_fb, closeness_n_fb), method="lm") + stat_cor(aes(degree_n_fb, closeness_n_fb), method = "pearson")
ggplot() + geom_point(aes(degree_n_fb, eigen_n_fb$vector)) + geom_smooth(aes(degree_n_fb, eigen_n_fb$vector), method="lm") + stat_cor(aes(degree_n_fb, eigen_n_fb$vector), method = "pearson")
