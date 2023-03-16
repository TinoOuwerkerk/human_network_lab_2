#call library
library(igraph)
library(RColorBrewer)
library(visNetwork)
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

#function to calculate centrality metrics
which.max(degree(Highschool, mode = "all"))
which.max(closeness(Highschool, normalized = TRUE))
which.max(betweenness(Highschool, directed = FALSE, normalized = TRUE))
which.max(eigen_centrality(Highschool)$vector)
#function to visualize the network (with interaction)
set.seed(100)
Highschool_interactive_layout<-visNetwork(data.frame(id=V(Highschool)$name),
                                          highschool_edges, main = "Highschool",submain="Can zoom in/out to check the IDs and
ties") %>%
  visIgraphLayout(layout = "layout_nicely",smooth = FALSE) %>%
  visNodes(shape="circle",label = TRUE) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), nodesIdSelection = T)
Highschool_interactive_layout
