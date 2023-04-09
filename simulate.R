library(igraph)
library(ggplot2)
library(reshape2)

source("networked-sir.R")
highschool_edge<-read.csv("Highschool_network_edge.csv",header=FALSE)
highschool_att<-read.csv("Highschool_network_att.csv",header = TRUE)
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

# Population size
N <- vcount(Highschool)

repeat {
  # Degree sequence drawn from geometric distribution
  # where all vertices have at least one edge
  dseq <- rgeom(n= N, prob = 0.2) + 1
  if(is_graphical(dseq)) {
    break
  }
}

# Parameters
X0 <- rep.int(0, N)
seed <- 5 # patient zero
beta <- 0.15 # infection rate
gamma <- 1
tmax <- 28
X0[seed]<-1
sim <- networked_sir(X0 = X0, G = Highschool, beta = beta, gamma = gamma, tmax = tmax)
dfm <- melt(sim$df, id.vars = "t")
ggplot(data = dfm, mapping = aes(t, value, color=variable)) +
  geom_line() +
  theme_light()
