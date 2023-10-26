
###############################################
# Figure 4
###############################################


# load phasetype R
library(PhaseTypeR)
library(igraph)
c<- 0.6
N=2


Y1= matrix(
  c(0, 1, 0, c/4, c/2, (1 - c), 1/(4*N), 1/(2*N), 1 - 1/N),ncol=3,byrow=T)

Y <- DPH(Y1,init_probs =  c(1/3,1/3,1/3))

Y_network <- phase_type_to_network(Y) %>% delete_vertices(1) 




# Change node names
vertex_attr(Y_network,'name') <- c(
  "1ind","1mp",
  "2mp","coal"
)
# Plot graph


plot(Y_network, edge.curved=.2, edge.color = 'blue', asp = 3/4,
     vertex.label.color = 'black', vertex.frame.color = NA,
     edge.arrow.size=0.5, vertex.size = 20,
     vertex.label.family="Helvetica")


 
 