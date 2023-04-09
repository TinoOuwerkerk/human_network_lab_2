library(igraph)

networked_sir <- function(X0, G, beta, gamma, tmax) {
	# Simulate the networked SIR model with arbitrary network G (igraph)
    # using the Gillespie stochastic simulation algorithm (SSA).
    # Parameters:
    # X0 : Initial configuration of state vector X0 in {0,1,2}^N,
    #   where {0,1,2} ~ {S,I,R} and N is the number of individuals
    # G : Network (igraph object) consisting of N nodes (individuals)
    # beta : per-edge infection rate
    # gamma : Recovery rate
    # tmax : Maximum time duration of simulation
    # Returns:
    # List with data frame [t,s(t),i(t),r(t)], network G,
    # initial state vector X0, final state vector X, and parameters beta, gamma

    N <- length(X0)  # Number of individuals
    if(N != vcount(G)) {
      error("Length of X0 has to equal number of vertices!")
    }

    lambda <- numeric(length = N)  # Allocate rate vector
    num_inf_neighbours <- numeric(length = N)  # Allocate number of infected neighbors vector

    # Set state X to initial state and time t to zero
    X <- X0;
    t <- 0

    # Keep track of number of individuals in each compartment.
    st <- length(which(X == 0))
    it <- length(which(X == 1))
    rt <- length(which(X == 2))
    new_inf <- c()
    # Initialize data matrix with observations
    df <- data.frame(t = 0, s = st, i = it, r = rt)
    dm <- as.matrix(df)  # Matrix is faster than data frame to update

    # Get adjacency list of G. This appears to be the fastest way
    # to find the neighbors of a particular individual
    adj_list <- igraph::as_adj_list(G)
    days_infected <- X0
    # Compute rates for all individuals
    for(i in 1:N) {
      if(X[i] == 0) {  # If susceptible
        neighbour_indices <- adj_list[[i]]
        num_inf_neighbours[i] <- length(which(X[neighbour_indices] == 1))
        lambda[i] <- beta#*num_inf_neighbours[i]
      } else if(X[i] == 1) {  # If infected
        lambda[i] <- gamma
      } else {
        lambda[i] <- 0
      }
    }
    
    while(t < tmax) {
      lambda0 <- sum(lambda)
      newly_infected <- 0
      t <- t + 1
      pos_j <- which(X==1)
      for(infect in pos_j){
        neighbour_indices <- adj_list[[infect]]  # Get neighbor indices of individual j
        for(j in neighbour_indices){ # get neighbours
          if(X[j]==0){ # check if they are susceptible
            if(lambda[j]>runif(1)){ # check if they get infected
              X[j]<-1 #infection event
              days_infected[j] <- 1 # count days of being infected
              lambda[j] <- gamma  # Update rate for individual j
              
              # Increment number of infected neighbors count for neighborhood
              num_inf_neighbours[neighbour_indices] <- num_inf_neighbours[neighbour_indices] + 1
              
              # Bookkeeping
              st <- st - 1
              it <- it + 1
              newly_infected <- newly_infected+1
            } 
          }
          }
        if(X[infect]==1) { # Recovery event
          days_infected[infect] <- days_infected[infect]+1
          if(days_infected[infect]>3){# check if  person is 3 days infected
            X[infect] <- 2  # Set individual j to recovered
            lambda[infect] <- 0  # Set individual j rate to zero
            
            # Decrement number of infected neighbors count
            num_inf_neighbours[neighbour_indices] <- num_inf_neighbours[neighbour_indices] - 1
            
            # Bookkeeping
            it <- it - 1
            rt <- rt + 1
          }
        }
        }
      # Update rates for susceptible neighbours of individual j
      update_indices <- neighbour_indices[X[neighbour_indices] == 0]
      lambda[update_indices] <- beta*num_inf_neighbours[update_indices]

      # Update data matrix
      dm <- rbind(dm,c(t,st,it,rt))
      # print(paste0('day:',t))
      # print(paste0('newly_infected:',newly_infected))
      # print(paste0('infected:',length(which(X==1))))
      # print(paste0('recovered:',length(which(X==2))))
      new_inf <- append(new_inf, newly_infected)
    }  # end while
    # Return fractions instead of counts
    dm[,c("s","i","r")] <- dm[,c("s","i","r")]/N
    df <- as.data.frame(dm)
    row.names(df) <- NULL  # Remove row names
    return(new_inf)
  }
