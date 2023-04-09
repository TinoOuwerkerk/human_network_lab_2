import matplotlib.pyplot as plt
from random import uniform, seed
import numpy as np
import pandas as pd
import csv
import time
from igraph import *

def compute_mean_by_index(list_of_lists):
    # Get the number of inner lists
    num_lists = len(list_of_lists)
    # Get the length of each inner list
    list_length = len(list_of_lists[0])
    # Initialize a list to store the means
    means = []

    # Iterate through the indices
    for i in range(list_length):
        # Initialize a sum variable for the current index
        sum_for_index = 0
        # Iterate through the inner lists
        for j in range(num_lists):
            # Add the value at the current index to the sum
            sum_for_index += list_of_lists[j][i]
        # Compute the mean for the current index
        mean_for_index = sum_for_index / num_lists
        # Append the mean to the result list
        means.append(mean_for_index)

    return means

def IC(g,S,p=0.5,mc=1000, timestamps = 28, Monte_Carlo = True):
    """
    Input:  
    g - graph object, 
    S - set of seed nodes(dtype list)
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps
    Monte_Carlo - Boolean, determines whether to make to random seed or not

    Output: 
    - average number of nodes activated in each Monte-Carlo simulation
    - average number of nodes influenced by the seed nodes in each timestamp
    """
    
    # Loop over the Monte-Carlo Simulations
    spread, sum_spread = [], []
    for i in range(mc):
        
        # Simulate propagation process      
        Active = S[:]
        sum = 0
        for _ in range(timestamps):
            new_ones = []
            # For each active node, find its neighbors that become activated
            for node in Active:
                if Monte_Carlo == True: 
                    np.random.seed(i)
                # Determine neighbors that become infected
                success = np.random.uniform(0,1,len(g.neighbors(node,mode="out"))) < p
                new_ones += list(np.extract(success, g.neighbors(node,mode="out")))

            new_active = list(set(new_ones) - set(Active))
            sum += len(new_active) 
            # in case the network is fully activated
            if new_active == []:
                break
            # Add newly activated nodes to the set of activated nodes
            Active = list(set(Active).union(set(new_active)))
            
        spread.append(len(Active))
        sum_spread.append(sum)
        
    return(np.mean(spread), np.mean(sum_spread))



def IC_immunized(g, S, immunized, p=0.15, mc=1000, timestamps = 28, full = False, Monte_Carlo = True):
    """
    Input:  
    g - graph object, 
    S - set of seed nodes(dtype list)
    immunized - set of immunized nodes(dtype list)
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps
    Monte_Carlo - Boolean, determines whether to make to random seed or not

    Output: 
    - average number of nodes activated in each Monte-Carlo simulation
    - average number of nodes influenced by the seed nodes in each timestamp
    """
    
    # Loop over the Monte-Carlo Simulations
    spread, sum_spread, infected_a_day = [], [], []
    for i in range(mc):
        S = list(set(S).difference(set(immunized)))
        # Simulate propagation process
        Active = S[:]
        sum = 0
        infected_per_day = []
        for _ in range(timestamps):
            new_ones = []
            # For each active node, find its neighbors that become activated
            for node in Active:
                if Monte_Carlo == True: 
                    np.random.seed(i)
                # Determine neighbors that become infected
                success = np.random.uniform(0,1, len(g.neighbors(node,mode="out"))) <= p
                new_ones += list(np.extract(success, g.neighbors(node,mode="out")))
                new_ones = list(set(new_ones).difference(set(immunized)))
                
            
            new_active = list(set(new_ones).difference(set(Active)))
            sum += len(Active)
            if new_active == [] and full == True:
                break
            infected_per_day.append(len(new_active))

            # Add newly activated nodes to the set of activated nodes
            Active = list(set(Active).union(set(new_active)))
            
        spread.append(len(Active))
        sum_spread.append(sum)
        infected_a_day.append(infected_per_day)
        
    return(np.mean(spread), np.mean(sum_spread), compute_mean_by_index(infected_a_day))

def greedy(g, k, p=0.15, mc=1000, timestamps = 28):
    """
    Input:  
    g - graph object
    k - number of seed nodes
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps
    Monte_Carlo - Boolean, determines whether to make to random seed or not

    Output: optimal seed set, resulting spread, time for each iteration
    """

    S, spread, timelapse, start_time = [], [], [], time.time()
    # setting list of nodes and making the node 107 the first node in odrer to compare with other methods
    l = [i for i in range(g.vcount())]
    l.remove(106)
    l.insert(0, 106)
    
    # Find k nodes with largest marginal gain
    for _ in range(k):

        # Loop over nodes that are not yet in seed set to find biggest marginal gain
        best_spread = 0
        #for j in set(range(g.vcount())) - set(S):
        for j in set(l) - set(S):

            # Get the spread
            s = IC(g, S + [j], p, mc, timestamps, Monte_Carlo = True)

            # Update the winning node and spread so far
            if s[1] > best_spread:
                best_spread, node = s[1], j

        # Add the selected node to the seed set
        S.append(node)
        
        # Add estimated spread and elapsed time
        spread.append([best_spread, s[1]])
        timelapse.append(time.time() - start_time)

    return(S,spread, timelapse)

def greedy_immunized(g,k,p=0.1,mc=1000):

    """
    g - graph object
    k - number of seed nodes
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps
    Monte_Carlo - Boolean, determines whether to make to random seed or not

    Output: optimal seed set, resulting spread, time for each iteration
    """

    immunaized_list, spread, timelapse, start_time = [],[], [], time.time()
    
    # Find k nodes with largest marginal gain
    for _ in range(k):

        # Loop over nodes that are not yet in seed set to find biggest marginal gain
        best_spread = IC_immunized(g, [106],[], p, mc, timestamps = 28, full = True, Monte_Carlo = True)[1]
        for j in set(range(g.vcount())) - set(immunaized_list):

            # Get the spread
            s = IC_immunized(g, [106], immunaized_list + [j], p, mc, timestamps=28, full = True, Monte_Carlo = True)

            # Update the winning node and spread so far
            if s[1] < best_spread and s[1] !=0 :
                best_spread, node = s[1], j

        # Add the selected node to the seed set
        if 'node' in locals():
            immunaized_list.append(node)
        
        # Add estimated spread and elapsed time
        spread.append(best_spread)
        timelapse.append(time.time() - start_time)

    return(immunaized_list, spread, timelapse)

def threshold_model(network, node_seed, threshold, n_day, mc=100, Monte_Carlo=False):
    spread, sum_spread, infected_per_day, infected_per_iteration = [], [], [], []
    for i in range(mc):
        if Monte_Carlo == True:
            np.random.seed(i)

        nNode = network.vcount()
        node_status = np.zeros(nNode, dtype=int) # start from a healthy population
        adj_matrix = np.array(network.get_adjacency().data)
        each_neighbors = {node: np.where(adj_matrix[node, :] > 0)[0] for node in range(nNode)} # get the neighbor list of each node
        infected_a_day,infected_List = [] , []
        

        for seed in node_seed:
            node_status[seed] = 1 # adopt (value=1), don't adopt (value=0)
        sum_of_ifected = 0
        for day in range(1, 29):  #n_day + 1):
            n_infected = 0

            for node in range(nNode):
                if node_status[node] == 0:
                    neighbours = each_neighbors[node]
                    n_neighbors = len(neighbours)
                    n_adopters = np.sum(node_status[neighbours] == 1)

                    #if n_adopters/n_neighbors > threshold[node]:
                    if n_adopters > threshold[node]:
                        node_status[node] = 1
                        n_infected += 1
            # if n_infected == 0:
            #     break
            infected_at_all = np.sum(node_status == 1)
            sum_of_ifected += np.sum(node_status == 1)
            infected_a_day.append(n_infected)
            infected_List.append(infected_at_all)

        # Collect all of the necessary information
        spread_per_simulation = np.sum(node_status == 1)
        spread.append(spread_per_simulation)
        sum_spread.append(sum_of_ifected)
        infected_per_day.append(infected_a_day)
        infected_per_iteration.append(infected_List)
    return np.mean(spread), np.mean(sum_spread), compute_mean_by_index(infected_per_day), compute_mean_by_index(infected_per_iteration)

def greedy_Th(g, k, threshold, mc=1000, timestamps = 28):
    """
    Input:  
    g - graph object
    k - number of seed nodes
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps
    Monte_Carlo - Boolean, determines whether to make to random seed or not

    Output: optimal seed set, resulting spread, time for each iteration
    """
    # Threshold_list = [0] * 15 + [0.1] * 5 + [0.2] * 5 + [0.3] + [0.4] * 2 + [0.5] * 2 + [0.6, 0.7, 0.8] + [0.9] * 4 + [1] * 85
    # # shuffle the list
    # np.random.shuffle(Threshold_list)
    S, spread, timelapse, start_time = [], [], [], time.time()
    # setting list of nodes and making the node 107 the first node in odrer to compare with other methods
    # l = [i for i in range(g.vcount())]
    # l.remove(106)
    # l.insert(0, 106)
    
    # Find k nodes with largest marginal gain
    for _ in range(k):

        # Loop over nodes that are not yet in seed set to find biggest marginal gain
        best_spread = 0
        for j in set(range(g.vcount())) - set(S):
        #for j in set(l) - set(S):

            # Get the spread
            s = threshold_model(g, S + [j],  threshold, 28, mc, Monte_Carlo = True)

            # Update the winning node and spread so far
            if s[1] > best_spread:
                best_spread, node = s[1], j

        # Add the selected node to the seed set
        S.append(node)
        
        # Add estimated spread and elapsed time
        spread.append(best_spread)
        timelapse.append(time.time() - start_time)

    return(S,spread, timelapse)

def greedy_immunized_Th(g, k, threshold=0.1 ,mc=1000):

    """
    g - graph object
    k - number of seed nodes
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps
    Monte_Carlo - Boolean, determines whether to make to random seed or not

    Output: optimal seed set, resulting spread, time for each iteration
    """

    immunaized_list, spread, timelapse, start_time = [],[], [], time.time()
    
    # Find k nodes with largest marginal gain
    for _ in range(k):

        # Loop over nodes that are not yet in seed set to find biggest marginal gain
        best_spread = threshold_model([106], g, threshold, [], mc, Monte_Carlo = True)[1]
        for j in set(range(g.vcount())) - set(immunaized_list):

            # Get the spread
            s = threshold_model([106], g, threshold, immunaized_list + [j], mc, Monte_Carlo = True)

            # Update the winning node and spread so far
            if s[1] < best_spread and s[1] !=0 :
                best_spread, node = s[1], j

        # Add the selected node to the seed set
        if 'node' in locals():
            immunaized_list.append(node)
        
        # Add estimated spread and elapsed time
        spread.append(best_spread)
        timelapse.append(time.time() - start_time)

    return(immunaized_list, spread, timelapse)


def celf(g,k,threshold,mc):  
    """
    Input:  graph object, number of seed nodes
    Output: optimal seed set, resulting spread, time for each iteration
    """
     
    # --------------------
    # Find the first node with greedy algorithm
    # --------------------
    
    # Calculate the first iteration sorted list
    start_time = time.time() 
    marg_gain = [threshold_model(g,[node], threshold, mc)[0] for node in range(g.vcount())]

    # Create the sorted list of nodes and their marginal gain 
    Q = sorted(zip(range(g.vcount()), marg_gain), key=lambda x: x[1],reverse=True)

    # Select the first node and remove from candidate list
    S, spread, SPREAD = [Q[0][0]], Q[0][1], [Q[0][1]]
    Q, timelapse = Q[1:], [time.time()-start_time]
    #  LOOKUPS  = [g.vcount()]
    
    # --------------------
    # Find the next k-1 nodes using the list-sorting procedure
    # --------------------
    
    for _ in range(k-1):    

        check, node_lookup = False, 0
        
        while not check:
            
            # Count the number of times the spread is computed
            node_lookup += 1
            
            # Recalculate spread of top node
            current = Q[0][0]
            
            # Evaluate the spread function and store the marginal gain in the list
            Q[0] = (current, threshold_model(g, S + [current], threshold, mc)[0] - spread)
            # print(Q)
            # Re-sort the list
            Q = sorted(Q, key = lambda x: x[1], reverse = True)

            # Check if previous top node stayed on top after the sort
            check = (Q[0][0] == current)

        # Select the next node
        spread += Q[0][1]
        S.append(Q[0][0])
        SPREAD.append(spread)
        # LOOKUPS.append(node_lookup)
        timelapse.append(time.time() - start_time)

        # Remove the selected node from the list
        Q = Q[1:]

    return(S,SPREAD,timelapse) #,LOOKUPS)


