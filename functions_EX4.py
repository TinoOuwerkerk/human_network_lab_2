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

def calculate_adopted_nei(node: int, node_status: list, each_neighbors: dict) -> float:
    """
    Calculate the percentage of adopted neighbors for a given node.

    Parameters:
    -----------
    node : int
        The index of the node for which to calculate the percentage of adopted neighbors.
    node_status : list
        A list representing the status (0 or 1) of each node in the network.
    each_neighbors : dict
        A dictionary where the keys are node indices and the values are lists of indices of neighbors of each node.

    Returns:
    --------
    float
        The percentage of adopted neighbors for the given node.
    """

    adopted_nei = [activated for activated in each_neighbors[node] if node_status[activated] == 1]
    return len(adopted_nei)

def th_model(seed: list, network: Graph, threshold: list, immunized, mc, Monte_Carlo = True) -> tuple:
    """
    Implementation of the Threshold Model using igraph in Python.

    Parameters:
    -----------
    node_seed : list
        The index of the seed node for starting the model.
    network : ig.Graph
        The input graph object representing the network.
    threshold : list
        The threshold value for determining whether a node gets influenced by its neighbors.

    Returns:
    --------
    tuple
        A tuple containing two elements:
        - A list of integers representing the total number of active people by the end of each day.
        - A list of lists of integers representing the indices of newly infected nodes for each day.
    """
    spread, sum_spread, infected_a_day= [], [], []
    for i in range(mc):
        if Monte_Carlo == True: 
            np.random.seed(i)
        # Prepare input for the 'calculate_adopted_nei' function
        adj_matrix = network.get_adjacency()
        adj_matrix = np.array(adj_matrix.data)
        each_neighbors = adj_matrix.nonzero()
        
        each_neighbors = {node: each_neighbors[1][each_neighbors[0] == node].tolist() for node in range(network.vcount())}

        n_node = network.vcount()
        node_status = [0] * n_node 
        neighbour_status = [0] * n_node
        new_infected = []
        day_total_infected = [0] * 28

        # Day 1
        day = 1
        sum_of_ifected = 0 
        infected_per_day = []
        for node in seed:
            node_status[node] = 1
        
        new_infected.append(seed)

        day_total_infected[day - 1] = sum(node_status)
        

        for day in range(2, 29):
            not_adopted = [node for node in range(n_node) if node_status[node] == 0]
            adopted = [node for node in range(n_node) if node_status[node] == 1]

            neighbour_status = {node : calculate_adopted_nei(node, node_status, each_neighbors) for node in range(n_node)}
            infected = []
            for node in adopted:
                for neighbour in g.neighbors(node):
                    if neighbour_status[neighbour] >= threshold[node]: 
                        infected.append(neighbour)
            
            new_infected.append(infected)
            for node in new_infected[day - 1]:
                node_status[node] = 1
                
            sum_of_ifected += len(adopted)
            day_total_infected[day - 1] = sum(node_status)
            infected_per_day.append(len(new_infected[day - 1]))

    spread.append(day_total_infected[-1])
    sum_spread.append(sum_of_ifected)
    infected_a_day.append(infected_per_day)

    return np.mean(spread), np.mean(sum_of_ifected), compute_mean_by_index(infected_a_day)

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
            s = th_model(S + [j], g, threshold, [], mc, Monte_Carlo = True)

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
        best_spread = th_model([106], g, threshold, [], mc, Monte_Carlo = True)[1]
        for j in set(range(g.vcount())) - set(immunaized_list):

            # Get the spread
            s = th_model([106], g, threshold, immunaized_list + [j], mc, Monte_Carlo = True)

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

