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

def IC(g,S,p=0.5,mc=1000, timestamps = 28):
    """
    Input:  
    g - graph object, 
    S - set of seed nodes(dtype list)
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps

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
            Active += new_active
            
        spread.append(len(Active))
        sum_spread.append(sum)
        
    return(np.mean(spread), np.mean(sum_spread))



def IC_immunized(g, S, immunized, p=0.15, mc=1000, timestamps = 28):
    """
    Input:  
    g - graph object, 
    S - set of seed nodes(dtype list)
    immunized - set of immunized nodes(dtype list)
    p - propagation probability
    mc - the number of Monte-Carlo simulations
    timestamps - the number of timestamps

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
                np.random.seed(i)

                # Determine neighbors that become infected
                success = np.random.uniform(0,1,len(g.neighbors(node,mode="out"))) < p
                new_ones += list(np.extract(success, g.neighbors(node,mode="out")))
                new_ones = list(set(new_ones).difference(set(immunized)))
            
            new_active = list(set(new_ones) - set(Active))
            infected_per_day.append(len(new_active))
            
            # in case the network is fully activated
            if new_active == []:
                break

            sum += len(Active) 
            # Add newly activated nodes to the set of activated nodes
            Active += new_active
            
        spread.append(len(Active))
        sum_spread.append(sum)
        infected_a_day.append(infected_per_day)
        
    return(np.mean(spread), np.mean(sum_spread), compute_mean_by_index(infected_a_day))


def greedy(g, k, p=0.15,mc=1000, timestamps = 28):
    """
    Input:  
    g - graph object
    k - number of seed nodes
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
            s = IC(g, S + [j], p, mc, timestamps)

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
    Input:  graph object, number of seed nodes
    Output: optimal seed set, resulting spread, time for each iteration
    """

    S, immunaized_list, spread, timelapse, start_time = [], [],[], [], time.time()
    
    # Find k nodes with largest marginal gain
    for _ in range(k):

        # Loop over nodes that are not yet in seed set to find biggest marginal gain
        best_spread = 122
        for j in set(range(g.vcount())) - set(immunaized_list):

            # Get the spread
            s = IC_immunized(g, [106], immunaized_list +  [j], p, mc)

            # Update the winning node and spread so far
            if s[1] < best_spread and s[1] !=0 :
                best_spread, node = s[1], j

        # Add the selected node to the seed set
        immunaized_list.append(node)
        
        # Add estimated spread and elapsed time
        spread.append(best_spread)
        timelapse.append(time.time() - start_time)

    return(immunaized_list, spread, timelapse)