#!/usr/bin/env python3

'''
Date: 16sep2026

This code combines sim_frag_clust_diff_mean_dists_24apr2024.py (cluster growth with
fragmentation based on an edge-degree threshold) with the cell-death logic from
sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py.

Fragmentation happens when the highest edge degree of the cluster gets higher than a
threshold. The code grows only one cluster at a time, and after fragmentation, one of
the two resulting components ('parent', the biggest one, 'propagule', the smallest one,
or a 'random' choice between the two) is kept growing while the other is discarded.

After each cell division, both the mother cell and its new daughter cell independently
have a probability (death_prob) of dying. Dead cells remain part of the network
structure (they are not removed), but they no longer divide.

Inputs:
	-i: first_dist_params: lognormal parameters (log_mean,log_std) of the first cell division
	-j: second_dist_params: lognormal parameters (log_mean,log_std) of subsequent (2+) cell divisions
	-n: number_sims: number of the simulation/network, used to name files differently when saving them
	-o: output_dir: output directory
	-g: generations: number of generations to end the simulation
	-e: edge_degree: edge degree threshold of when fragmentation should happen
	-c: cluster_to_keep: cluster kept after fragmentation ('parent', 'propagule', 'random')
	-p: death_prob: probability of cell death after each division (0-1)

Outputs:
	-degree distribution file (saved as probabilities)
	-difference doubling time between mother and daughter cell
	-fragmentation information (edges removed, size at fracture, proportion of the biggest
	 fragment, identity of the discarded cluster, and how many of the 2 cells forming the
	 removed edge were alive at the moment of fragmentation: 'both', 'one' or 'none')
	-networks diameter, number of mother cells with 2+ undivided daughter cells, and number
	 of filamentous branches, all three computed on the unbroken network (the network as it
	 was right before fragmentation)
	-updates to fragmentation, saves how many division events happened between fragmentations
	-Network files (graphml) saved every multiple-of-10 generation, for the first network to reach it
	 (node attribute "alive" is included automatically since it is stored on the graph nodes)

Modifications from the original sim_frag_clust_diff_mean_dists_24apr2024.py:
16sep2026:
- Added -p/--death_prob input and cell death logic (ported from
  sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py): each cell (mother and
  daughter) independently survives a division event with probability (1 - death_prob).
  A cell that dies is flagged alive=False and is no longer scheduled for future
  divisions, but it is NOT removed from the network.
- Added the "alive_cells" column to fragmentation_inf.csv. Within break_cluster, right
  before the fragmenting edge is removed, the alive status of the two nodes in
  edge_to_remove is checked and summarized as 'both', 'one' or 'none'.
- Added "cases_mother_with_undivided_cells" and "num_filamentous_branches" columns to
  networks_diameter.csv (ported from sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py),
  computed at the same point and on the same network as the diameter (the unbroken
  network, right before fragmentation).
- Added a retry-on-extinction mechanism: since cells can now die, it's possible (though
  unlikely at typical death_prob values) for a whole cluster to die out before reaching
  -g generations. Rather than buffering everything in memory like the no-fragmentation
  death-probability script does, this version keeps writing directly to disk as it goes,
  but bookmarks each output file's write position (and which network-saving flags/graphml
  files were touched) at the start of every attempt. If an attempt goes extinct early, all
  files are truncated back to their bookmarked position, any graphml files written during
  that attempt are deleted, and the same simulation number is retried from scratch.
'''

import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
from random import choice
import random
import numpy as np
import time

#only used for testing purposes
from itertools import permutations

parser = argparse.ArgumentParser()
parser.add_argument('-i','--first_dist_params',dest="first_dist_params",required=True) #mean and variance of lognormal distribution for first division
parser.add_argument('-j','--second_dist_params',dest="second_dist_params",required=True) #mean and variance of lognormal distribution for second or more divisions
parser.add_argument('-n','--number_sims',dest="number_sims",required=True) #number of simulations it is going to run
parser.add_argument('-o','--output_dir',dest="output_dir",required=True) #output directory
parser.add_argument('-g','--generations',dest="generations",required=True) #Number of generations
parser.add_argument('-e','--edge_degree',dest="edge_degree",required=True) #edge degree threshold
parser.add_argument('-c','--cluster_to_keep',dest="cluster_to_keep",required=True) # 'parent', 'propagule', 'random'
parser.add_argument('-p','--death_prob',dest="death_prob",required=True) #probability of cell death after each division (0-1)
args = parser.parse_args()

### Global variables ####
list_flags_saved_network=[]

#### Functions for simulating networks growth ########

#sampling the doubling time of the next division using a log normal distribution
def sample_doub_t(first_div_params, second_div_params, division_cells):
	if(division_cells==0):
		dist_mu=first_div_params[0]
		dist_sigma=first_div_params[1]
	else:
		dist_mu=second_div_params[0]
		dist_sigma=second_div_params[1]
	
	return(round(float(np.random.lognormal(dist_mu, dist_sigma, 1)[0]), 4))

# Modified function to add cells to ordered list by time (similar to no-fragmentation version)
def add_cells_list(input_time, cell_id, ordered_list):
	if not ordered_list:
		ordered_list.append([input_time, [cell_id]])
		return ordered_list
	
	for i, el in reversed(list(enumerate(ordered_list))):
		if input_time == el[0]:
			ordered_list[i][-1].append(cell_id)
			return ordered_list
		elif input_time > el[0]:
			ordered_list.insert(i+1, [input_time, [cell_id]])
			return ordered_list
	
	ordered_list.insert(0, [input_time, [cell_id]])
	return ordered_list

# Modified function to divide all cells at the current timepoint. Now also applies the
# cell death probability check to the mother and the daughter cell right after the
# division event (ported from sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py).
# A cell that dies (alive=False) is NOT scheduled for any future division, but it stays
# in the network as a node.
def divide_cells(cluster_population, cells_to_divide, cont_ids, curr_time):
	first_div_params = input_variables['first_div_params']
	second_div_params = input_variables['second_div_params']
	sim_number = input_variables['sim_number']
	death_prob = input_variables['death_prob']
	
	list_ids = cells_to_divide[0][1].copy()  # Make a copy to avoid modifying the original list
	
	for mother_id in list_ids:
		# Get mother cell identity before division
		temp_mother_identity = cluster_population.nodes[mother_id]['identity']
		
		# Create daughter cell
		daughter_id = cont_ids
		cluster_population.add_node(daughter_id, number_divisions=0, identity=temp_mother_identity, alive=True)
		cluster_population.add_edge(mother_id, daughter_id)
		temp_daughter_time = sample_doub_t(first_div_params, second_div_params, 0)

		# Check if daughter cell survives (cell death probability check)
		daughter_survives = random.random() >= death_prob
		cluster_population.nodes[daughter_id]["alive"] = daughter_survives

		# Only schedule the daughter's next division if it survives
		if daughter_survives:
			cells_to_divide = add_cells_list(curr_time + temp_daughter_time, daughter_id, cells_to_divide)
		cont_ids += 1

		# Update mother cell
		cluster_population.nodes[mother_id]["number_divisions"] += 1
		temp_mother_div = cluster_population.nodes[mother_id]["number_divisions"]
		temp_mother_time = sample_doub_t(first_div_params, second_div_params, temp_mother_div)

		# Check if mother cell survives (cell death probability check)
		mother_survives = random.random() >= death_prob
		cluster_population.nodes[mother_id]["alive"] = mother_survives

		# Only schedule the mother's next division if it survives
		if mother_survives:
			cells_to_divide = add_cells_list(curr_time + temp_mother_time, mother_id, cells_to_divide)

		# Save doubling information (kept for both alive and dead cells, same as before)
		dict_files['sampled_times'].write(sim_number+','+str(temp_mother_div)+','+str(temp_mother_time)+'\n')
		dict_files['sampled_times'].write(sim_number+','+str(0)+','+str(temp_daughter_time)+'\n')
		
		# Save difference in doubling time
		cont_gen = input_variables.get('current_generation', 0)
		dict_files['diff_d_t'].write(sim_number+","+str(cont_gen)+","+str(abs(temp_mother_time-temp_daughter_time))+"\n")

	# Remove the processed timepoint
	cells_to_divide.pop(0)
	
	return [cont_ids, cells_to_divide]

#Function to save the degree distribution of the network, the degree distribution is saved as a
#probability distribution to be normalized for network size
def save_degree_distribution_to_csv(network, num_generations):
	# Get the degree of all nodes
	degrees = dict(network.degree())
	
	# Create a DataFrame from the degrees dictionary
	df_degrees = pd.DataFrame(degrees.items(), columns=['Node', 'Degree'])
	
	# Create a DataFrame for the degree distribution
	degree_counts=df_degrees['Degree'].value_counts().reset_index()
	degree_counts.columns=['Degree', 'Frequency']
	
	# Calculate the probability distribution
	total_nodes = len(df_degrees)
	degree_counts['Probability']=degree_counts['Frequency'] / total_nodes
	
	# Add static columns to the dataframe
	degree_counts['num_generations']=num_generations
	degree_counts['sim_number']=input_variables['sim_number']

	for i in range(len(degree_counts)):
		dict_files['degree_dist'].write(",".join([str(j) for j in degree_counts.iloc[i]])+"\n")

#Function to remove the next doubling times of the cells of the parent propagule that is being
#discarded
def remove_cells_by_component(list_cells, component):
	positions_to_remove = []
	
	# Iterate through the list of lists and save the positions to remove
	for i, cell_data in enumerate(list_cells):
		# Check if any cell in the timepoint belongs to the component
		cells_in_timepoint = cell_data[1]
		remaining_cells = [cell for cell in cells_in_timepoint if cell not in component]
		
		if len(remaining_cells) == 0:
			# All cells in this timepoint should be removed
			positions_to_remove.append(i)
		elif len(remaining_cells) < len(cells_in_timepoint):
			# Some cells should be removed, update the list
			list_cells[i][1] = remaining_cells
	
	# Remove empty timepoints in reverse order to avoid index issues
	for pos in reversed(positions_to_remove):
		list_cells.pop(pos)
	
	return list_cells

#Counts mother cells that have 2 or more undivided (degree-1) daughter cells
#(ported from sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py)
def count_mothers_with_multiple_undivided(network):
	mothers_with_multiple_undivided = 0
	for node in network.nodes():
		neighbors = list(network.neighbors(node))
		undivided_daughters = sum(1 for neighbor in neighbors if network.degree(neighbor) == 1)
		if undivided_daughters >= 2:
			mothers_with_multiple_undivided += 1
	return mothers_with_multiple_undivided

#Identifies nodes that are part of a filamentous branch (chains of degree<=2 nodes)
#(ported from sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py)
def identify_filament_nodes(graph, min_length=3):
	filament_nodes = set()
	visited = set()
	
	def traverse_filament(node):
		path = [node]
		current = node
		visited.add(current)
		while True:
			neighbors = list(graph.neighbors(current))
			unvisited_neighbors = [n for n in neighbors if n not in visited]
			if len(unvisited_neighbors) == 1 and graph.degree(unvisited_neighbors[0]) <= 2:
				next_node = unvisited_neighbors[0]
				path.append(next_node)
				visited.add(next_node)
				current = next_node
			else:
				break
		return path if len(path) >= min_length else []
	
	for node, degree in dict(graph.degree()).items():
		if degree == 1 and node not in visited:
			filament = traverse_filament(node)
			filament_nodes.update(filament)
	
	return filament_nodes

#Counts the number of filamentous branches in the network
#(ported from sim_clust_no_fragmentation_lognormal_death_prob_5aug2025.py)
def count_filamentous_branches(graph, min_length=3):
	filament_nodes = identify_filament_nodes(graph, min_length)
	visited = set()
	branch_count = 0
	
	for node in filament_nodes:
		if node not in visited and graph.degree(node) == 1:
			# Start from a terminal node and traverse the filament
			current = node
			path = [current]
			visited.add(current)
			
			while True:
				neighbors = list(graph.neighbors(current))
				unvisited_neighbors = [n for n in neighbors if n not in visited and n in filament_nodes]
				if len(unvisited_neighbors) == 1:
					next_node = unvisited_neighbors[0]
					path.append(next_node)
					visited.add(next_node)
					current = next_node
				else:
					break
			
			if len(path) >= min_length:
				branch_count += 1
	
	return branch_count

#function to save intermediate networks
def save_intermediate_networks(network, temp_generation_cluster):
	output_dir=input_variables["output_dir"]
	sim_number=input_variables["sim_number"]
	network_dir=input_variables["network_dir"]

	global list_flags_saved_network

	#save cluster before if it is the first cluster of this generation
	if(temp_generation_cluster%10==0):
		if(list_flags_saved_network[int(temp_generation_cluster/10)]==False):
			nx.write_graphml(network, os.path.join(network_dir,str(sim_number)+"_"+str(temp_generation_cluster)+".graphml"))
			
			list_flags_saved_network[int(temp_generation_cluster/10)]=True

	#calculate network diameter for each network and append it into the file
	temp_diameter_network=str(nx.diameter(network))

	#calculate number of mother cells with 2+ undivided daughter cells, and number of
	#filamentous branches, both on this same (unbroken) network used for the diameter
	mothers_with_undivided = count_mothers_with_multiple_undivided(network)
	num_filamentous_branches = count_filamentous_branches(network)

	dict_files['diameter'].write(sim_number+","+str(temp_generation_cluster)+","+temp_diameter_network+","+str(mothers_with_undivided)+","+str(num_filamentous_branches)+"\n")

	save_degree_distribution_to_csv(network, temp_generation_cluster)

# Checks the alive status of the two nodes that form the edge that is about to be
# removed during fragmentation, and summarizes it as 'both', 'one' or 'none' alive.
# This is checked on the intact cluster, before the edge is removed or node identities
# are reassigned to 'parent'/'propagule'.
def get_alive_status_of_edge(network, edge_to_remove):
	node1, node2 = edge_to_remove
	alive_node1 = network.nodes[node1]['alive']
	alive_node2 = network.nodes[node2]['alive']

	num_alive = sum([alive_node1, alive_node2])

	if num_alive == 2:
		return 'both'
	elif num_alive == 1:
		return 'one'
	else:
		return 'none'

def save_fragmentation_summary(propagule, edge_to_remove, temp_generation_cluster, size_network, curr_identity, alive_cells):
	output_dir=input_variables["output_dir"]
	sim_number=input_variables["sim_number"]

	len_propagule=len(propagule)
	prop_propagule=str(len_propagule/size_network)

	dict_files['fragmentation'].write(sim_number+","+str(temp_generation_cluster)+","+str(edge_to_remove[0])+"-"+str(edge_to_remove[1])+","+str(size_network)+","+str(len_propagule)+","+prop_propagule+","+curr_identity+","+alive_cells+"\n")

#Function to calculate the highest edge degree of a network, it returns the edge and highest edge degree value
def calculate_max_edge_degree(network):
	#initialize result variables
	max_edge_degree = -1
	max_edge = None

	for edge in network.edges():
		node1, node2 = edge
		degree_node1 = network.degree(node1)
		degree_node2 = network.degree(node2)
		edge_degree = degree_node1 + degree_node2 - 2

		if (edge_degree > max_edge_degree):
			max_edge_degree = edge_degree
			max_edge = edge

	return([max_edge, max_edge_degree])

#Function that does the fragmentation of the cluster and saves information about the cluster
def break_cluster(snowflake, cells_to_divide, temp_generation_cluster, cont_updates, edge_to_remove, max_num_generations):
	num_generations=input_variables["num_generations"]
	output_dir=input_variables['output_dir']
	sim_number=input_variables["sim_number"]
	cluster_to_keep=input_variables['cluster_to_keep']

	#calculate amount of cells in the cluster
	size_network=snowflake.number_of_nodes()

	#Find edge to remove according to fracture proportion
	#creating copy of the cluster
	unbroken_cluster=snowflake.copy()

	#Check the alive status of the 2 nodes forming the edge that is about to be removed,
	#before the edge is removed or identities are reassigned
	alive_cells=get_alive_status_of_edge(snowflake, edge_to_remove)

	#removing the edge in snowflake variable
	snowflake.remove_edge(*edge_to_remove) #this is updated to the variable outside the function

	curr_identity=snowflake.nodes[list(snowflake.nodes())[0]]['identity']

	#Selecting the component which is going to be removed
	if(cluster_to_keep=='parent'):
		component_to_remove=sorted(nx.connected_components(snowflake), key=len)[0]
		for node in snowflake.nodes():
			snowflake.nodes[node]['identity']='parent'
	elif(cluster_to_keep=='propagule'):
		component_to_remove=sorted(nx.connected_components(snowflake), key=len)[-1]
		for node in snowflake.nodes():
			snowflake.nodes[node]['identity']='propagule'
	else:
		if(random.choice(['parent', 'propagule'])=='propagule'):
			component_to_remove=sorted(nx.connected_components(snowflake), key=len)[-1]
			for node in snowflake.nodes():
				snowflake.nodes[node]['identity']='propagule'
		else:
			component_to_remove=sorted(nx.connected_components(snowflake), key=len)[0]
			for node in snowflake.nodes():
				snowflake.nodes[node]['identity']='parent'

	#obtain the propagule
	propagule = sorted(nx.connected_components(snowflake), key=len)[0]

	#save fragmentation summary information
	save_fragmentation_summary(propagule, edge_to_remove, temp_generation_cluster, size_network, curr_identity, alive_cells)

	#Removing cells from the selected component
	cells_to_divide=remove_cells_by_component(cells_to_divide, component_to_remove)
	snowflake.remove_nodes_from([node for node in snowflake.nodes() if node in component_to_remove])

	#save network properties before fracture
	save_intermediate_networks(unbroken_cluster, temp_generation_cluster)

	#Save how many updates have passed
	dict_files['snowflakes_update'].write(sim_number+","+str(temp_generation_cluster)+","+str(cont_updates)+"\n")

	# to finish the simulation the cells_to_divide variable is cleared
	if(temp_generation_cluster>max_num_generations):
		cells_to_divide=[]

	return(cells_to_divide)

#Main loop where the growth of one cluster is simulated until it reaches the desired amount of
#generations
def simulate_one_cluster_growth(input_variables):
	#defining general use variables
	output_dir=input_variables['output_dir']
	num_generations=input_variables['num_generations']
	sim_number=input_variables['sim_number']
	edge_degree_threshold=input_variables['edge_degree_threshold']
	first_div_params=input_variables['first_div_params']
	second_div_params=input_variables['second_div_params']

	cont_ids=1

	#creating graph and adding first cell
	snowflake=nx.Graph() #creating graph
	snowflake.add_node(cont_ids,number_divisions=0, identity='parent', alive=True) #adding the first cell
	
	# Initialize cells to divide with new structure: [[time, [cell_ids]], ...]
	temp_next_doub = sample_doub_t(first_div_params, second_div_params, 0)
	cells_to_divide = []
	cells_to_divide = add_cells_list(temp_next_doub, cont_ids, cells_to_divide)
	cont_ids += 1

	t_sim = 0
	cont_updates = 1
	cont_gen = 0

	# Flag to track whether the simulation actually reached the target number of
	# generations, or whether it died out early (extinction). It only ever becomes
	# True inside break_cluster's own "reached the end" check below, mirroring the
	# condition that clears cells_to_divide.
	reached_target = False

	# Store current generation in input_variables for access in divide_cells
	input_variables['current_generation'] = cont_gen

	while len(cells_to_divide) > 0:
		# Update simulation time
		curr_time = cells_to_divide[0][0]
		t_sim = curr_time
		
		# Update current generation in input_variables
		input_variables['current_generation'] = cont_gen

		# Divide all cells at the current timepoint
		cont_ids, cells_to_divide = divide_cells(snowflake, cells_to_divide, cont_ids, curr_time)

		# Check for fragmentation after all cells at this timepoint have divided
		edge_to_remove, edge_degree = calculate_max_edge_degree(snowflake)

		if(edge_degree >= edge_degree_threshold):
			# Test fragmentation if enabled
			if(input_variables['test_frag']):
				is_fractured_edge_between_highest_degree_nodes(snowflake, edge_to_remove, cont_gen)

			# Break cluster and update cells_to_divide
			if(cont_gen > num_generations):
				reached_target = True
			cells_to_divide = break_cluster(snowflake, cells_to_divide, cont_gen, cont_updates, edge_to_remove, num_generations)
			cont_gen += 1

		cont_updates += 1

	# True if the target number of generations was reached; False if the cluster died
	# out (all cells stopped dividing) before getting there
	return reached_target

# Rest of the functions remain the same as in original code...

'''
Steps to find if the nodes of highest degree are forming the edge that is getting fractured:
1.- Find nodes of highest edge degree
2.- Get all edges combinations of the highest degree nodes
3.- Compare all possible pairs to the actual edge deleted
'''
def is_fractured_edge_between_highest_degree_nodes(network, edge_to_remove, cont_gen):
	sim_number=input_variables["sim_number"]

	#1
	highest_degree_nodes=get_nodes_with_highest_degrees(network)

	#2
	list_possible_edges=pairwise_permutations(highest_degree_nodes)

	#3
	temp_flag=str(edge_in_list(edge_to_remove, list_possible_edges))

	#write answer in input file
	dict_files['test_frag'].write(sim_number+","+str(cont_gen)+','+temp_flag+'\n')

#1.- Get list of nodes of the highest degree
def get_nodes_with_highest_degrees(graph):
	# Get a dictionary with nodes as keys and degrees as values
	degree_dict = dict(graph.degree())

	# Sort the dictionary by degree in descending order
	sorted_degree = sorted(degree_dict.items(), key=lambda x: x[1], reverse=True)

	# Find the highest degree value
	highest_degree = sorted_degree[0][1]

	# Initialize lists for highest and second-highest degree nodes
	highest_degree_nodes = []
	second_highest_degree_nodes = []
	second_highest_degree=-1

	# Iterate through the sorted_degree list
	for node, degree in sorted_degree:
		if degree == highest_degree:
			highest_degree_nodes.append(node)
			if len(highest_degree_nodes) >= 2:
				break
		elif not second_highest_degree_nodes and degree < highest_degree:
			second_highest_degree_nodes.append(node)
			second_highest_degree=degree
		elif degree == second_highest_degree:
			second_highest_degree_nodes.append(node)

	# If there are at least 2 nodes in highest_degree_nodes, return only them
	if len(highest_degree_nodes) >= 2:
		return (highest_degree_nodes)
	# If there are nodes with the second-highest degree, return them along with the highest_degree_nodes
	elif second_highest_degree_nodes:
		return(highest_degree_nodes + second_highest_degree_nodes)
	# If there are no nodes with the second-highest degree, return highest_degree_nodes
	else:
		return(highest_degree_nodes)

#2.- Get all edges combinations of the highest degree nodes
def pairwise_permutations(numbers):
	pairs = list(permutations(numbers, 2))
	return pairs + [(y, x) for (x, y) in pairs]

#3.- Compare all possible pairs to the actual edge deleted
def edge_in_list(edge, list_edges):
	for temp_edge in list_edges:
		if(temp_edge==edge):
			return(True) #if it was found return True
	
	return(False) #If it wasn't found return False

#Function to initialize all output file, a dictionary with all the opened files is returned for it
#to be accesed by all functions
def initialize_files(input_variables):
	output_dir=input_variables['output_dir']

	dict_files={}

	#Initializing results file if they don't exist
	diff_d_t_file=os.path.join(output_dir,"diff_doub_t.csv")
	if(os.path.exists(diff_d_t_file)):
		dict_files['diff_d_t']=open(diff_d_t_file, "a")
	else:
		dict_files['diff_d_t']=open(diff_d_t_file, "w")
		dict_files['diff_d_t'].write('sim_number,generation_number,diff_minutes\n')

	diam_file=os.path.join(output_dir, "networks_diameter.csv")
	if(os.path.exists(diam_file)):
		dict_files['diameter']=open(diam_file, "a")
	else:
		dict_files['diameter']=open(diam_file, "w")
		dict_files['diameter'].write("sim_number,generation,diameter,cases_mother_with_undivided_cells,num_filamentous_branches\n")

	fragmentation_file=os.path.join(output_dir, "fragmentation_inf.csv")
	if(os.path.exists(fragmentation_file)):
		dict_files['fragmentation']=open(fragmentation_file, "a")
	else:
		dict_files['fragmentation']=open(fragmentation_file, "w")
		dict_files['fragmentation'].write("sim_number,generation,removed_edge,cluster_size,size_propagule,proportion_propagule,identity,alive_cells\n")

	degree_dist_file=os.path.join(output_dir,"degree_distribution.csv")
	if(os.path.exists(degree_dist_file)):
		dict_files['degree_dist']=open(degree_dist_file, "a")
	else:
		dict_files['degree_dist']=open(degree_dist_file, "w")
		dict_files['degree_dist'].write("Degree,Frequency,Probability,num_generations,sim_number\n")

	sampled_times_file=os.path.join(output_dir,"sampled_times.csv")
	if(os.path.exists(sampled_times_file)):
		dict_files['sampled_times']=open(sampled_times_file, "a")
	else:
		dict_files['sampled_times']=open(sampled_times_file, "w")
		dict_files['sampled_times'].write("sim_number,number_divisions,minutes"+"\n")

	updates_file=os.path.join(output_dir, "updates_to_fragmentation.csv")
	if(os.path.exists(updates_file)):
		dict_files['snowflakes_update']=open(updates_file, "a")
	else:
		dict_files['snowflakes_update']=open(updates_file, "w")
		dict_files['snowflakes_update'].write("sim_number,generation,updates\n")

	#temporary file
	if(input_variables['test_frag']):
		test_file=os.path.join(output_dir, "test_nodes_frag.csv")
		if(os.path.exists(test_file)):
			dict_files['test_frag']=open(test_file, "a")
		else:
			dict_files['test_frag']=open(test_file, "w")
			dict_files['test_frag'].write("sim_number,generation,flag\n")

	return(dict_files)

#function to close all files opened in the dictionary
def close_all_files(dict_files):
	for file_obj in dict_files.values():
		file_obj.close()

#### MAIN ####

#initializing variables
#doubling time distribution parameters
first_div_params=[float(elem) for elem in args.first_dist_params.split(",")]
second_div_params=[float(elem) for elem in args.second_dist_params.split(",")]

output_dir=args.output_dir
num_generations=int(args.generations)
num_iterations=int(args.number_sims)
edge_degree_threshold=int(args.edge_degree)

cluster_to_keep=args.cluster_to_keep.lower()
if(not cluster_to_keep in ['parent', 'propagule', 'random']):
	print("Error no valid option given in argument -c\nAll valid options are 'parent', 'offspring', 'random'")
	exit()

#death probability, must be a value between 0 and 1
death_prob=float(args.death_prob)
if(death_prob < 0 or death_prob > 1):
	print("Error no valid value given in argument -p\ndeath_prob must be a number between 0 and 1")
	exit()

#creating output directory if it doesn't exist, and stopping the program execution if the 
#directory exists
if(not os.path.exists(output_dir)):
	os.mkdir(output_dir)
else:
	print("Output directory already exists, stopping program execution to avoid overwriting data")
	exit()
network_dir=os.path.join(output_dir, 'network_dir')
if(not os.path.exists(network_dir)):
	os.mkdir(network_dir)

#For testing if the edge of the fragmentation is formed by the nodes of highest edge degree
test_edge_of_fragmentation=False

#Dictionary containing input variables and parameters of the simulation
#this variable is used as a global variable
input_variables={"output_dir":output_dir, "network_dir":network_dir, "num_generations":num_generations, 
"edge_degree_threshold":edge_degree_threshold, "test_frag":test_edge_of_fragmentation, 
"first_div_params":first_div_params, "second_div_params":second_div_params, "cluster_to_keep":cluster_to_keep,
"death_prob":death_prob}

#Initializing all output files. Note: this variable is accesed as a global variable
dict_files=initialize_files(input_variables)

#this list is to save if the network has already been saved in the network files,
#only the first network and the final networks are going to be saved
list_flags_saved_network=[False for i in range(num_generations+1) if i%10==0]

#creating and writing information of the program in the log file
log_file=open(os.path.join(output_dir,"log.txt"), "a")
log_file.write("sim_frag_clust_death_prob_16sep2026.py\nInputs received:\n")
log_file.write("-i: "+args.first_dist_params+"\n")
log_file.write("-j: "+args.second_dist_params+"\n")
log_file.write("-n: "+args.number_sims+"\n")
log_file.write("-o: "+args.output_dir+"\n")
log_file.write("-g: "+args.generations+"\n")
log_file.write("-e: "+args.edge_degree+"\n")
log_file.write("-c: "+args.cluster_to_keep+"\n")
log_file.write("-p: "+args.death_prob+"\n")
log_file.write("\nExecution times:\n")

for i in range(1, num_iterations+1):
	input_variables["sim_number"]=str(i) #changing variable type as it is only used as a character for saving results of the simulation

	#counting execution time
	start_time = time.time()

	#Keep retrying this simulation number until it reaches the target number of
	#generations. Each output file's write position, and which network-saving flags
	#are already set, are bookmarked before every attempt so that a failed
	#(extinct) attempt can be fully erased before retrying.
	attempt = 1
	while True:
		file_bookmarks = {}
		for file_name, file_obj in dict_files.items():
			file_obj.flush()
			file_bookmarks[file_name] = file_obj.tell()
		flags_snapshot = list_flags_saved_network.copy()

		#grow the cluster of this simulation
		reached_target = simulate_one_cluster_growth(input_variables)

		if reached_target:
			break

		#Extinction before reaching -g generations: erase everything written to disk
		#during this failed attempt
		for file_name, file_obj in dict_files.items():
			file_obj.flush()
			file_obj.seek(file_bookmarks[file_name])
			file_obj.truncate()

		#Remove any graphml network files written during the failed attempt, and
		#reset the corresponding "already saved" flags so the retry can save them
		for idx in range(len(list_flags_saved_network)):
			if (not flags_snapshot[idx]) and list_flags_saved_network[idx]:
				failed_network_file = os.path.join(network_dir, str(i)+"_"+str(idx*10)+".graphml")
				if os.path.exists(failed_network_file):
					os.remove(failed_network_file)
				list_flags_saved_network[idx] = False

		attempt += 1

	end_time=time.time()
	elapsed_time=end_time-start_time

	if attempt > 1:
		log_file.write("Simulation "+str(i)+" completed in "+str(elapsed_time)+" seconds after "+str(attempt)+" attempts.\n")
	else:
		log_file.write("Simulation "+str(i)+" completed in "+str(elapsed_time)+" seconds.\n")

log_file.write("\nProgram finished execution without any errors.")
log_file.close()

#closing all files from the dictionary
close_all_files(dict_files)