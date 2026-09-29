#!/usr/bin/env python3

'''
Date: 28sep2026

This code combines sim_frag_clust_death_prob_16sep2026.py (cluster growth with an
edge-degree fragmentation threshold, lognormal doubling times, and cell death) with the
"keep all clusters" fragmentation logic of
sim_frag_clust_edge_degree_all_clusters_exponential_9mar2024.py.

Fragmentation happens when the highest edge degree of a cluster gets higher than a
threshold. Unlike sim_frag_clust_death_prob_16sep2026.py (which grows only one cluster
at a time and discards one side after every fragmentation), this version keeps BOTH
resulting components after every fragmentation event, so the population branches out
into more and more independent clusters over time, all of which keep growing/dividing/
fragmenting on their own until each one individually reaches -g generations (number of
fragmentation events undergone by that cluster's lineage), at which point that specific
cluster stops dividing.

After each cell division, both the mother cell and its new daughter cell independently
have a probability (death_prob) of dying. Dead cells remain part of the network
structure (they are not removed), but they no longer divide.

Inputs:
	-i: first_dist_params: lognormal parameters (log_mean,log_std) of the first cell division
	-j: second_dist_params: lognormal parameters (log_mean,log_std) of subsequent (2+) cell divisions
	-n: number_sims: number of the simulation/network, used to name files differently when saving them
	-o: output_dir: output directory
	-g: generations: number of (per-cluster) fragmentation generations to end each cluster's growth
	-e: edge_degree: edge degree threshold of when fragmentation should happen
	-p: death_prob: probability of cell death after each division (0-1)

Outputs:
	-difference in doubling time between mother and daughter cell
	-fragmentation information (edges removed, size at fracture, proportion of the propagule,
	 identity of the cluster before fracturing, how many of the 2 cells forming the removed
	 edge were alive at the moment of fragmentation: 'both', 'one' or 'none', and cluster_id/
	 parent_id so lineages can be tracked)
	-networks diameter, number of mother cells with 2+ undivided daughter cells, and number
	 of filamentous branches, all three computed on the unbroken network (the network as it
	 was right before that cluster's fragmentation), plus the cluster_id
	-sampled doubling times
	-Network files (graphml) saved every multiple-of-10 generation, for the first cluster
	 (across the whole run) to reach it (node attribute "alive" is included automatically
	 since it is stored on the graph nodes)

Not produced (not needed for this version):
	-updates_to_fragmentation.csv ("snowflakes_update")
	-degree_distribution.csv ("degree_dist")

Differences from sim_frag_clust_death_prob_16sep2026.py:
- Removed the -c/--cluster_to_keep input entirely: there is no "discard one side" step.
  After every fragmentation event both resulting components keep dividing.
- All clusters live together in a single networkx graph (cluster_population), the same
  way sim_frag_clust_edge_degree_all_clusters_exponential_9mar2024.py does it. Each node
  carries cluster_id/parent_id/cluster_generation attributes (ported from that script) in
  addition to alive/number_divisions/identity (ported from the death_prob script).
- Fragmentation is now checked across every currently-existing cluster each timepoint
  (fragmentation_edge_degree_population, ported from the all_clusters script), instead of
  a single break_cluster call for one lineage. Cell-death bookkeeping (alive flag, only
  scheduling surviving cells for future division) and the alive_cells ('both'/'one'/'none')
  fragmentation-edge check from the death_prob script were folded into this per-cluster loop.
- "generation" keeps the same meaning it already had in both source scripts: the number of
  fragmentation events a given cluster's lineage has undergone. Because a cluster is only
  removed from the simulation once its OWN lineage reaches -g generations, different
  clusters can be at different generations at the same simulated time. The diff_doub_t.csv
  "generation_number" column reflects the generation of the specific cluster the dividing
  mother cell belonged to at the time of division (there is no single shared generation
  counter any more, since many clusters coexist).
- Retry-on-extinction, but only checked before the first fragmentation: dying out entirely
  is only realistically possible in the earliest generations, while the population is still
  a single small cluster that hasn't yet reached the edge-degree threshold once. Once at
  least one fragmentation event has happened, the run is kept as-is even if some later
  lineage happens to die out from cell death before reaching -g generations (that lineage's
  branch just stops a bit short; there are still plenty of other fragmentations from the run
  to study, so it isn't worth discarding and re-running the whole thing over it). So
  reached_target is simply whether at least one fragmentation event occurred anywhere in the
  population before cells_to_divide emptied out: reached_target = any_fragmentation_happened.
  Exactly like the death_prob script, every output file's write position (and which
  network-saving flags/graphml files were touched) is bookmarked at the start of each
  attempt, so a failed (pre-first-fragmentation extinction) attempt can be fully erased on
  disk before the same simulation number is retried from scratch.
- Dropped the degree-distribution and snowflakes-update outputs and their files entirely, per
  request. Dropped the now-unused pandas import (only used for degree distribution) and the
  cluster_to_keep-related "from random import choice" import.
'''

import os
import argparse
import networkx as nx
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
parser.add_argument('-g','--generations',dest="generations",required=True) #Number of (per-cluster) generations
parser.add_argument('-e','--edge_degree',dest="edge_degree",required=True) #edge degree threshold
parser.add_argument('-p','--death_prob',dest="death_prob",required=True) #probability of cell death after each division (0-1)
args = parser.parse_args()

### Global variables ####
list_flags_saved_network=[]
list_saved_network_filepaths=[]

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

#function to add cells to ordered list by time
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

#Function to divide all cells at the current timepoint. Carries cluster_generation/
#cluster_id/parent_id/identity from mother to daughter (ported from the all_clusters
#script) and applies the cell death probability check to the mother and the daughter
#cell right after the division event (ported from sim_frag_clust_death_prob_16sep2026.py).
#A cell that dies (alive=False) is NOT scheduled for any future division, but it stays
#in the network as a node.
def divide_cells(cluster_population, cells_to_divide, cont_ids):
	first_div_params = input_variables['first_div_params']
	second_div_params = input_variables['second_div_params']
	sim_number = input_variables['sim_number']
	death_prob = input_variables['death_prob']

	curr_time = cells_to_divide[0][0]
	list_ids = cells_to_divide[0][1].copy()  # Make a copy to avoid modifying the original list

	for mother_id in list_ids:
		# Get mother cell identity/lineage info before division
		temp_mother_identity = cluster_population.nodes[mother_id]['identity']
		mother_clust_gen = cluster_population.nodes[mother_id]['cluster_generation']
		mother_clust_id = cluster_population.nodes[mother_id]['cluster_id']
		temp_parent_id = cluster_population.nodes[mother_id]['parent_id']

		# Create daughter cell
		daughter_id = cont_ids
		cluster_population.add_node(daughter_id, number_divisions=0, identity=temp_mother_identity,
			cluster_generation=mother_clust_gen, cluster_id=mother_clust_id, parent_id=temp_parent_id, alive=True)
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

		# Save difference in doubling time, tagged with the generation of the cluster this
		# mother cell belonged to at the time of division
		dict_files['diff_d_t'].write(sim_number+","+str(mother_clust_gen)+","+str(abs(temp_mother_time-temp_daughter_time))+"\n")

	# Remove the processed timepoint
	cells_to_divide.pop(0)

	return [cont_ids, cells_to_divide]

#Function to remove the next doubling times of the cells belonging to a cluster that just
#got pruned (reached -g generations)
def remove_cells_by_component(cells_to_divide, component_to_remove):
	updated_list = []

	for time_entry, cells in cells_to_divide:
		filtered_cells = [cell for cell in cells if cell not in component_to_remove]
		if filtered_cells:
			updated_list.append([time_entry, filtered_cells])

	return updated_list

#Counts mother cells that have 2 or more undivided (degree-1) daughter cells
#(ported from sim_frag_clust_death_prob_16sep2026.py)
def count_mothers_with_multiple_undivided(network):
	mothers_with_multiple_undivided = 0
	for node in network.nodes():
		neighbors = list(network.neighbors(node))
		undivided_daughters = sum(1 for neighbor in neighbors if network.degree(neighbor) == 1)
		if undivided_daughters >= 2:
			mothers_with_multiple_undivided += 1
	return mothers_with_multiple_undivided

#Identifies nodes that are part of a filamentous branch (chains of degree<=2 nodes)
#(ported from sim_frag_clust_death_prob_16sep2026.py)
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
#(ported from sim_frag_clust_death_prob_16sep2026.py)
def count_filamentous_branches(graph, min_length=3):
	filament_nodes = identify_filament_nodes(graph, min_length)
	visited = set()
	branch_count = 0

	for node in filament_nodes:
		if node not in visited and graph.degree(node) == 1:
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

#function to save network properties (diameter, mothers with undivided cells, filamentous
#branches) and, every multiple-of-10 generation, the graphml file for the first cluster
#(across the whole run) to reach it
def save_intermediate_networks(network, temp_generation_cluster, temp_cluster_id):
	sim_number=input_variables["sim_number"]
	network_dir=input_variables["network_dir"]

	global list_flags_saved_network
	global list_saved_network_filepaths

	#save network if it is the first cluster to reach this generation multiple of 10
	if(temp_generation_cluster%10==0):
		idx=int(temp_generation_cluster/10)
		if(list_flags_saved_network[idx]==False):
			filepath=os.path.join(network_dir, str(sim_number)+"_"+str(temp_generation_cluster)+"_"+str(temp_cluster_id)+".graphml")
			nx.write_graphml(network, filepath)

			list_flags_saved_network[idx]=True
			list_saved_network_filepaths[idx]=filepath

	#calculate network diameter for each network and append it into the file
	temp_diameter_network=str(nx.diameter(network))

	#calculate number of mother cells with 2+ undivided daughter cells, and number of
	#filamentous branches, both on this same (unbroken) network used for the diameter
	mothers_with_undivided = count_mothers_with_multiple_undivided(network)
	num_filamentous_branches = count_filamentous_branches(network)

	dict_files['diameter'].write(sim_number+","+str(temp_generation_cluster)+","+temp_diameter_network+","+str(temp_cluster_id)+","+str(mothers_with_undivided)+","+str(num_filamentous_branches)+"\n")

# Checks the alive status of the two nodes that form the edge that is about to be
# removed during fragmentation, and summarizes it as 'both', 'one' or 'none' alive.
# This is checked on the intact (pre-fracture) cluster.
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

#function to save information about the fragmentation: size at fracture, size of the
#propagule, proportion of cells the propagule is keeping, identity of the cluster before
#it fractured, alive status of the fragmenting edge, and cluster_id/parent_id so lineages
#can be tracked across the branching population
def save_fragmentation_summary(propagule, edge_to_remove, temp_generation_cluster, size_network, curr_identity, alive_cells, cluster_id, parent_id):
	sim_number=input_variables["sim_number"]

	len_propagule=len(propagule)
	prop_propagule=str(len_propagule/size_network)

	dict_files['fragmentation'].write(sim_number+","+str(temp_generation_cluster)+","+str(edge_to_remove[0])+"-"+str(edge_to_remove[1])+","+str(size_network)+","+str(len_propagule)+","+prop_propagule+","+curr_identity+","+alive_cells+","+str(cluster_id)+","+str(parent_id)+"\n")

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

#Function that checks every currently-existing cluster for fragmentation and, for every
#one that is over the edge-degree threshold, fragments it and keeps BOTH resulting
#components growing (ported/merged from sim_frag_clust_edge_degree_all_clusters_exponential_9mar2024.py's
#fragmentation_edge_degree_population, with the death_prob script's alive_cells check and
#diameter-file metrics folded in). A cluster is only removed from cluster_population (and
#from cells_to_divide) once its own lineage reaches num_generations. Also returns a flag
#saying whether at least one cluster fragmented during this call, used by
#simulate_population_growth to decide whether extinction should trigger a retry.
def fragmentation_edge_degree_population(cluster_population, edge_deg_threshold, cont_cluster_ids, num_generations, cells_to_divide):
	#create a copy of all networks/clusters that currently exist
	list_clusters=[cluster_population.subgraph(c).copy() for c in nx.connected_components(cluster_population)]

	any_fragmented=False

	for temp_cluster in list_clusters:
		#calculate the highest edge degree and the edge with the highest edge degree
		edge_to_remove, edge_degree=calculate_max_edge_degree(temp_cluster)

		if(edge_degree>=edge_deg_threshold):
			any_fragmented=True

			size_network=temp_cluster.number_of_nodes()

			list_temp_cluster_nodes=list(temp_cluster.nodes())

			#Check the alive status of the 2 nodes forming the edge that is about to be
			#removed, on the intact (pre-fracture) cluster
			alive_cells=get_alive_status_of_edge(temp_cluster, edge_to_remove)

			#Test fragmentation if enabled
			if(input_variables['test_frag']):
				is_fractured_edge_between_highest_degree_nodes(temp_cluster, edge_to_remove)

			#removing the edge in cluster_population
			cluster_population.remove_edge(*edge_to_remove)

			temp_cluster_id=temp_cluster.nodes[list_temp_cluster_nodes[0]]['cluster_id']

			#update cluster generations in all nodes of the cluster that is fragmenting
			#(both what will become the propagule and the rest of the cluster)
			temp_clust_generations=cluster_population.nodes[list_temp_cluster_nodes[0]]['cluster_generation']
			for temp_node in list_temp_cluster_nodes:
				cluster_population.nodes[temp_node]['cluster_generation']+=1
				cluster_population.nodes[temp_node]['parent_id']=temp_cluster_id

			#identify the propagule (smaller resulting component) and give it a new cluster_id
			fractured_temp_cluster=temp_cluster.copy()
			fractured_temp_cluster.remove_edge(*edge_to_remove)
			propagule = min(nx.connected_components(fractured_temp_cluster), key=len)
			for temp_node in list(propagule):
				cluster_population.nodes[temp_node]['cluster_id']=cont_cluster_ids
			cont_cluster_ids+=1

			#save fragmentation summary information (identity/parent_id read from temp_cluster,
			#i.e. as they were right before this fragmentation event)
			curr_identity=temp_cluster.nodes[list_temp_cluster_nodes[0]]['identity']
			temp_parent_id=temp_cluster.nodes[list_temp_cluster_nodes[0]]['parent_id']

			save_fragmentation_summary(propagule, edge_to_remove, temp_clust_generations, size_network, curr_identity, alive_cells, temp_cluster_id, temp_parent_id)

			#save network diameter/mothers-with-undivided/filamentous-branches on the
			#unbroken (pre-fracture) network
			save_intermediate_networks(temp_cluster, temp_clust_generations, temp_cluster_id)

			#Update identity of the two resulting clusters
			for node in temp_cluster.nodes():
				if node in propagule:
					cluster_population.nodes[node]['identity']='propagule'
				else:
					cluster_population.nodes[node]['identity']='parent'

			#Check whether this (now-split) lineage has reached the target number of
			#generations; if so, prune it entirely from the simulation
			temp_clust_generations+=1 #+1 because we are asking if the upcoming generation should grow or not
			if(temp_clust_generations>=num_generations):
				cells_to_divide=remove_cells_by_component(cells_to_divide, temp_cluster.nodes())
				cluster_population.remove_nodes_from([node for node in cluster_population.nodes() if node in temp_cluster.nodes()])

	return(cont_cluster_ids, cells_to_divide, any_fragmented)

#Main loop where the growth of the whole population (starting from a single cell, and
#branching out into more and more independent clusters via fragmentation) is simulated
#until every cluster currently in existence has either been pruned (reached num_generations)
#or the whole population has died out
def simulate_population_growth(input_variables):
	num_generations=input_variables['num_generations']
	edge_degree_threshold=input_variables['edge_degree_threshold']
	first_div_params=input_variables['first_div_params']
	second_div_params=input_variables['second_div_params']

	cont_ids=1

	#creating graph and adding first cell
	cluster_population=nx.Graph()
	cluster_population.add_node(cont_ids, number_divisions=0, identity='parent',
		cluster_generation=0, cluster_id=1, parent_id=1, alive=True)

	#sample doubling time and initialize cells to divide
	temp_next_doub = sample_doub_t(first_div_params, second_div_params, 0)
	cells_to_divide = []
	cells_to_divide = add_cells_list(temp_next_doub, cont_ids, cells_to_divide)
	cont_ids += 1

	cont_cluster_ids=2

	#Tracks whether at least one fragmentation event has happened anywhere in the
	#population yet. Extinction (cells_to_divide emptying out from cell death) is only
	#treated as a failed attempt if it happens before the first fragmentation - once
	#the population has fragmented at least once, the run is kept as-is even if some
	#later lineage happens to die out early, since there are still plenty of other
	#fragmentation events from the run to study.
	any_fragmentation_happened = False

	while len(cells_to_divide) > 0:
		#dividing all the cells that are going to divide at this timepoint
		cont_ids, cells_to_divide = divide_cells(cluster_population, cells_to_divide, cont_ids)

		#check highest edge degree for every currently-existing cluster, and fragment
		#(keeping both resulting components) any that are over the threshold
		cont_cluster_ids, cells_to_divide, fragmented_this_round = fragmentation_edge_degree_population(
			cluster_population, edge_degree_threshold, cont_cluster_ids, num_generations, cells_to_divide)

		if fragmented_this_round:
			any_fragmentation_happened = True

	#Only a pre-first-fragmentation extinction counts as a failed attempt that needs a
	#retry. Once any_fragmentation_happened is True the attempt is always kept, even if
	#cells_to_divide ended up emptying out because a later lineage died out early.
	reached_target = any_fragmentation_happened

	return reached_target

'''
Steps to find if the nodes of highest degree are forming the edge that is getting fractured:
1.- Find nodes of highest edge degree
2.- Get all edges combinations of the highest degree nodes
3.- Compare all possible pairs to the actual edge deleted
'''
def is_fractured_edge_between_highest_degree_nodes(network, edge_to_remove):
	sim_number=input_variables["sim_number"]

	#1
	highest_degree_nodes=get_nodes_with_highest_degrees(network)

	#2
	list_possible_edges=pairwise_permutations(highest_degree_nodes)

	#3
	temp_flag=str(edge_in_list(edge_to_remove, list_possible_edges))

	#write answer in input file
	dict_files['test_frag'].write(sim_number+','+temp_flag+'\n')

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

#Function to initialize all output files, a dictionary with all the opened files is returned for it
#to be accessed by all functions
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
		dict_files['diameter'].write("sim_number,generation,diameter,cluster_id,cases_mother_with_undivided_cells,num_filamentous_branches\n")

	fragmentation_file=os.path.join(output_dir, "fragmentation_inf.csv")
	if(os.path.exists(fragmentation_file)):
		dict_files['fragmentation']=open(fragmentation_file, "a")
	else:
		dict_files['fragmentation']=open(fragmentation_file, "w")
		dict_files['fragmentation'].write("sim_number,generation,removed_edge,cluster_size,size_propagule,proportion_propagule,identity,alive_cells,cluster_id,parent_id\n")

	sampled_times_file=os.path.join(output_dir,"sampled_times.csv")
	if(os.path.exists(sampled_times_file)):
		dict_files['sampled_times']=open(sampled_times_file, "a")
	else:
		dict_files['sampled_times']=open(sampled_times_file, "w")
		dict_files['sampled_times'].write("sim_number,number_divisions,minutes"+"\n")

	#temporary file
	if(input_variables['test_frag']):
		test_file=os.path.join(output_dir, "test_nodes_frag.csv")
		if(os.path.exists(test_file)):
			dict_files['test_frag']=open(test_file, "a")
		else:
			dict_files['test_frag']=open(test_file, "w")
			dict_files['test_frag'].write("sim_number,flag\n")

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
"first_div_params":first_div_params, "second_div_params":second_div_params,
"death_prob":death_prob}

#Initializing all output files. Note: this variable is accessed as a global variable
dict_files=initialize_files(input_variables)

#these lists track, for every multiple-of-10 generation, whether a graphml network has
#already been saved for it (across the whole run) and, if so, its filepath (so a failed
#attempt can delete it during rollback)
list_flags_saved_network=[False for i in range(num_generations+1) if i%10==0]
list_saved_network_filepaths=[None for i in range(num_generations+1) if i%10==0]

#creating and writing information of the program in the log file
log_file=open(os.path.join(output_dir,"log.txt"), "a")
log_file.write("sim_frag_clust_all_clusters_death_prob_28sep2026.py\nInputs received:\n")
log_file.write("-i: "+args.first_dist_params+"\n")
log_file.write("-j: "+args.second_dist_params+"\n")
log_file.write("-n: "+args.number_sims+"\n")
log_file.write("-o: "+args.output_dir+"\n")
log_file.write("-g: "+args.generations+"\n")
log_file.write("-e: "+args.edge_degree+"\n")
log_file.write("-p: "+args.death_prob+"\n")
log_file.write("\nExecution times:\n")

for i in range(1, num_iterations+1):
	input_variables["sim_number"]=str(i) #changing variable type as it is only used as a character for saving results of the simulation

	#counting execution time
	start_time = time.time()

	#Keep retrying this simulation number until every cluster in the population reaches
	#its target number of generations. Each output file's write position, and which
	#network-saving flags/filepaths are already set, are bookmarked before every attempt
	#so that a failed (extinct) attempt can be fully erased before retrying.
	attempt = 1
	while True:
		file_bookmarks = {}
		for file_name, file_obj in dict_files.items():
			file_obj.flush()
			file_bookmarks[file_name] = file_obj.tell()
		flags_snapshot = list_flags_saved_network.copy()

		#grow the whole population (all clusters) for this simulation
		reached_target = simulate_population_growth(input_variables)

		if reached_target:
			break

		#Extinction before every cluster reached -g generations: erase everything written
		#to disk during this failed attempt
		for file_name, file_obj in dict_files.items():
			file_obj.flush()
			file_obj.seek(file_bookmarks[file_name])
			file_obj.truncate()

		#Remove any graphml network files written during the failed attempt, and reset the
		#corresponding "already saved" flags so the retry can save them
		for idx in range(len(list_flags_saved_network)):
			if (not flags_snapshot[idx]) and list_flags_saved_network[idx]:
				failed_network_file = list_saved_network_filepaths[idx]
				if failed_network_file and os.path.exists(failed_network_file):
					os.remove(failed_network_file)
				list_flags_saved_network[idx] = False
				list_saved_network_filepaths[idx] = None

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