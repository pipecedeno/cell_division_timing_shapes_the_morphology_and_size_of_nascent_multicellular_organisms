#!/usr/bin/env python3

'''
Date: 27may2025

This code simulates snowflake yeast growth without fragmentation, growing clusters until they 
reach a maximum size. For each simulation, it records:
- Network diameter
- Number of mother cells with 2+ undivided daughter cells
- Max edge degree

The first simulation saves all networks after each division event.

Inputs:
	-doubling_t: file of empirical data from where the cell doubling times are going to be extracted from
	-number_sims: number of simulations
	-output_dir: output directory 
	-max_clust_size: maximum number of cells the clusters are allowed to have

Outputs:
	-diameter: columns "sim_number,diameter,cases_mother_with_undivided_cells,max_edge_degree,num_nodes"
	-diff_d_t: columns "sim_number,diff_minutes"
	-graphml files for all networks in first simulation

'''

import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import time

parser = argparse.ArgumentParser()
parser.add_argument('-d','--doubling_t',dest="doubling_t",required=True)
parser.add_argument('-n','--number_sims',dest="number_sims",required=True)
parser.add_argument('-o','--output_dir',dest="output_dir",required=True)
parser.add_argument('-m','--max_clust_size',dest="max_clust_size",required=True)
args = parser.parse_args()

def load_doubling_time(path_to_file, column):
	temp_file = pd.read_csv(path_to_file, header=0)
	grouped_data = temp_file.groupby('division_number')['minutes'].apply(list)
	doub_t_dist = dict(zip(grouped_data.index, grouped_data.values))
	return doub_t_dist

def sample_doub_t(dict_doub_t_dist, divisions_cell, curr_time):
	keys_num_div = dict_doub_t_dist.keys()
	if divisions_cell in dict_doub_t_dist:
		sampling_dist = divisions_cell
	else:
		sampling_dist = max(keys_num_div)
	return curr_time + random.choice(dict_doub_t_dist[sampling_dist])

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

# This function adds all the cells that are going to divide for that specific time
def divide_cells(cluster_population, cells_to_divide, dict_doub_t_dist, cont_ids):
	curr_time = cells_to_divide[0][0]
	list_ids = cells_to_divide[0][1].copy()  # Make a copy to avoid modifying the original list 
	# (this is the nodes that are going to divide this time)
	max_size = input_variables['max_clust_size']
	reached_max_size = False

	for mother_id in list_ids:
		# Check if adding a new cell would exceed max size
		if cluster_population.number_of_nodes() >= max_size:
			reached_max_size = True
			break

		# Create daughter cell
		daughter_id = cont_ids
		cluster_population.add_node(daughter_id, number_divisions=0)
		cluster_population.add_edge(mother_id, daughter_id)
		temp_daughter_time = sample_doub_t(dict_doub_t_dist, 0, curr_time)
		
		# Only add daughter's next division time if we haven't reached max size
		if cluster_population.number_of_nodes() < max_size:
			cells_to_divide = add_cells_list(temp_daughter_time, daughter_id, cells_to_divide)
		cont_ids += 1

		# Divide mother cell
		cluster_population.nodes[mother_id]["number_divisions"] += 1
		temp_mother_div = cluster_population.nodes[mother_id]["number_divisions"]
		temp_mother_time = sample_doub_t(dict_doub_t_dist, temp_mother_div, curr_time)
		
		# Only add mother's next division time if we haven't reached max size
		if cluster_population.number_of_nodes() < max_size:
			cells_to_divide = add_cells_list(temp_mother_time, mother_id, cells_to_divide)

		# Save doubling time difference
		dict_files['diff_d_t'].write(input_variables["sim_number"] + "," + 
								   str(abs(temp_mother_time-temp_daughter_time)) + "\n")
		
		# Check if we've reached max size after adding cells
		if cluster_population.number_of_nodes() >= max_size:
			reached_max_size = True
			break

	cells_to_divide.pop(0)
	
	# If we reached max size, clear all future divisions
	if reached_max_size:
		cells_to_divide.clear()
		
	return [cont_ids, cells_to_divide]

def count_mothers_with_multiple_undivided(network):
	mothers_with_multiple_undivided = 0
	for node in network.nodes():
		neighbors = list(network.neighbors(node))
		undivided_daughters = sum(1 for neighbor in neighbors if network.degree(neighbor) == 1)
		if undivided_daughters >= 2:
			mothers_with_multiple_undivided += 1
	return mothers_with_multiple_undivided

# Function to calculate and return the edge with the highest edge degree
def calculate_max_edge_degree(network):
	# Initialize result variables
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

def save_network_parameters(network, sim_number, division_step=None):
	"""Save network diameter, mothers with undivided cells, number of nodes, and max edge degree"""
	if network.number_of_nodes() > 1:  # Need at least 2 nodes for diameter
		network_diameter = nx.diameter(network)
	else:
		network_diameter = 0
	
	mothers_with_undivided = count_mothers_with_multiple_undivided(network)
	num_nodes = network.number_of_nodes()
	
	# Calculate max edge degree
	if network.number_of_edges() > 0:  # Need at least 1 edge for max edge degree
		max_edge, max_edge_degree = calculate_max_edge_degree(network)
	else:
		max_edge_degree = -1
	
	# Save parameters
	dict_files['diameter'].write(f"{sim_number},{network_diameter},{mothers_with_undivided},{num_nodes},{max_edge_degree}\n")

# function to save the graphml files of the network for each step only for the first simulation
def save_network_if_needed(network, sim_number, division_step=None):
	"""Save network as graphml if it's the first simulation"""
	if int(sim_number) == 1:
		if division_step is not None:
			num_nodes = network.number_of_nodes()
			network_path = os.path.join(input_variables['network_dir'], f'network_sim1_step{division_step}_n{num_nodes}.graphml')
		else:
			network_path = os.path.join(input_variables['network_dir'], f'network_sim1_final.graphml')
		nx.write_graphml(network, network_path)

def simulate_one_cluster_growth(input_variables, dict_doub_t_dist):
	max_clust_size = input_variables['max_clust_size']
	sim_number = input_variables['sim_number']
	
	# Initialize cluster with first cell
	cont_ids = 1
	snowflake = nx.Graph()
	snowflake.add_node(cont_ids, number_divisions=0)
	
	# Save initial network parameters (1 node)
	save_network_parameters(snowflake, sim_number)
	save_network_if_needed(snowflake, sim_number, 0)
	
	# Sample first doubling time
	temp_next_doub = sample_doub_t(dict_doub_t_dist, 0, 0)
	cells_to_divide = []
	cells_to_divide = add_cells_list(temp_next_doub, cont_ids, cells_to_divide)
	cont_ids += 1

	curr_time = cells_to_divide[0][0]
	division_step = 1
	
	# Grow cluster until max size
	while len(cells_to_divide) > 0 and snowflake.number_of_nodes() < max_clust_size:
		# Update current time
		curr_time = cells_to_divide[0][0]
		
		# Divide cells at current time
		cont_ids, cells_to_divide = divide_cells(snowflake, cells_to_divide, dict_doub_t_dist, cont_ids)
		
		# Save network parameters after each division event
		save_network_parameters(snowflake, sim_number)
		
		# Save network if it's the first simulation
		save_network_if_needed(snowflake, sim_number, division_step)
		
		division_step += 1
		
		# Check if we've reached max size
		if snowflake.number_of_nodes() >= max_clust_size:
			break

def initialize_files(input_variables):
	output_dir = input_variables['output_dir']
	dict_files = {}
	
	# Initialize diameter file
	diam_file = os.path.join(output_dir, "network_information.csv")
	if os.path.exists(diam_file):
		dict_files['diameter'] = open(diam_file, "a")
	else:
		dict_files['diameter'] = open(diam_file, "w")
		dict_files['diameter'].write("sim_number,diameter,cases_mother_with_undivided_cells,num_nodes,max_edge_degree\n")
	
	# Initialize doubling time difference file
	diff_d_t_file = os.path.join(output_dir, "diff_doub_t.csv")
	if os.path.exists(diff_d_t_file):
		dict_files['diff_d_t'] = open(diff_d_t_file, "a")
	else:
		dict_files['diff_d_t'] = open(diff_d_t_file, "w")
		dict_files['diff_d_t'].write("sim_number,diff_minutes\n")
	
	return dict_files

def close_all_files(dict_files):
	for file_obj in dict_files.values():
		file_obj.close()

#### MAIN ####

# Load doubling time distributions
dict_doub_t_dist = load_doubling_time(args.doubling_t, "minutes")
output_dir = args.output_dir
max_clust_size = int(args.max_clust_size)
num_iterations = int(args.number_sims)

# Create output directory
if not os.path.exists(output_dir):
	os.mkdir(output_dir)
else:
	print("Output directory already exists, stopping program execution to avoid overwriting data")
	exit()

# Create network directory
network_dir = os.path.join(output_dir, 'network_dir')
if not os.path.exists(network_dir):
	os.mkdir(network_dir)

# Initialize input variables dictionary
input_variables = {
	"dict_doub_t_dist": dict_doub_t_dist,
	"output_dir": output_dir,
	"network_dir": network_dir,
	"max_clust_size": max_clust_size
}

# Initialize output files
dict_files = initialize_files(input_variables)

# Create and write to log file
log_file = open(os.path.join(output_dir, "log.txt"), "a")
log_file.write("sim_clust_no_fragmentation_27may2024.py\nInputs received:\n")
log_file.write("-d: " + args.doubling_t + "\n")
log_file.write("-n: " + args.number_sims + "\n")
log_file.write("-o: " + args.output_dir + "\n")
log_file.write("-m: " + args.max_clust_size + "\n")
log_file.write("\nExecution times:\n")

# Run simulations
for i in range(1, num_iterations + 1):
	input_variables["sim_number"] = str(i)
	
	start_time = time.time()
	simulate_one_cluster_growth(input_variables, dict_doub_t_dist)
	end_time = time.time()
	elapsed_time = end_time - start_time
	
	log_file.write(f"Simulation {i} completed in {elapsed_time} seconds.\n")
	print(f"{i}/{num_iterations}", end='\r')

log_file.write("\nProgram finished execution without any errors.")
log_file.close()

# Close all output files
close_all_files(dict_files)