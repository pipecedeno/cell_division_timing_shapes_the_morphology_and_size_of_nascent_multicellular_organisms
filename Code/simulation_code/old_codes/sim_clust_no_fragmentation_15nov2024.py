#!/usr/bin/env python3

'''
Date: 15nov2024

This code simulates snowflake yeast growth without fragmentation, growing clusters until they 
reach a maximum size. For each simulation, it records:
- Network diameter
- Degree distribution
- Number of mother cells with 2+ undivided daughter cells

The first 5 networks are saved as graphml files.

Inputs:
	-doubling_t: file of empirical data from where the cell doubling times are going to be extracted from
	-number_sims: number of simulations
	-output_dir: output directory 
	-max_clust_size: maximum number of cells the clusters are allowed to have

Outputs:
	-degree_dist file: columns "Degree,Frequency,Probability,sim_number"
	-diameter: columns "sim_number,diameter,cases_mother_with_undivided_cells"
	-diff_d_t: columns "sim_number,diff_minutes"
	-graphml files for first 5 networks

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

# Number of networks to save (hardcoded)
NETWORKS_TO_SAVE = 5

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

def save_degree_distribution_to_csv(network):
	# Get the degree of all nodes
	degrees = dict(network.degree())
	df_degrees = pd.DataFrame(degrees.items(), columns=['Node', 'Degree'])
	
	# Create degree distribution
	degree_counts = df_degrees['Degree'].value_counts().reset_index()
	degree_counts.columns = ['Degree', 'Frequency']
	
	# Calculate probability distribution
	total_nodes = len(df_degrees)
	degree_counts['Probability'] = degree_counts['Frequency'] / total_nodes
	
	# Add simulation number
	degree_counts['sim_number'] = input_variables['sim_number']

	for i in range(len(degree_counts)):
		dict_files['degree_dist'].write(",".join([str(j) for j in degree_counts.iloc[i]]) + "\n")

def count_mothers_with_multiple_undivided(network):
	mothers_with_multiple_undivided = 0
	for node in network.nodes():
		neighbors = list(network.neighbors(node))
		undivided_daughters = sum(1 for neighbor in neighbors if network.degree(neighbor) == 1)
		if undivided_daughters >= 2:
			mothers_with_multiple_undivided += 1
	return mothers_with_multiple_undivided

def save_network_if_needed(network, sim_number):
	"""Save network as graphml if it's one of the first NETWORKS_TO_SAVE simulations"""
	if int(sim_number) <= NETWORKS_TO_SAVE:
		network_path = os.path.join(input_variables['network_dir'], f'network_{sim_number}.graphml')
		nx.write_graphml(network, network_path)

def simulate_one_cluster_growth(input_variables, dict_doub_t_dist):
	max_clust_size = input_variables['max_clust_size']
	sim_number = input_variables['sim_number']
	
	# Initialize cluster with first cell
	cont_ids = 1
	snowflake = nx.Graph()
	snowflake.add_node(cont_ids, number_divisions=0)
	
	# Sample first doubling time
	temp_next_doub = sample_doub_t(dict_doub_t_dist, 0, 0)
	cells_to_divide = []
	cells_to_divide = add_cells_list(temp_next_doub, cont_ids, cells_to_divide)
	cont_ids += 1

	curr_time = cells_to_divide[0][0]
	
	# Grow cluster until max size
	while len(cells_to_divide) > 0 and snowflake.number_of_nodes() < max_clust_size:
		# Update current time
		curr_time = cells_to_divide[0][0]
		
		# Divide cells at current time
		cont_ids, cells_to_divide = divide_cells(snowflake, cells_to_divide, dict_doub_t_dist, cont_ids)
		
		if snowflake.number_of_nodes() >= max_clust_size:
			# Save final network properties
			network_diameter = nx.diameter(snowflake)
			mothers_with_undivided = count_mothers_with_multiple_undivided(snowflake)
			
			# Save diameter and mothers with undivided cells
			dict_files['diameter'].write(f"{sim_number},{network_diameter},{mothers_with_undivided}\n")
			
			# Save degree distribution
			# save_degree_distribution_to_csv(snowflake)
			
			# Save network if it's one of the first 5
			save_network_if_needed(snowflake, sim_number)
			break

def initialize_files(input_variables):
	output_dir = input_variables['output_dir']
	dict_files = {}
	
	# Initialize degree distribution file
	# degree_dist_file = os.path.join(output_dir, "degree_distribution.csv")
	# if os.path.exists(degree_dist_file):
	# 	dict_files['degree_dist'] = open(degree_dist_file, "a")
	# else:
	# 	dict_files['degree_dist'] = open(degree_dist_file, "w")
	# 	dict_files['degree_dist'].write("Degree,Frequency,Probability,sim_number\n")
	
	# Initialize diameter file
	diam_file = os.path.join(output_dir, "networks_diameter.csv")
	if os.path.exists(diam_file):
		dict_files['diameter'] = open(diam_file, "a")
	else:
		dict_files['diameter'] = open(diam_file, "w")
		dict_files['diameter'].write("sim_number,diameter,cases_mother_with_undivided_cells\n")
	
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
log_file.write("sim_clust_no_fragmentation_15nov2024.py\nInputs received:\n")
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



