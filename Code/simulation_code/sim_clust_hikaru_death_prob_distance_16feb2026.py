#!/usr/bin/env python3

'''
Date: 19aug2025

This code simulates snowflake yeast growth without fragmentation, growing clusters until they 
reach a maximum size. For each simulation, it records:
- Network diameter
- Number of mother cells with 2+ undivided daughter cells
- Max edge degree
- Number of filamentous branches
- Distance from basal cell (node 1) for all nodes
- Edge degree distribution (node degree frequencies)

The first simulation saves the network after completion.

Inputs:
	-i: doubling_time_means: comma-separated mean doubling times for normal distribution (e.g., "40,90,60")
	    First value = first division, second value = second division, last value = third and subsequent divisions
	-j: doubling_time_sd: standard deviation for normal distribution (will be scaled for each mean)
	-p: death_prob: probability of cell death after each division (0-1)
	-number_sims: number of simulations
	-output_dir: output directory 
	-max_clust_size: maximum number of cells the clusters are allowed to have

Outputs:
	-network_information.csv: columns "sim_number,diameter,cases_mother_with_undivided_cells,max_edge_degree,num_nodes,num_filamentous_branches"
	-distance_basal_cell.csv: columns "sim_number,distance_basal,count"
	-edge_degree_dist.csv: columns "sim_number,degree,frequency"
	-graphml files for final networks in first simulation

Modifications:
29may2025:
Made the simulation to save the network properties every 5 cells added instead of doing it every time all the cells off a specific time
are added.

Modified to use log-normal distribution parameters instead of empirical data file and added filamentous branches counting.

5Aug2025:
Added cell death probability parameter and safety feature to restart simulations that don't reach target size.
Removed diff_doub_t.csv output file.

19aug2025:
- Added distance_basal_cell.csv output with distance analysis from basal cell (node 1)
- Changed to only save network information and files at final network size
- Removed every-5-nodes saving behavior

2sep2025:
- Added edge_degree_dist.csv output with node degree distribution analysis

16feb2026:
- Changed input parameters to accept normal distribution parameters (mean, sd) which are converted to lognormal
- Changed to support multiple doubling time distributions (not limited to 2)
- Standard deviation is scaled by (mean_i / mean_last) for each distribution

'''

import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import random
import numpy as np
import time
import math

parser = argparse.ArgumentParser()
parser.add_argument('-i','--doubling_time_means',dest="doubling_time_means",required=True, 
                    help="Comma-separated mean doubling times (e.g., '40,90,60')")
parser.add_argument('-j','--doubling_time_sd',dest="doubling_time_sd",required=True,
                    help="Standard deviation for doubling time normal distribution")
parser.add_argument('-p','--death_prob',dest="death_prob",required=True)
parser.add_argument('-n','--number_sims',dest="number_sims",required=True)
parser.add_argument('-o','--output_dir',dest="output_dir",required=True)
parser.add_argument('-m','--max_clust_size',dest="max_clust_size",required=True)
args = parser.parse_args()

def normal_to_lognormal(mean, sd):
	"""
	Convert normal distribution parameters to lognormal distribution parameters.
	
	Based on the R code:
	variance <- sd^2
	meanlog <- log(mean^2 / sqrt(variance + mean^2))
	sdlog <- sqrt(log(1 + variance / mean^2))
	
	Args:
		mean: Mean of the normal distribution
		sd: Standard deviation of the normal distribution
	
	Returns:
		tuple: (meanlog, sdlog) - parameters for lognormal distribution
	"""
	variance = sd ** 2
	meanlog = math.log(mean ** 2 / math.sqrt(variance + mean ** 2))
	sdlog = math.sqrt(math.log(1 + variance / mean ** 2))
	
	return meanlog, sdlog

def sample_doub_t(lognormal_params_list, divisions_cell, curr_time):
	"""Sample doubling time from log-normal distribution based on division number"""
	# Determine which distribution to use based on division number
	if divisions_cell < len(lognormal_params_list):
		dist_mu, dist_sigma = lognormal_params_list[divisions_cell]
	else:
		# Use the last distribution for all subsequent divisions
		dist_mu, dist_sigma = lognormal_params_list[-1]
	
	doubling_time = round(float(np.random.lognormal(dist_mu, dist_sigma, 1)[0]), 4)
	return curr_time + doubling_time

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

# Function to identify filament nodes
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

def count_filamentous_branches(graph, min_length=3):
	"""Count the number of filamentous branches in the network"""
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

def calculate_distance_from_basal_cell(network, basal_cell_id=1):
	"""Calculate distance from basal cell (node 1) to all other nodes and return counts by distance"""
	if basal_cell_id not in network.nodes():
		return {}
	
	# Calculate shortest path lengths from basal cell to all other nodes
	distances = nx.single_source_shortest_path_length(network, basal_cell_id)
	
	# Count nodes at each distance
	distance_counts = {}
	for node, distance in distances.items():
		if distance not in distance_counts:
			distance_counts[distance] = 0
		distance_counts[distance] += 1
	
	return distance_counts

def calculate_edge_degree_distribution(network):
	"""Calculate the degree distribution of the network"""
	degree_counts = {}
	
	# Get degree for each node
	for node in network.nodes():
		degree = network.degree(node)
		if degree not in degree_counts:
			degree_counts[degree] = 0
		degree_counts[degree] += 1
	
	return degree_counts

# This function adds all the cells that are going to divide for that specific time
def divide_cells(cluster_population, cells_to_divide, lognormal_params_list, death_prob, cont_ids):
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
		temp_daughter_time = sample_doub_t(lognormal_params_list, 0, curr_time)
		
		# Check if daughter cell survives (cell death probability check)
		daughter_survives = random.random() >= death_prob
		
		# Only add daughter's next division time if it survives and we haven't reached max size
		if daughter_survives and cluster_population.number_of_nodes() < max_size:
			cells_to_divide = add_cells_list(temp_daughter_time, daughter_id, cells_to_divide)
		cont_ids += 1

		# Divide mother cell
		cluster_population.nodes[mother_id]["number_divisions"] += 1
		temp_mother_div = cluster_population.nodes[mother_id]["number_divisions"]
		temp_mother_time = sample_doub_t(lognormal_params_list, temp_mother_div, curr_time)
		
		# Check if mother cell survives
		mother_survives = random.random() >= death_prob
		
		# Only add mother's next division time if it survives and we haven't reached max size
		if mother_survives and cluster_population.number_of_nodes() < max_size:
			cells_to_divide = add_cells_list(temp_mother_time, mother_id, cells_to_divide)

	cells_to_divide.pop(0)
	return cont_ids, cells_to_divide

def save_network_parameters(network, sim_number, results_table):
	"""Save network parameters to memory table"""
	diameter = nx.diameter(network)
	
	# Count mothers with 2+ undivided daughters
	cases_mother_with_undivided_cells = 0
	for node in network.nodes():
		neighbors = list(network.neighbors(node))
		undivided_daughters = sum(1 for n in neighbors if network.nodes[n]["number_divisions"] == 0)
		if undivided_daughters >= 2:
			cases_mother_with_undivided_cells += 1
	
	# Get max edge degree
	degrees = [network.degree(node) for node in network.nodes()]
	max_edge_degree = max(degrees) if degrees else 0
	
	# Count filamentous branches
	num_filamentous_branches = count_filamentous_branches(network, min_length=3)
	
	num_nodes = network.number_of_nodes()
	
	# Store in memory table instead of writing directly to file
	results_table.append([sim_number, diameter, cases_mother_with_undivided_cells, num_nodes, max_edge_degree, num_filamentous_branches])

def save_distance_basal_cell(network, sim_number, distance_results_table):
	"""Save distance from basal cell analysis to memory"""
	distance_counts = calculate_distance_from_basal_cell(network, basal_cell_id=1)
	
	# Store each distance and its count
	for distance, count in distance_counts.items():
		distance_results_table.append([sim_number, distance, count])

def save_edge_degree_distribution(network, sim_number, degree_results_table):
	"""Save edge degree distribution analysis to memory"""
	degree_counts = calculate_edge_degree_distribution(network)
	
	# Store each degree and its frequency
	for degree, frequency in degree_counts.items():
		degree_results_table.append([sim_number, degree, frequency])

# function to save the graphml files of the network for the first simulation
def save_network_if_needed(network, sim_number, network_files):
	"""Store network as graphml data in memory if it's the first simulation"""
	if int(sim_number) == 1:
		num_nodes = network.number_of_nodes()
		filename = f'network_sim1_final_n{num_nodes}.graphml'
		
		# Store network data in memory using BytesIO instead of StringIO
		import io
		buffer = io.BytesIO()
		nx.write_graphml(network, buffer)
		network_files[filename] = buffer.getvalue().decode('utf-8')  # Convert bytes to string for storage
		buffer.close()

def simulate_one_cluster_growth(input_variables, lognormal_params_list, death_prob):
	max_clust_size = input_variables['max_clust_size']
	sim_number = input_variables['sim_number']
	
	# Initialize memory storage for results
	results_table = []
	distance_results_table = []
	degree_results_table = []
	network_files = {}
	
	# Initialize cluster with first cell
	cont_ids = 1
	snowflake = nx.Graph()
	snowflake.add_node(cont_ids, number_divisions=0)
	
	# Sample first doubling time
	temp_next_doub = sample_doub_t(lognormal_params_list, 0, 0)
	cells_to_divide = []
	cells_to_divide = add_cells_list(temp_next_doub, cont_ids, cells_to_divide)
	cont_ids += 1

	curr_time = cells_to_divide[0][0]
	
	# Grow cluster until max size
	while len(cells_to_divide) > 0 and snowflake.number_of_nodes() < max_clust_size:
		# Update current time
		curr_time = cells_to_divide[0][0]
		
		# Divide cells at current time
		cont_ids, cells_to_divide = divide_cells(snowflake, cells_to_divide, lognormal_params_list, death_prob, cont_ids)
		
		# Check if we've reached max size
		if snowflake.number_of_nodes() >= max_clust_size:
			break
	
	# Only save network parameters and files when final size is reached
	if snowflake.number_of_nodes() >= max_clust_size:
		# Save final network parameters to memory table
		save_network_parameters(snowflake, sim_number, results_table)
		
		# Save distance from basal cell analysis
		save_distance_basal_cell(snowflake, sim_number, distance_results_table)
		
		# Save edge degree distribution analysis
		save_edge_degree_distribution(snowflake, sim_number, degree_results_table)
		
		# Save network if it's the first simulation
		save_network_if_needed(snowflake, sim_number, network_files)
	
	# Return whether the simulation reached the target size and the results
	success = snowflake.number_of_nodes() >= max_clust_size
	return success, results_table, distance_results_table, degree_results_table, network_files

def initialize_files(input_variables):
	output_dir = input_variables['output_dir']
	dict_files = {}
	
	# Initialize network information file
	network_info_file = os.path.join(output_dir, "network_information.csv")
	if os.path.exists(network_info_file):
		dict_files['network_info'] = open(network_info_file, "a")
	else:
		dict_files['network_info'] = open(network_info_file, "w")
		dict_files['network_info'].write("sim_number,diameter,cases_mother_with_undivided_cells,num_nodes,max_edge_degree,num_filamentous_branches\n")
	
	# Initialize distance basal cell file
	distance_file = os.path.join(output_dir, "distance_basal_cell.csv")
	if os.path.exists(distance_file):
		dict_files['distance_basal'] = open(distance_file, "a")
	else:
		dict_files['distance_basal'] = open(distance_file, "w")
		dict_files['distance_basal'].write("sim_number,distance_basal,count\n")
	
	# Initialize edge degree distribution file
	degree_file = os.path.join(output_dir, "edge_degree_dist.csv")
	if os.path.exists(degree_file):
		dict_files['edge_degree_dist'] = open(degree_file, "a")
	else:
		dict_files['edge_degree_dist'] = open(degree_file, "w")
		dict_files['edge_degree_dist'].write("sim_number,degree,frequency\n")
	
	return dict_files

def write_results_to_disk(results_table, distance_results_table, degree_results_table, network_files, dict_files, input_variables):
	"""Write stored results to disk files"""
	# Write network information to CSV
	for row in results_table:
		dict_files['network_info'].write(f"{row[0]},{row[1]},{row[2]},{row[3]},{row[4]},{row[5]}\n")
	
	# Write distance basal cell information to CSV
	for row in distance_results_table:
		dict_files['distance_basal'].write(f"{row[0]},{row[1]},{row[2]}\n")
	
	# Write edge degree distribution information to CSV
	for row in degree_results_table:
		dict_files['edge_degree_dist'].write(f"{row[0]},{row[1]},{row[2]}\n")
	
	# Write network files to disk (only for first simulation)
	for filename, content in network_files.items():
		network_path = os.path.join(input_variables['network_dir'], filename)
		with open(network_path, 'w', encoding='utf-8') as f:
			f.write(content)

def close_all_files(dict_files):
	for file_obj in dict_files.values():
		file_obj.close()

#### MAIN ####

# Parse doubling time means (comma-separated)
mean_values = [float(elem) for elem in args.doubling_time_means.split(",")]
base_sd = float(args.doubling_time_sd)
death_prob = float(args.death_prob)

# Validate that we have at least one mean value
if len(mean_values) == 0:
	print("Error: At least one doubling time mean must be provided")
	exit()

# Calculate scaled standard deviations
# SD is scaled by (mean_i / mean_last) for each distribution
mean_last = mean_values[-1]
scaled_sds = [base_sd * (mean_i / mean_last) for mean_i in mean_values]

# Convert normal distribution parameters to lognormal parameters
lognormal_params_list = []
for mean_val, sd_val in zip(mean_values, scaled_sds):
	meanlog, sdlog = normal_to_lognormal(mean_val, sd_val)
	lognormal_params_list.append((meanlog, sdlog))

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
	"output_dir": output_dir,
	"network_dir": network_dir,
	"max_clust_size": max_clust_size
}

# Initialize output files
dict_files = initialize_files(input_variables)

# Create and write to log file
log_file = open(os.path.join(output_dir, "log.txt"), "a")
log_file.write("sim_clust_no_fragmentation_lognormal_death_prob_5aug2025_updated.py\nInputs received:\n")
log_file.write("-i (means): " + args.doubling_time_means + "\n")
log_file.write("-j (base SD): " + args.doubling_time_sd + "\n")
log_file.write("-p: " + args.death_prob + "\n")
log_file.write("-n: " + args.number_sims + "\n")
log_file.write("-o: " + args.output_dir + "\n")
log_file.write("-m: " + args.max_clust_size + "\n")
log_file.write("\nProcessed parameters:\n")
log_file.write(f"Normal distribution means: {mean_values}\n")
log_file.write(f"Base standard deviation: {base_sd}\n")
log_file.write(f"Scaled standard deviations: {scaled_sds}\n")
log_file.write(f"Lognormal parameters (meanlog, sdlog): {lognormal_params_list}\n")
log_file.write("\nExecution times:\n")

# Run simulations
successful_simulations = 0
for i in range(1, num_iterations + 1):
	input_variables["sim_number"] = str(i)
	
	start_time = time.time()
	
	# Keep trying until simulation reaches target size
	attempt = 1
	while True:
		success, results_table, distance_results_table, degree_results_table, network_files = simulate_one_cluster_growth(input_variables, lognormal_params_list, death_prob)
		if success:
			# Write results to disk only if simulation was successful
			write_results_to_disk(results_table, distance_results_table, degree_results_table, network_files, dict_files, input_variables)
			successful_simulations += 1
			break
		else:
			attempt += 1
			# No need to reset files - just try the simulation again
	
	end_time = time.time()
	elapsed_time = end_time - start_time
	
	if attempt > 1:
		log_file.write(f"Simulation {i} completed in {elapsed_time} seconds after {attempt} attempts.\n")
	else:
		log_file.write(f"Simulation {i} completed in {elapsed_time} seconds.\n")
	print(f"{i}/{num_iterations}", end='\r')

log_file.write(f"\nProgram finished execution without any errors. {successful_simulations}/{num_iterations} simulations completed successfully.\n")
log_file.close()

# Close all output files
close_all_files(dict_files)