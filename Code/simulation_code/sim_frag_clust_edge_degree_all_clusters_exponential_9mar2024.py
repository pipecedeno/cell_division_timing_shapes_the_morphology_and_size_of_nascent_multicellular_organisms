#!/usr/bin/env python3

'''
Date: 22feb2024

For this program fragmentation is going to happen when the highest edge degree
of the cluster gets higher than the threshold. And all clusters are going to allow to be
grown until the amount of cluster generations is reached. 

This code is a modification of sim_frag_clust_edge_degree_20oct2023.py.

###
Notes:
One of the biases that this program has when there are many cells that divide at the same time,
is that the first cell that appeared in the list with that time is the one that is going to divide
first, and also because of how the elements are added to the list, if the mother and the daughter
cell sample the same doubling time, the daughter is going to divide first because that is the cell
that gets added first to the list.
###

Inputs:
	-e: edge degree threshold of when fragmentation should happen
	-d: doubling time distributions, csv file that has the doubling times to be sampled from
	-n: number of the simulation/network, is to name the files differently when saving them.
	-o: output directory
	-g: number of generations to end the simulation

Outputs:
	-degree distribution file (saved as probabilities)
	-difference doubling time between mother a daughter cell
	-fragmentation information (edges removed, size at fracture, proportion of the propagule)
	-networks diameter
	-Network files are only saved for the first network to reach multiples of 10

Modification history:
Date: 1 March 2024
Changed the code so that now the propagule that is going to be deleted can either be the parent, propagule or
selected at random. Added a new column to the fragmentation file so that it can be known which is going to be
the next propagule in case it is of interest to filter those values when doing the parent offspring regressions.
Also the code for cell divisions and cell list is the same as in growth_phase_clusters_v0_9feb2024.py, which is
the code of the growth phase of the settling selection.
Note: the divide_cells was slightly changed and also the sample_doub_t as there are not 2 strains from which to
sample doubling times.

Date: 4March2024
Updated the code so that now it saves the identity of the cluster in the fragmentation file and the nodes are
going to store a new property called identity, which can be either 'parent' or 'propagule'. This will be useful
when analyzing the heritability of the network properties, as I could filter and only keep the prent-offspring
regressions, which are all were the offspring has 'propagule' in the identity column.

Date: 13November2024
Added the function count_mothers_with_multiple_undivided which allows to count how many mother cells have two or
more undivided daughter cells, this informations is saved in the fragmentation file in the column:
cases_mother_with_undivided_cells


'''

import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import random
# import netlsd
# import numpy as np
# from scipy.stats import skew, kurtosis
import time

#only used for testing purposes
from itertools import permutations

parser = argparse.ArgumentParser()
parser.add_argument('-d','--doubling_t',dest="doubling_t",required=True) #csv file of the doubling time distributions
parser.add_argument('-n','--sim_number',dest="sim_number",required=True) #number of simulations it is going to run
parser.add_argument('-o','--output_dir',dest="output_dir",required=True) #output directory
parser.add_argument('-g','--generations',dest="generations",required=True) #Number of generations
parser.add_argument('-e','--edge_degree',dest="edge_degree",required=True) #edge degree threshold
# parser.add_argument('-c','--cluster_to_keep',dest="cluster_to_keep",required=True) # 'parent', 'propagule', 'random'
args = parser.parse_args()

### Global variables ####
# list_flags_saved_network=[]




#### Functions for simulating networks growth ########


#This function returns a dictionary with where each key is the number of division and the value associated to that key
#is the observed doubling times
def load_doubling_time(path_to_file, column):
	temp_file=pd.read_csv(path_to_file, header=0)

	# Group the data by 'division_number' and extract the 'minutes' values into a list
	grouped_data = temp_file.groupby('division_number')['minutes'].apply(list)

	# Create a dictionary where the key is the unique value of 'division_number' and the value is the list of 'minutes'
	doub_t_dist = dict(zip(grouped_data.index, grouped_data.values))

	return(doub_t_dist)


#sampling the doubling time of the next division
def sample_doub_t(dict_doub_t_dist, divisions_cell, curr_time):
	keys_num_div=dict_doub_t_dist.keys()

	#selecting from which distribution the answer is going to be sampled
	#if the number of times the cell has divided is in the dictionary then that distribution is going
	#to be used for the sampling of the doubling time, if it is not in the dictionary then the last
	#distribution is going to be used
	if(divisions_cell in dict_doub_t_dist):
		sampling_dist=divisions_cell
	else:
		sampling_dist=max(keys_num_div)

	#using a uniform distribution to sample from the list
	return(curr_time+random.choice(dict_doub_t_dist[sampling_dist]))



#This function adds the new doubling time to the list of doubling times, because it is expected
#that the new doubling time is bigger than the elements of the list, it is going to iterate it
#in reverse to make this more efficient. The element input is formed of [time until doubling, id 
#of cell], and ordered list is the already ordered list of doublings.
# def add_to_ordered_list(element, ordered_list):
# 	# Find the index to insert the element by iterating in reverse order
# 	for i, el in reversed(list(enumerate(ordered_list))):
# 		if element[0] >= el[0]:
# 			ordered_list.insert(i + 1, element)
# 			return(ordered_list)
# 	# If the element's number is smaller than all elements, insert it at the beginning
# 	ordered_list.insert(0, element)
# 	return(ordered_list)

# #This function is to substract the first value to the whole list, it will return the new updated
# #list of values
# def subtract_time(number, lst):
# 	modified_lst = []
# 	for element in lst:
# 		time = element[0] - number
# 		modified_element = [time] + element[1:]
# 		modified_lst.append(modified_element)
# 	return(modified_lst)

'''
input_time: doubling time of the cell
cell_id: id of the dividing cell
ordered_list: list of the doubling times, normally called cells_to_divide

cells_to_divide entries:
-time next division
-id of the cell
(I don't need to know the id of the cluster as they all are going to be in the same nx.graph() )

cells_to_divide=[[55, [1, 3, 4]], [60, [5, 7]], [65, [6]]]
'''
def add_cells_list(input_time, cell_id, ordered_list):
	# If the list is empty, insert new element at the beginning
	if not ordered_list:
		ordered_list.append([input_time, [cell_id]])
		return ordered_list
	
	# Find were to add the new time
	for i, el in reversed(list(enumerate(ordered_list))):
		#if the time already exists then insert the cell_id in the list of ids
		if (input_time == el[0]):
			ordered_list[i][-1].append(cell_id)
			return(ordered_list)
		
		#if the time is bigger than one of the elements then insert the new time and the id 
		#after the current position
		elif (input_time>el[0]):
			ordered_list.insert(i+1, [input_time, [cell_id]])
			return(ordered_list)
	
	# If the input's number is smaller than all elements, insert it at the beginning
	ordered_list.insert(0, [input_time, [cell_id]])
	return(ordered_list)


#Function to save the degree distribution of the network, the degree distribution is saved as a
#probability distribution to be normalized for network size
def save_degree_distribution_to_csv(network, num_generations, temp_cluster_id):
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
	degree_counts['cluster_id']=temp_cluster_id

	for i in range(len(degree_counts)):
		dict_files['degree_dist'].write(",".join([str(j) for j in degree_counts.iloc[i]])+"\n")


#Function to remove the next doubling times of the cells of the parent propagule that is being
#discarded
# def remove_cells_by_component(list_cells, component):
# 	positions_to_remove = []
	
# 	# Iterate through the list of lists and save the positions to remove
# 	for i, cell_data in enumerate(list_cells):
# 		if(cell_data[1] in component):
# 			positions_to_remove.append(i)
	
# 	# Remove the entries from the list of lists in reverse order to avoid index issues
# 	for pos in reversed(positions_to_remove):
# 		list_cells.pop(pos)
	
# 	return(list_cells)

def remove_cells_by_component(cells_to_divide, component_to_remove):
	# Initialize the updated list for storing cells to divide without the removed component
	updated_list = []

	# Iterate through each time entry and the corresponding cells scheduled to divide
	for time_entry, cells in cells_to_divide:
		# Filter out cells that are in the component to remove
		filtered_cells = [cell for cell in cells if cell not in component_to_remove]

		# Only add to the updated list if there are cells left to divide at this time
		if(filtered_cells):
			updated_list.append([time_entry, filtered_cells])

	return(updated_list)

#function to save intermediate networks
def save_intermediate_networks(network, temp_generation_cluster):

	# output_dir=input_variables["output_dir"]
	sim_number=input_variables["sim_number"]
	# network_dir=input_variables["network_dir"]

	# global list_flags_saved_network

	# #save cluster before if it is the first cluster of this generation
	# if(temp_generation_cluster%10==0):
	# 	if(list_flags_saved_network[int(temp_generation_cluster/10)]==False):
	# 		nx.write_graphml(network, os.path.join(network_dir,str(sim_number)+"_"+str(temp_generation_cluster)+".graphml"))
			
	# 		list_flags_saved_network[int(temp_generation_cluster/10)]=True

	list_temp_cluster_nodes=list(network.nodes())

	#calculate network diameter for each network and append it into the file
	temp_diameter_network=str(nx.diameter(network))
	temp_cluster_id=network.nodes[list_temp_cluster_nodes[0]]['cluster_id']
	dict_files['diameter'].write(sim_number+","+str(temp_generation_cluster)+","+temp_diameter_network+','+str(temp_cluster_id)+"\n")

	save_degree_distribution_to_csv(network, temp_generation_cluster, temp_cluster_id)


def count_mothers_with_multiple_undivided(network):
	"""
	Count the number of mother cells that have 2 or more undivided daughter cells
	in the network.
	
	Parameters:
	network (networkx.Graph): The network to analyze
	
	Returns:
	int: Number of mother cells with 2 or more undivided daughter cells
	"""
	mothers_with_multiple_undivided = 0
	
	for node in network.nodes():
		# Get all neighbors (daughters) of the current node
		neighbors = list(network.neighbors(node))
		
		# Count how many of these neighbors have degree 1 (undivided cells)
		undivided_daughters = sum(1 for neighbor in neighbors if network.degree(neighbor) == 1)
		
		# If the node has 2 or more undivided daughters, increment the counter
		if undivided_daughters >= 2:
			mothers_with_multiple_undivided += 1
			
	return(mothers_with_multiple_undivided)

#function to save information about the fragmentation, size at fracture, size of the parent after
#fracture and the percentage/proportion of cells the parent is keeping
def save_fragmentation_summary(propagule, propagule_network, edge_to_remove, temp_generation_cluster, size_network, curr_identity, cluster_id, temp_parent_id):
	# output_dir=input_variables["output_dir"]
	sim_number=input_variables["sim_number"]

	len_propagule=len(propagule)
	prop_propagule=str(len_propagule/size_network)

	cont_mothers_with_undivided=count_mothers_with_multiple_undivided(propagule_network)

	dict_files['fragmentation'].write(str(sim_number)+","+str(temp_generation_cluster)+","+str(edge_to_remove[0])+"-"+str(edge_to_remove[1])+","+str(size_network)+","+str(len_propagule)+","+str(prop_propagule)+","+str(cont_mothers_with_undivided)+","+curr_identity+','+str(cluster_id)+','+str(temp_parent_id)+"\n")


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


def fragmentation_edge_degree_population(cluster_population, edge_deg_threshold, cont_cluster_ids, num_generations, cells_to_divide):
	# create a copy of all networks 
	list_clusters=[cluster_population.subgraph(c).copy() for c in nx.connected_components(cluster_population)]

	for temp_cluster in list_clusters:
		#calculate the highest edge degree and the edge with the highest edge degree
		edge_to_remove, edge_degree=calculate_max_edge_degree(temp_cluster)

		if(edge_degree>=edge_deg_threshold):

			size_network=temp_cluster.number_of_nodes()
			cluster_population.remove_edge(*edge_to_remove)
			
			list_temp_cluster_nodes=list(temp_cluster.nodes())
			temp_cluster_id=temp_cluster.nodes[list_temp_cluster_nodes[0]]['cluster_id'] #To update who is the parent of the clusters

			#update cluster generations in all nodes of the cluster that is fracturing (both in the propagule and in the parent)
			temp_clust_generations=cluster_population.nodes[list_temp_cluster_nodes[0]]['cluster_generation']
			for temp_node in list_temp_cluster_nodes:
				cluster_population.nodes[temp_node]['cluster_generation']+=1
				#make the parent_id the cluster_id for all cells
				cluster_population.nodes[temp_node]['parent_id']=temp_cluster_id

			#update cluster ids for the propagule
			#use temp cluster to find the nodes of the propagule after fracture
			fractured_temp_cluster=temp_cluster.copy()
			fractured_temp_cluster.remove_edge(*edge_to_remove)
			propagule = min(nx.connected_components(fractured_temp_cluster), key=len)
			for temp_node in list(propagule):
				cluster_population.nodes[temp_node]['cluster_id']=cont_cluster_ids
			
			cont_cluster_ids+=1


			#save fragmentation summary information
			curr_identity=temp_cluster.nodes[list_temp_cluster_nodes[0]]['identity']
			temp_cluster_id=temp_cluster.nodes[list_temp_cluster_nodes[0]]['cluster_id']
			temp_parent_id=temp_cluster.nodes[list_temp_cluster_nodes[0]]['parent_id'] #to save the information of who the parent of the cluster is

			propagule_network=temp_cluster.subgraph(propagule).copy()

			save_fragmentation_summary(propagule, propagule_network, edge_to_remove, temp_clust_generations, size_network, curr_identity, temp_cluster_id, temp_parent_id)

			#Removing cells from the selected component
			# cells_to_divide=remove_cells_by_component(cells_to_divide, component_to_remove)
			# cluster_population.remove_nodes_from([node for node in cluster_population.nodes() if node in component_to_remove])

			#save network diameter and degree distribution
			save_intermediate_networks(temp_cluster, temp_clust_generations)

			#Update indentity of the clusters
			for node in temp_cluster.nodes():
				if node in propagule:
					cluster_population.nodes[node]['identity']='propagule'
				else:
					cluster_population.nodes[node]['identity']='parent'


			#Update cluster generations
			temp_clust_generations+=1 #+1 because we are asking if the upcoming generation should grow or not
			if(temp_clust_generations>=num_generations):
				cells_to_divide=remove_cells_by_component(cells_to_divide, temp_cluster.nodes())
				cluster_population.remove_nodes_from([node for node in cluster_population.nodes() if node in temp_cluster.nodes()])



	return(cont_cluster_ids, cells_to_divide)



'''
In the function divide cells all the cells that are dividing in the current time point are going to divided
-cluster_population is the networkx Graph object with the networks of all the clusters in the population
-cells_to_divide is the list from where the current time and the list of cells that are going to be dividing this time
	is going to be extracted from
Note: the function doesn't return the cluster_population variable as it receives a pointer to the object and the real
variable is being modified
'''
def divide_cells(cluster_population, cells_to_divide, dict_doub_t_dist, cont_ids, cont_gen):

	sim_number=input_variables['sim_number']

	curr_time=cells_to_divide[0][0]
	list_ids=cells_to_divide[0][1] #ids of dividing cells

	cont_mother_id_bigger=0

	for mother_id in list_ids:

		temp_mother_identity=cluster_population.nodes[mother_id]['identity']
		mother_clust_gen=cluster_population.nodes[mother_id]['cluster_generation']
		mother_clust_id=cluster_population.nodes[mother_id]['cluster_id']
		temp_parent_id=cluster_population.nodes[mother_id]['parent_id']

		#create daughter cell
		daughter_id=cont_ids
		cluster_population.add_node(cont_ids, number_divisions=0, identity=temp_mother_identity, cluster_generation=mother_clust_gen, cluster_id=mother_clust_id, parent_id=temp_parent_id)
		cluster_population.add_edge(mother_id, daughter_id)
		temp_daughter_time=sample_doub_t(dict_doub_t_dist, 0, curr_time) #using 0 instead of cluster_population.nodes[daughter_id]["number_divisions"]
		cells_to_divide=add_cells_list(temp_daughter_time, daughter_id, cells_to_divide)
		cont_ids+=1

		if(daughter_id<mother_id):
			cont_mother_id_bigger+=1

		#divide mother cell
		cluster_population.nodes[mother_id]["number_divisions"]+=1 #adding that the mother divided one more time (because of the step before)
		temp_mother_div=cluster_population.nodes[mother_id]["number_divisions"]
		temp_mother_time=sample_doub_t(dict_doub_t_dist, temp_mother_div, curr_time)
		cells_to_divide=add_cells_list(temp_mother_time, mother_id, cells_to_divide)

		dict_files['diff_d_t'].write(sim_number+","+str(cont_gen)+","+str(abs(temp_mother_time-temp_daughter_time))+"\n")

	if(cont_mother_id_bigger>0):
		print('Mother id is bigger, '+str(cont_mother_id_bigger))

	#to get rid of the current time used
	cells_to_divide.pop(0)

	return([cont_ids, cells_to_divide])


#Main loop where the growth of one cluster is simulated until it reaches the desired amount of
#generations
def simulate_one_cluster_growth(input_variables, dict_doub_t_dist):

	#defining general use variables
	output_dir=input_variables['output_dir']
	num_generations=input_variables['num_generations']
	edge_degree_threshold=input_variables['edge_degree_threshold']

	cont_ids=1

	#creating graph and adding first cell
	cluster_population=nx.Graph() #creating graph

	cluster_population.add_node(cont_ids,number_divisions=0, identity='parent', cluster_generation=0, cluster_id=1, parent_id=1) #adding the first cell
	
	#sample doubling time
	temp_next_doub=sample_doub_t(dict_doub_t_dist, 0, 0)

	#add value to the list in order
	cells_to_divide=[]
	cells_to_divide=add_cells_list(temp_next_doub, cont_ids, cells_to_divide)

	cont_ids+=1

	t_sim=0

	cont_updates=1

	cluster_fract_size=-1

	cont_gen=0

	cont_cluster_ids=2

	reached_max_generations=False

	#update time so that the next round cells can divide
	# curr_time=cells_to_divide[0][0]

	# while reached_max_generations==False:
	while len(cells_to_divide)>0:

		# curr_time=cells_to_divide[0][0]
		
		#dividing all the cells that are going to divide this timepoint
		cont_ids, cells_to_divide=divide_cells(cluster_population, cells_to_divide, dict_doub_t_dist, cont_ids, cont_gen) #it gets a reference to the network variable so it doesn't need to return the variable
		#note: save doubling times need to be added to divide cells

		#check highest edge degree for all clusters to trigger fragmentation if any goes beyond the threshold
		cont_cluster_ids, cells_to_divide=fragmentation_edge_degree_population(cluster_population, edge_degree_threshold, cont_cluster_ids, num_generations, cells_to_divide)

		cont_updates+=1

		#delete the times in cells to divide of the current time point
		# cells_to_divide=update_cells_to_divide(cells_to_divide)
		# cells_to_divide.pop(0)

		#update current time and population size
		# curr_time=get_smallest_time(cells_to_divide)
		# last_time=curr_time
		# curr_time=cells_to_divide[0][0] #select the next smallest time


	#save final networks and their information
	# save_final_networks(snowflake, cont_gen, cont_updates)


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
		dict_files['diameter'].write("sim_number,generation,diameter,cluster_id\n")

	fragmentation_file=os.path.join(output_dir, "fragmentation_inf.csv")
	if(os.path.exists(fragmentation_file)):
		dict_files['fragmentation']=open(fragmentation_file, "a")
	else:
		dict_files['fragmentation']=open(fragmentation_file, "w")
		dict_files['fragmentation'].write("sim_number,generation,removed_edge,cluster_size,size_propagule,proportion_propagule,cases_mother_with_undivided_cells,identity,cluster_id,parent_id\n")

	degree_dist_file=os.path.join(output_dir,"degree_distribution.csv")
	if(os.path.exists(degree_dist_file)):
		dict_files['degree_dist']=open(degree_dist_file, "a")
	else:
		dict_files['degree_dist']=open(degree_dist_file, "w")
		dict_files['degree_dist'].write("Degree,Frequency,Probability,num_generations,sim_number,cluster_id\n")

	# updates_file=os.path.join(output_dir, "updates_to_fragmentation.csv")
	# if(os.path.exists(updates_file)):
	# 	dict_files['snowflakes_update']=open(updates_file, "a")
	# else:
	# 	dict_files['snowflakes_update']=open(updates_file, "w")
	# 	dict_files['snowflakes_update'].write("sim_number,generation,updates\n")

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

'''
To do:


-Correct the file of fragmentation to save the identity of the cluster as 'parent' or 'propagule',
	This new property is going to be saved in the nodes as identity and it is going to be changed depending
	on which cluster is selected to continue growing.
	The initial value for the cluster is 'parent'

-Make sure that I actually got rid off all network similarity stuff

-Get rid off the cluster_to_keep stuff

Copy the loop from the growth phase simulations:
	-The loop needs to be until there is no more entries in cells_to_divide
	-divide cells function
	-fragmentation_edge_degree_population
		I can make it so that I break the cluster but after saving their fragmentation information I
		am going to check if they have reached the maximum amount of cluster generations, and if that
		is the case I can remove those nodes from the snowflake graph variable and the times from the
		cells_to_divide list, that way I won't need to use the save_final_networks function and just use
		save_intermediate_networks all the time. And with this I should also get information from one last 
		doubling which I wasn't able to get before
	-Keep other parts of the loop 



From sim_frag_pop_edge_degree_5oct2023.py:
-Copy how they check the number of generations after each fragmentation
	Add to get rid off the clusters that have finished growing and won't continue growing
	Repurpose the function remove_cells_by_component which now needs to be used after
	the cluster has reached the maximum amount of generations.


-Check the ram usage of running 8 generations and 9, and see if it is possible to run it in parallel
	Maybe to run it in parallel I need to keep the results in RAM but they are going to be to many
	threads writing in the same output file.
	Maybe the best option of parallelization is to do it in a bash script


Note:
-Check cluster generations in the population simulations just to get an idea of how many cluster 
generations are happening

'''

#initializing variables
#Loading doubling time distributions
dict_doub_t_dist=load_doubling_time(args.doubling_t, "minutes")

output_dir=args.output_dir

num_generations=int(args.generations)

# num_iterations=int(args.number_sims)
sim_number=args.sim_number

edge_degree_threshold=int(args.edge_degree)

# cluster_to_keep=args.cluster_to_keep.lower()
# if(not cluster_to_keep in ['parent', 'propagule', 'random']):
# 	print("Error no valid option given in argument -c\nAll valid options are 'parent', 'offspring', 'random'")
# 	exit()

#creating output directory if it doesn't exist, and stopping the program execution if the 
#directory exists
# if(not os.path.exists(output_dir)):
# 	os.mkdir(output_dir)
# else:
# 	print("Output directory already exists, stopping program execution to avoid overwriting data")
# 	exit()
# network_dir=os.path.join(output_dir, 'network_dir')
# if(not os.path.exists(network_dir)):
# 	os.mkdir(network_dir)


#For testing if the edge of the fragmentation is formed by the nodes of highest edge degree
test_edge_of_fragmentation=True


#Creating temporary directory for the output of this program
simulation_directory=os.path.join(output_dir, 'sim_'+sim_number)
os.mkdir(simulation_directory)
#this directory is going to be used as the output directory


#Dictionary containing input variables and parameters of the simulation
#this variable is used as a global variable
input_variables={"dict_doub_t_dist":dict_doub_t_dist,
	"output_dir":simulation_directory,
	"num_generations":num_generations, 
	"edge_degree_threshold":edge_degree_threshold,
	"test_frag":test_edge_of_fragmentation,
	"sim_number":sim_number}
#"cluster_to_keep":cluster_to_keep
#"network_dir":network_dir

#Initializing all output files. Note: this variable is accesed as a global variable
dict_files=initialize_files(input_variables)



#this list is to save if the network has already been saved in the network files,
#only the first network and the final networks are going to be saved
#note: this variable is used as a global variable to avoid having to pass it in all functions
#and return it after being modified
#file name (sim_number)_(generation).graphml
# list_flags_saved_network=[False for i in range(num_generations+1) if i%10==0]


simulate_one_cluster_growth(input_variables, dict_doub_t_dist)


### Loop before making the code in parallel
#creating and writing information of the program in the log file
# log_file=open(os.path.join(output_dir,"log.txt"), "a")
# log_file.write("sim_frag_clust_edge_degree_20oct2023.py\nInputs received:\n")
# log_file.write("-d: "+args.doubling_t+"\n")
# log_file.write("-n: "+args.number_sims+"\n")
# log_file.write("-o: "+args.output_dir+"\n")
# log_file.write("-g: "+args.generations+"\n")
# log_file.write("-e: "+args.edge_degree+"\n")
# log_file.write("\nExecution times:\n")
# for i in range(1, num_iterations+1):
# 	input_variables["sim_number"]=str(i) #changing variable type as it is only used as a character for saving results of the simulation
# 	#counting execution time
# 	start_time = time.time()
# 	#grow the cluster of this simulation
# 	simulate_one_cluster_growth(input_variables, dict_doub_t_dist)
# 	end_time=time.time()
# 	elapsed_time=end_time-start_time
# 	log_file.write("Simulation "+str(i)+" completed in "+str(elapsed_time)+" seconds.\n")
# 	print(str(i)+'/'+str(num_iterations), end='\r')
# log_file.write("\nProgram finished execution without any errors.")
# log_file.close()

#closing all files from the dictionary
close_all_files(dict_files)



