#!/usr/bin/env python3

'''
Date: 22feb2024

For this program fragmentation is going to happen when the highest edge degree
of the cluster gets higher than the threshold. The code will grow only one cluster at a time,
depending on the input of the user.



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
	-feature values file (netsimile feature vectors)
	-fragmentation information (edges removed, size at fracture, proportion of the biggest fragment)
	-heat feature values (netlsd feature vector)
	-networks diameter
	-updates to fragmentation, it saves how many cells divided (updates) between fragmentations
	-Network files are only saved for the first network to reach multiples of 10


'''

import os
import pandas as pd
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import random
import netlsd
import numpy as np
from scipy.stats import skew, kurtosis
import time

#only used for testing purposes
from itertools import permutations

parser = argparse.ArgumentParser()
parser.add_argument('-d','--doubling_t',dest="doubling_t",required=True) #csv file of the doubling time distributions
parser.add_argument('-n','--number_sims',dest="number_sims",required=True) #number of simulations it is going to run
parser.add_argument('-o','--output_dir',dest="output_dir",required=True) #output directory
parser.add_argument('-g','--generations',dest="generations",required=True) #Number of generations
parser.add_argument('-e','--edge_degree',dest="edge_degree",required=True) #edge degree threshold
parser.add_argument('-c','--cluster_to_keep',dest="cluster_to_keep",required=True) # 'parent', 'propagule', 'random'
args = parser.parse_args()

### Global variables ####
list_flags_saved_network=[]


####### Functions to calculate NetSimile ########


#Function to compute feature extraction from netsimile library in netrd
#This function calculates the distributions for each network from each of the properties
#analyzed by netsimile
def feature_extraction(G):
	"""Node feature extraction.

	Parameters
	----------

	G (nx.Graph): a networkx graph.

	Returns
	-------

	node_features (float): the Nx7 matrix of node features."""

	# necessary data structures
	node_features = np.zeros(shape=(G.number_of_nodes(), 7))
	node_list = sorted(G.nodes())
	node_degree_dict = dict(G.degree())
	node_clustering_dict = dict(nx.clustering(G))
	egonets = {n: nx.ego_graph(G, n) for n in node_list}

	# node degrees
	degs = [node_degree_dict[n] for n in node_list]

	# clustering coefficient
	clusts = [node_clustering_dict[n] for n in node_list]

	# average degree of neighborhood
	neighbor_degs = [
		np.mean([node_degree_dict[m] for m in egonets[n].nodes if m != n])
		if node_degree_dict[n] > 0
		else 0
		for n in node_list
	]

	# average clustering coefficient of neighborhood
	neighbor_clusts = [
		np.mean([node_clustering_dict[m] for m in egonets[n].nodes if m != n])
		if node_degree_dict[n] > 0
		else 0
		for n in node_list
	]

	# number of edges in the neighborhood
	neighbor_edges = [
		egonets[n].number_of_edges() if node_degree_dict[n] > 0 else 0
		for n in node_list
	]

	# number of outgoing edges from the neighborhood
	# the sum of neighborhood degrees = 2*(internal edges) + external edges
	# node_features[:,5] = node_features[:,0] * node_features[:,2] - 2*node_features[:,4]
	neighbor_outgoing_edges = [
		len(
			[
				edge
				for edge in set.union(*[set(G.edges(j)) for j in egonets[i].nodes])
				if not egonets[i].has_edge(*edge)
			]
		)
		for i in node_list
	]

	# number of neighbors of neighbors (not in neighborhood)
	neighbors_of_neighbors = [
		len(
			set([p for m in G.neighbors(n) for p in G.neighbors(m)])
			- set(G.neighbors(n))
			- set([n])
		)
		if node_degree_dict[n] > 0
		else 0
		for n in node_list
	]

	# assembling the features
	node_features[:, 0] = degs
	node_features[:, 1] = clusts
	node_features[:, 2] = neighbor_degs
	node_features[:, 3] = neighbor_clusts
	node_features[:, 4] = neighbor_edges
	node_features[:, 5] = neighbor_outgoing_edges
	node_features[:, 6] = neighbors_of_neighbors

	return np.nan_to_num(node_features)


#This function calculates the distribution summary statistics for eahc of the properties tracked 
#for each network
def graph_signature(node_features):
	signature_vec = np.zeros(7 * 5)

	# for each of the 7 features
	for k in range(7):
		# find the mean
		signature_vec[k * 5] = node_features[:, k].mean()
		# find the median
		signature_vec[k * 5 + 1] = np.median(node_features[:, k])
		# find the std
		signature_vec[k * 5 + 2] = node_features[:, k].std()
		# find the skew
		signature_vec[k * 5 + 3] = skew(node_features[:, k])
		# find the kurtosis
		signature_vec[k * 5 + 4] = kurtosis(node_features[:, k])

	return(signature_vec)







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
	dict_files['diameter'].write(sim_number+","+str(temp_generation_cluster)+","+temp_diameter_network+"\n")

	save_degree_distribution_to_csv(network, temp_generation_cluster)

	#calculate parent feature values of netlsd and netsimile
	if(input_variables['compute_similarity']):
		temp_heat=netlsd.heat(network)
		#saving heat values of network
		dict_files['heat_val'].write(sim_number+","+str(temp_generation_cluster)+","+",".join([str(i) for i in temp_heat])+"\n")

		temp_node_features=feature_extraction(network)
		temp_feat_vec=graph_signature(temp_node_features)
		#saving feature values
		dict_files['feature_val'].write(sim_number+","+str(temp_generation_cluster)+","+",".join([str(i) for i in temp_feat_vec])+"\n")



#function to save information about the fragmentation, size at fracture, size of the parent after
#fracture and the percentage/proportion of cells the parent is keeping
def save_fragmentation_summary(propagule, edge_to_remove, temp_generation_cluster, size_network, curr_identity):
	output_dir=input_variables["output_dir"]
	sim_number=input_variables["sim_number"]

	len_propagule=len(propagule)
	prop_propagule=str(len_propagule/size_network)

	dict_files['fragmentation'].write(sim_number+","+str(temp_generation_cluster)+","+str(edge_to_remove[0])+"-"+str(edge_to_remove[1])+","+str(size_network)+","+str(len_propagule)+","+prop_propagule+","+curr_identity+"\n")


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
	save_fragmentation_summary(propagule, edge_to_remove, temp_generation_cluster, size_network, curr_identity)

	#Removing cells from the selected component
	cells_to_divide=remove_cells_by_component(cells_to_divide, component_to_remove)
	snowflake.remove_nodes_from([node for node in snowflake.nodes() if node in component_to_remove])

	#save network properties before fracture
	save_intermediate_networks(unbroken_cluster, temp_generation_cluster)

	#Save how many updates have passed
	dict_files['snowflakes_update'].write(sim_number+","+str(temp_generation_cluster)+","+str(cont_updates)+"\n")

	if(temp_generation_cluster>max_num_generations):
		cells_to_divide=remove_cells_by_component(cells_to_divide, snowflake.nodes())

	return(cells_to_divide)


#Function to save properties of the final network
def save_final_networks(network, cont_gen, cont_updates):
	network_dir=input_variables["network_dir"]
	sim_number=input_variables["sim_number"]


	global list_flags_saved_network

	#save cluster before if it is the first cluster of this generation
	if(cont_gen%10==0):
		if(list_flags_saved_network[int(cont_gen/10)]==False):
			nx.write_graphml(network, os.path.join(network_dir,str(sim_number)+"_"+str(cont_gen)+".graphml"))
			
			list_flags_saved_network[int(cont_gen/10)]=True

	temp_diameter_network=str(nx.diameter(network))
	dict_files['diameter'].write(sim_number+","+str(cont_gen)+","+temp_diameter_network+"\n")

	save_degree_distribution_to_csv(network, cont_gen)

	#calculate parent feature values of netlsd and netsimile
	if(input_variables['compute_similarity']):
		temp_heat=netlsd.heat(network)
		#saving heat values of network
		dict_files['heat_val'].write(sim_number+","+str(cont_gen)+","+",".join([str(i) for i in temp_heat])+"\n")

		temp_node_features=feature_extraction(network)
		temp_feat_vec=graph_signature(temp_node_features)
		#saving feature values
		dict_files['feature_val'].write(sim_number+","+str(cont_gen)+","+",".join([str(i) for i in temp_feat_vec])+"\n")

	#Save how many updates have passed
	dict_files['snowflakes_update'].write(sim_number+","+str(cont_gen)+","+str(cont_updates)+"\n")


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

		#create daughter cell
		daughter_id=cont_ids
		cluster_population.add_node(cont_ids, number_divisions=0, identity=temp_mother_identity)
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
	snowflake=nx.Graph() #creating graph

	snowflake.add_node(cont_ids,number_divisions=0, identity='parent') #adding the first cell
	
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

	reached_max_generations=False

	#update time so that the next round cells can divide
	curr_time=cells_to_divide[0][0]

	# while reached_max_generations==False:
	while len(cells_to_divide)>0:

		#update current time and population size
		last_time=curr_time
		curr_time=cells_to_divide[0][0] #select the next smallest time

		#dividing all the cells that are going to divide this timepoint
		cont_ids, cells_to_divide=divide_cells(snowflake, cells_to_divide, dict_doub_t_dist, cont_ids, cont_gen) #it gets a reference to the network variable so it doesn't need to return the variable
		#note: save doubling times need to be added to divide cells


		# calculate first order edge degree
		edge_to_remove, edge_degree=calculate_max_edge_degree(snowflake)

		if(edge_degree>=edge_degree_threshold):

			#find if the edge to remove is formed by the highest degree nodes
			if(input_variables['test_frag']):
				is_fractured_edge_between_highest_degree_nodes(snowflake, edge_to_remove, cont_gen)

			#Calling function to break clusters, update list of cells_to_divide, and save intermediate networks
			cells_to_divide=break_cluster(snowflake, cells_to_divide, cont_gen, cont_updates, edge_to_remove, num_generations)

			cont_gen+=1


			# if(cont_gen<num_generations):

			# 	#find if the edge to remove is formed by the highest degree nodes
			# 	if(input_variables['test_frag']):
			# 		is_fractured_edge_between_highest_degree_nodes(snowflake, edge_to_remove, cont_gen)

			# 	#Calling function to break clusters, update list of cells_to_divide, and save intermediate networks
			# 	cells_to_divide=break_cluster(snowflake, cells_to_divide, cont_gen, cont_updates, edge_to_remove)

			# 	cont_gen+=1

			# else:
				
			# 	reached_max_generations=True

		cont_updates+=1

		#delete the times in cells to divide of the current time point
		# cells_to_divide=update_cells_to_divide(cells_to_divide)
		# cells_to_divide.pop(0)



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
		dict_files['diameter'].write("sim_number,generation,diameter\n")

	fragmentation_file=os.path.join(output_dir, "fragmentation_inf.csv")
	if(os.path.exists(fragmentation_file)):
		dict_files['fragmentation']=open(fragmentation_file, "a")
	else:
		dict_files['fragmentation']=open(fragmentation_file, "w")
		dict_files['fragmentation'].write("sim_number,generation,removed_edge,cluster_size,size_propagule,proportion_propagule,identity\n")

	degree_dist_file=os.path.join(output_dir,"degree_distribution.csv")
	if(os.path.exists(degree_dist_file)):
		dict_files['degree_dist']=open(degree_dist_file, "a")
	else:
		dict_files['degree_dist']=open(degree_dist_file, "w")
		dict_files['degree_dist'].write("Degree,Frequency,Probability,num_generations,sim_number\n")

	if(input_variables['compute_similarity']):
		heat_values_file=os.path.join(output_dir,"heat_values.csv")
		if(os.path.exists(heat_values_file)):
			dict_files['heat_val']=open(heat_values_file, "a")
		else:
			dict_files['heat_val']=open(heat_values_file, "w")
			dict_files['heat_val'].write("sim_number,generation,"+",".join(["val_"+str(i) for i in range(1, 251)])+"\n")

		feat_values_file=os.path.join(output_dir,"feature_values.csv")
		if(os.path.exists(feat_values_file)):
			dict_files['feature_val']=open(feat_values_file, "a")
		else:
			dict_files['feature_val']=open(feat_values_file, "w")
			dict_files['feature_val'].write("sim_number,generation,"+",".join(["val_"+str(i) for i in range(1, 36)])+"\n")

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
#Loading doubling time distributions
dict_doub_t_dist=load_doubling_time(args.doubling_t, "minutes")

output_dir=args.output_dir

num_generations=int(args.generations)

num_iterations=int(args.number_sims)

edge_degree_threshold=int(args.edge_degree)

cluster_to_keep=args.cluster_to_keep.lower()
if(not cluster_to_keep in ['parent', 'propagule', 'random']):
	print("Error no valid option given in argument -c\nAll valid options are 'parent', 'offspring', 'random'")
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
test_edge_of_fragmentation=True


#This flag will be used to turn off network similarity calculations as it is no longer necesary to obtain them
compute_similarity=False

#Dictionary containing input variables and parameters of the simulation
#this variable is used as a global variable
input_variables={"dict_doub_t_dist":dict_doub_t_dist,
	"output_dir":output_dir,
	"network_dir":network_dir,
	"num_generations":num_generations, 
	"edge_degree_threshold":edge_degree_threshold,
	"test_frag":test_edge_of_fragmentation,
	"compute_similarity":compute_similarity,
	"cluster_to_keep":cluster_to_keep}

#Initializing all output files. Note: this variable is accesed as a global variable
dict_files=initialize_files(input_variables)



#this list is to save if the network has already been saved in the network files,
#only the first network and the final networks are going to be saved
#note: this variable is used as a global variable to avoid having to pass it in all functions
#and return it after being modified
#file name (sim_number)_(generation).graphml
list_flags_saved_network=[False for i in range(num_generations+1) if i%10==0]

#creating and writing information of the program in the log file
log_file=open(os.path.join(output_dir,"log.txt"), "a")
log_file.write("sim_frag_clust_edge_degree_20oct2023.py\nInputs received:\n")
log_file.write("-d: "+args.doubling_t+"\n")
log_file.write("-n: "+args.number_sims+"\n")
log_file.write("-o: "+args.output_dir+"\n")
log_file.write("-g: "+args.generations+"\n")
log_file.write("-e: "+args.edge_degree+"\n")
log_file.write("-c: "+args.cluster_to_keep+"\n")
log_file.write("\nExecution times:\n")

for i in range(1, num_iterations+1):

	input_variables["sim_number"]=str(i) #changing variable type as it is only used as a character for saving results of the simulation

	#counting execution time
	start_time = time.time()

	#grow the cluster of this simulation
	simulate_one_cluster_growth(input_variables, dict_doub_t_dist)

	end_time=time.time()
	elapsed_time=end_time-start_time


	log_file.write("Simulation "+str(i)+" completed in "+str(elapsed_time)+" seconds.\n")

	print(str(i)+'/'+str(num_iterations), end='\r')


log_file.write("\nProgram finished execution without any errors.")
log_file.close()

#closing all files from the dictionary
close_all_files(dict_files)




