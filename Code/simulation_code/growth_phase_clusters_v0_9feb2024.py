#!/usr/bin/env python3

'''

Date: 9Feb2024
This code is part of the population simulations used to calculate the selection rate of the growth phase and the 
settling phase as if it was a competition experiment in the lab.

inputs:
-doubling time strain 1
-doubling time strain 2
-starting concentrations (format strain1:strain2, example: 100:100)
-temp_dir: directory where the files of the edges are going to be saved and the ids of the clusters that survive settling selection
-output_dir: Here the files time_registry and selection_proportions
-edge degree threshold
-Current selectling selection round (just to know if it is necessary to initialize the files or not, or
	if it needs to seed the population, or load it from a file)

Output:
-final_population_#: csv file with the edges, the time they were added and the id of the cluster
To have the time of when a edge was added I need to keep track of the daughters of each cell
and save the time when they were added (so a list of 2 entries [daughter_cell_id, time])
Cluster id needs to be updated every time a cluster fractures
	strain, node1, node2, time_grown, cluster_id, cluster_generation,

-time_registry: csv file with population at each time point
	sim_number, strain, num_clusters, time, transfer, num_cells, total_clusters, total_cells

-selection_proportions: csv file of proportions of clusters after the growth and settling round
	sim_number, transfer, experiment_phase, end, strain1, strain2, prop_cluster_pop1, prop_cluster_pop2, prop_cell_pop1, prop_cell_pop2,
	total_clusters, total_cells 
	experiment_phase [growth/settling]
	end [True/False] Is it the end of the phase? or the beginning?


Things to keep track:
-Cluster generations (this can only be tracked in the nodes)
-amount of cell per time of the experiment


Parameters to define:
-population capacity
-duration of the growth phase (24 hours)


Instead of iterating through the list of cells, iterate through time
each cell needs to have as a parameter from which strain they are, and that value needs to be used to access
a dictionary with both dictionaries of doubling times

Main loop
1.- while curr_pop<pop_capacity and curr_t<24 hours
2.- divide all cells in that time point
3.- check max edge degree of each network and break the edges that go above the threshold
4.- add 5 minutes to timer and modify the cells_to_divide list

cells_to_divide entries:
-time next division
-id of the cell
(I don't need to know the id of the cluster as they all are going to be in the same nx.graph() )

cells_to_divide=[[55, [1, 3, 4]], [60, [5, 7]], [65, [6]]]

Information stored in the nodes:
	-number_divisions
	-cluster generation
	-strain (this is important to know from which distribution to sample)
	-time added
	-cluster_id (this is important to when loading the edges information to matlab)
	-transfer added
	Note: transfer added is necessary to be able to build the clusters correctly in matlab as the time is not enough
	because time is restarted each new growth phase

Question:
-should cluster generation start from 1 or 0?

Possible output:
At the end of the simulation save the strain and the cluster generation of each cluster,
this may be interesting as the clusters that fracture faster may be going through more generations


Modification history:
Added carrying capacity as an input

Date: 16Apr2024
Tried the option of using a named tuple instead of using the dictionary to use less RAM in the 
information stored in the nodes

'''

import os
import pandas as pd
import argparse
import networkx as nx
import random
import numpy as np

import time
from collections import namedtuple

import psutil
import sys



parser = argparse.ArgumentParser()
parser.add_argument('-f','--strain1_dt',dest="strain1_dt",required=True) #csv file of the doubling time distributions for strain 1
parser.add_argument('-s','--strain2_dt',dest="strain2_dt",required=True) #csv file of the doubling time distributions for strain 2
parser.add_argument('-n','--strain_names',dest="strain_names",required=True) #name of the strains (format: strain1:strain2)
parser.add_argument('-t','--temp_dir',dest="temp_dir",required=True) #temp directory
parser.add_argument('-o','--output_dir',dest="output_dir",required=True) #output directory
parser.add_argument('-e','--edge_degree',dest="edge_degree",required=True) #edge degree threshold
parser.add_argument('-r','--curr_selec_round',dest="curr_selec_round",required=True) #current selection round
parser.add_argument('-c','--pop_concentration',dest="pop_concentration",required=True) #starting concentration of the populations (format example: 100:100)
parser.add_argument('-i','--sim_number',dest="sim_number",required=True) #simulation number
parser.add_argument('-k','--carrying_cap',dest="carrying_cap",required=True) #carrying capacity
args = parser.parse_args()


#### FUNCTIONS #####

#This function returns a dictionary with where each key is the number of division and the value associated to that key
#is the observed doubling times
def load_doubling_time(path_to_file, column):
	temp_file=pd.read_csv(path_to_file, header=0)

	# Group the data by 'division_number' and extract the 'minutes' values into a list
	grouped_data = temp_file.groupby('division_number')['minutes'].apply(list)

	# Create a dictionary where the key is the unique value of 'division_number' and the value is the list of 'minutes'
	doub_t_dist = dict(zip(grouped_data.index, grouped_data.values))

	return(doub_t_dist)


'''
sampling the doubling time of the next division
-dict_doub_t_strains: should be a dictionaty with the entries of both strains of the simulations,
	and each entry of the strain should have a dictionary with the doubling times of the first and second
	doubling at least and a list from where to sample a time
-divisions_cell: how many times the cell for which the time is going to be sampled has divided
-curr_time: current time of the simulation, just to add it to the sampled time
-strain_name: strain of the cell that is going to divide

'''
def sample_doub_t(dict_doub_t_strains, divisions_cell, curr_time, strain_name):

	#selecting the distribution for one of the strains
	dict_doub_t_dist=dict_doub_t_strains[strain_name]

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




'''
This code is where the clusters are going to be loaded from the previously grown clusters to start the simulation
the Graph object is going to be initialized to save the networks and the list of cells_to_divide is also going
to be created.
Here the ids of both strains which were previously grown are non overlapping and that is why we don't need to modify 
them at all to avoid any errors of repeating ids
Also the ids of nodes and clusters are going to be obtained from the second strain loaded

'''
def initialize_population(strain_1_name, strain_2_name, pop_concentration, strain_doubling_times):
	#I need to restart the ids
	
	concentrations=[int(elem) for elem in pop_concentration.split(':')]

	clusters=nx.Graph(graph_attr={'name': 'adjacency_list'})
	cells_to_divide=[]


	#process strain 1 clusters
	clusters, cells_to_divide, _, _=sample_from_cluster_population(clusters, cells_to_divide, concentrations[0], strain_1_name)

	#process strain 2 clusters
	clusters, cells_to_divide, cont_ids, cont_cluster_ids=sample_from_cluster_population(clusters, cells_to_divide, concentrations[1], strain_2_name)

	return([clusters, cells_to_divide, cont_ids, cont_cluster_ids])


'''
This function is going to create a network from one of the already grown strains. It is going to sample a cluster_id
and create that network and it is going to repeat the process until the amount of cells goes above the number of 
desired cells (this could be modified for the code to actually try to start at that amount if the sizes of the clusters
are calculated)
'''
def sample_from_cluster_population(clusters, cells_to_divide, max_concentration, strain_name):

	temp_dir=input_variables['temp_dir']

	#load edges of all networks
	network_edges=pd.read_csv(os.path.join(temp_dir, 'initial_population_'+strain_name+'.csv'))

	#load all nodes information
	nodes_information=pd.read_csv(os.path.join(temp_dir, 'initial_node_inf_'+strain_name+'.csv'))
	
	unique_cluster_ids=list(np.unique(nodes_information.cluster_id))

	cont_size_population=0

	sampled_ids=[]

	log_file=open('log.txt', 'a')
	log_file.write('Loading networks for '+strain_name+'\n')

	while cont_size_population<=max_concentration:

		#sample a cluster id at random
		temp_id=random.sample(unique_cluster_ids, 1)[0]

		#checking that the same cluster is not added twice
		if(not temp_id in sampled_ids):

			#add nodes to the network
			temp_nodes=nodes_information[nodes_information['cluster_id']==temp_id]

			for index, row in temp_nodes.iterrows():
				clusters.add_node(row['node_id'], strain=row['strain'], number_divisions=row['number_divisions'], 
					cluster_generation=row['cluster_generation'], time_added=row['time_added'], cluster_id=row['cluster_id'],
					transfer_added=row['transfer_added'])

				#sample doubling time
				temp_next_doub=sample_doub_t(strain_doubling_times, row['number_divisions'], 0, strain_name) #current time is 0
				#add value to the list in order
				cells_to_divide=add_cells_list(temp_next_doub, row['node_id'], cells_to_divide)

			## Using namedtuples
			# for index, row in temp_nodes.iterrows():
			# 	node_info = NodeInfo(strain=row['strain'], number_divisions=row['number_divisions'],
			# 						 cluster_generation=row['cluster_generation'], time_added=row['time_added'],
			# 						 cluster_id=row['cluster_id'], transfer_added=row['transfer_added'])
			# 	clusters.add_node(row['node_id'], info=node_info)

			# 	# sample doubling time
			# 	temp_next_doub = sample_doub_t(strain_doubling_times, node_info.number_divisions, 0, strain_name)
			# 	# add value to the list in order
			# 	cells_to_divide = add_cells_list(temp_next_doub, row['node_id'], cells_to_divide)

			# add edges to the network
			temp_edges=network_edges[network_edges['cluster_id']==temp_id]

			for index, row in temp_edges.iterrows():
				clusters.add_edge(row['node1'], row['node2'])

			sampled_ids.append(temp_id)

			cont_size_population+=len(temp_nodes)

			log_file.write('cluster_id: '+str(temp_id)+' Number of cells: '+str(len(temp_nodes))+'\n')

	log_file.write('Total amount of cells added: '+str(cont_size_population)+'/'+str(max_concentration)+'\n')
	log_file.write('Number of clusters added: '+str(len(sampled_ids))+'\n')
	log_file.close()

	cont_ids=max(list(nodes_information.node_id))+1
	cont_cluster_ids=max(list(nodes_information.cluster_id))+1

	return([clusters, cells_to_divide, cont_ids, cont_cluster_ids])





def load_surviving_population(curr_selec_round, temp_dir, strain_doubling_times):

	#load edges of all networks
	network_edges=pd.read_csv(os.path.join(temp_dir, 'final_population_'+str(curr_selec_round-1)+'.csv'))
	cont_cluster_id=max(network_edges.cluster_id)+1

	#load all nodes information
	nodes_information=pd.read_csv(os.path.join(temp_dir, 'node_inf_'+str(curr_selec_round-1)+'.csv'))
	cont_ids=max(nodes_information.node_id)+1

	#read surviving networks ids
	file=open(os.path.join(temp_dir, 'surviving_ids_'+str(curr_selec_round-1)+'.txt'), 'r')
	list_ids=[int(line.rstrip('\n').split(',')[0]) for line in file.readlines()]
	file.close()

	#iterating through the surviving ids and make the networks for those clusters that survived and sample their starting doubling times
	clusters=nx.Graph(graph_attr={'name': 'adjacency_list'})
	cells_to_divide=[]

	for temp_id in list_ids:
		#add nodes to the network
		temp_nodes=nodes_information[nodes_information['cluster_id']==temp_id]

		for index, row in temp_nodes.iterrows():
			clusters.add_node(row['node_id'], strain=row['strain'], number_divisions=row['number_divisions'], 
				cluster_generation=row['cluster_generation'], time_added=row['time_added'], cluster_id=row['cluster_id'],
				transfer_added=row['transfer_added'])

			#sample doubling time
			temp_next_doub=sample_doub_t(strain_doubling_times, row['number_divisions'], 0, strain_2_name) #current time is 0
			#add value to the list in order
			cells_to_divide=add_cells_list(temp_next_doub, row['node_id'], cells_to_divide)

		## Using namedtuples
		# for index, row in temp_nodes.iterrows():
		# 	node_info = NodeInfo(strain=row['strain'], number_divisions=row['number_divisions'],
		# 						 cluster_generation=row['cluster_generation'], time_added=row['time_added'],
		# 						 cluster_id=row['cluster_id'], transfer_added=row['transfer_added'])
		# 	clusters.add_node(row['node_id'], info=node_info)

		# 	# sample doubling time
		# 	temp_next_doub = sample_doub_t(strain_doubling_times, node_info.number_divisions, 0, strain_2_name)
		# 	# add value to the list in order
		# 	cells_to_divide = add_cells_list(temp_next_doub, row['node_id'], cells_to_divide)

		# add edges to the network
		temp_edges=network_edges[network_edges['cluster_id']==temp_id]

		for index, row in temp_edges.iterrows():
			clusters.add_edge(row['node1'], row['node2'])

	return([clusters, cells_to_divide, cont_ids, cont_cluster_id])


'''
In the function divide cells all the cells that are dividing in the current time point are going to divided
-cluster_population is the networkx Graph object with the networks of all the clusters in the population
-cells_to_divide is the list from where the current time and the list of cells that are going to be dividing this time
	is going to be extracted from
Note: the function doesn't return the cluster_population variable as it receives a pointer to the object and the real
variable is being modified
'''
def divide_cells(cluster_population, cells_to_divide, dict_doub_t_dist, cont_ids):

	curr_selec_round=input_variables['curr_selec_round']

	curr_time=cells_to_divide[0][0]
	list_ids=cells_to_divide[0][1] #ids of dividing cells

	cont_mother_id_bigger=0

	for mother_id in list_ids:
		mother_strain=cluster_population.nodes[mother_id]['strain']
		mother_clust_gen=cluster_population.nodes[mother_id]['cluster_generation']
		mother_clust_id=cluster_population.nodes[mother_id]['cluster_id']

		#create daughter cell
		daughter_id=cont_ids
		cluster_population.add_node(cont_ids, strain=mother_strain, number_divisions=0, cluster_generation=mother_clust_gen,
			time_added=curr_time, cluster_id=mother_clust_id, transfer_added=curr_selec_round)
		cluster_population.add_edge(mother_id, daughter_id)
		temp_daughter_time=sample_doub_t(dict_doub_t_dist, 0, curr_time, mother_strain) #using 0 instead of cluster_population.nodes[daughter_id]["number_divisions"]
		cells_to_divide=add_cells_list(temp_daughter_time, daughter_id, cells_to_divide)
		cont_ids+=1

		if(daughter_id<mother_id):
			cont_mother_id_bigger+=1

		#divide mother cell
		cluster_population.nodes[mother_id]["number_divisions"]+=1 #adding that the mother divided one more time (because of the step before)
		temp_mother_div=cluster_population.nodes[mother_id]["number_divisions"]
		temp_mother_time=sample_doub_t(dict_doub_t_dist, temp_mother_div, curr_time, mother_strain)
		cells_to_divide=add_cells_list(temp_mother_time, mother_id, cells_to_divide)

	## Using namedtuples
	# for mother_id in list_ids:
	# 	mother_info = cluster_population.nodes[mother_id]['info']

	# 	# create daughter cell
	# 	daughter_id = cont_ids
	# 	daughter_info = NodeInfo(strain=mother_info.strain, number_divisions=0, cluster_generation=mother_info.cluster_generation,
	# 							 time_added=curr_time, cluster_id=mother_info.cluster_id, transfer_added=curr_selec_round)
	# 	cluster_population.add_node(cont_ids, info=daughter_info)
	# 	cluster_population.add_edge(mother_id, daughter_id)
	# 	temp_daughter_time = sample_doub_t(dict_doub_t_dist, 0, curr_time, mother_info.strain)
	# 	cells_to_divide = add_cells_list(temp_daughter_time, daughter_id, cells_to_divide)
	# 	cont_ids += 1

	# 	# divide mother cell
	# 	updated_mother_info = mother_info._replace(number_divisions=mother_info.number_divisions + 1)
	# 	cluster_population.nodes[mother_id]['info'] = updated_mother_info
	# 	temp_mother_time = sample_doub_t(dict_doub_t_dist, updated_mother_info.number_divisions, curr_time, mother_info.strain)
	# 	cells_to_divide = add_cells_list(temp_mother_time, mother_id, cells_to_divide)

	if(cont_mother_id_bigger>0):
		print('Mother id is bigger, '+str(cont_mother_id_bigger))

	return([cont_ids, cells_to_divide])


# Function to calculate the highest edge degree of a network, it returns an edge at random if there are ties
def calculate_max_edge_degree(network):
	# Initialize result variables
	max_edge_degree = -1
	max_edges = []  # List to store edges with the highest edge degree

	for edge in network.edges():
		node1, node2 = edge
		degree_node1 = network.degree(node1)
		degree_node2 = network.degree(node2)
		edge_degree = degree_node1 + degree_node2 - 2

		# Check if the current edge has a higher degree than the current max
		if edge_degree > max_edge_degree:
			max_edge_degree = edge_degree
			max_edges = [edge]  # Reset the list with the new highest degree edge
		elif edge_degree == max_edge_degree:
			max_edges.append(edge)  # Add edge to list if it ties with the current max

	# Return one of the edges at random if there are ties, or the single edge if it's the only one
	if len(max_edges) > 1:
		return [random.choice(max_edges), max_edge_degree]
	else:
		# If there's only one edge with the highest value, return that edge
		return [max_edges[0], max_edge_degree] if max_edges else [None, max_edge_degree]


def fragmentation_edge_degree_population(cluster_population, edge_deg_threshold, cont_cluster_ids):
	# create a copy of all networks 
	list_clusters=[cluster_population.subgraph(c).copy() for c in nx.connected_components(cluster_population)]

	for temp_cluster in list_clusters:
		#calculate the highest edge degree and the edge with the highest edge degree
		edge_to_remove, edge_degree=calculate_max_edge_degree(temp_cluster)

		if(edge_degree>=edge_deg_threshold):

			cluster_population.remove_edge(*edge_to_remove)

			#update cluster generations in all nodes of the cluster that is fracturing (both in the propagule and in the parent)
			for temp_node in temp_cluster.nodes():
				cluster_population.nodes[temp_node]['cluster_generation']+=1

			## Using namedtuple
			# for temp_node in temp_cluster.nodes():
			# 	node_info = cluster_population.nodes[temp_node]['info']
			# 	updated_info = node_info._replace(cluster_generation=node_info.cluster_generation + 1)
			# 	cluster_population.nodes[temp_node]['info'] = updated_info

			#update cluster ids for the propagule
			#use temp cluster to find the nodes of the propagule after fracture
			fractured_temp_cluster=temp_cluster.copy()
			fractured_temp_cluster.remove_edge(*edge_to_remove)
			smallest_component = min(nx.connected_components(fractured_temp_cluster), key=len)
			for temp_node in list(smallest_component):
				cluster_population.nodes[temp_node]['cluster_id']=cont_cluster_ids

			## Using namedtuples
			# for temp_node in list(smallest_component):
			# 	node_info = cluster_population.nodes[temp_node]['info']
			# 	updated_info = node_info._replace(cluster_id=cont_cluster_ids)
			# 	cluster_population.nodes[temp_node]['info'] = updated_info
			
			cont_cluster_ids+=1

	return(cont_cluster_ids)


def write_output_information_timepoint(curr_time, cluster_population, dict_files):
	#time registry file
	#sim_number, strain, num_clusters, time, transfer, num_cells, total_clusters, total_cells
	strain_1_name=input_variables['strain_1_name']
	strain_2_name=input_variables['strain_2_name']

	curr_selec_round=input_variables['curr_selec_round']
	sim_number=input_variables['sim_number']

	cont_cells_strain1=0
	cont_cells_strain2=0

	cont_clusters_strain1=0
	cont_clusters_strain2=0


	for cluster in [list(i) for i in nx.connected_components(cluster_population)]:
		cluster_size=len(cluster)
		temp_clust_strain=cluster_population.nodes[cluster[0]]['strain']

	## Using namedtuple
	# for cluster in [list(i) for i in nx.connected_components(cluster_population)]:
	# 	cluster_size = len(cluster)
	# 	temp_clust_strain = cluster_population.nodes[cluster[0]]['info'].strain

		if(temp_clust_strain==strain_1_name):
			cont_cells_strain1+=cluster_size
			cont_clusters_strain1+=1
		else:
			cont_cells_strain2+=cluster_size
			cont_clusters_strain2+=1

	total_cells=cont_cells_strain1+cont_cells_strain2
	total_clusters=cont_clusters_strain1+cont_clusters_strain2

	dict_files['time_registry'].write(str(sim_number)+','+strain_1_name+','+str(cont_clusters_strain1)+','+
		str(curr_time)+','+str(curr_selec_round)+','+str(cont_cells_strain1)+','+str(total_clusters)+','+str(total_cells)+'\n')

	dict_files['time_registry'].write(str(sim_number)+','+strain_2_name+','+str(cont_clusters_strain2)+','+
		str(curr_time)+','+str(curr_selec_round)+','+str(cont_cells_strain2)+','+str(total_clusters)+','+str(total_cells)+'\n')


def write_proportion_populations(initial_population, final_population, dict_files, final_time):
	#write proportions file
	# sim_number, transfer, experiment_phase, time, strain1, strain2, prop_cluster_pop1, prop_cluster_pop2, prop_cell_pop1, prop_cell_pop2,
	#total_clusters, total_cells
	curr_selec_round=input_variables['curr_selec_round']
	sim_number=input_variables['sim_number']

	strain_1_name=input_variables['strain_1_name']
	strain_2_name=input_variables['strain_2_name']


	#calculate values for the initial population
	ini_cont_cells_strain1=0
	ini_cont_cells_strain2=0

	ini_cont_clusters_strain1=0
	ini_cont_clusters_strain2=0

	for cluster in [list(i) for i in nx.connected_components(initial_population)]:
		cluster_size=len(cluster)

		temp_clust_strain=initial_population.nodes[cluster[0]]['strain']

		if(temp_clust_strain==strain_1_name):
			ini_cont_cells_strain1+=cluster_size
			ini_cont_clusters_strain1+=1
		else:
			ini_cont_cells_strain2+=cluster_size
			ini_cont_clusters_strain2+=1

	ini_total_cells=ini_cont_cells_strain1+ini_cont_cells_strain2
	ini_total_clusters=ini_cont_clusters_strain1+ini_cont_clusters_strain2

	#calculate values for the final population
	fin_cont_cells_strain1=0
	fin_cont_cells_strain2=0

	fin_cont_clusters_strain1=0
	fin_cont_clusters_strain2=0

	for cluster in [list(i) for i in nx.connected_components(final_population)]:
		cluster_size=len(cluster)

		temp_clust_strain=final_population.nodes[cluster[0]]['strain']

		if(temp_clust_strain==strain_1_name):
			fin_cont_cells_strain1+=cluster_size
			fin_cont_clusters_strain1+=1
		else:
			fin_cont_cells_strain2+=cluster_size
			fin_cont_clusters_strain2+=1

	fin_total_cells=fin_cont_cells_strain1+fin_cont_cells_strain2
	fin_total_clusters=fin_cont_clusters_strain1+fin_cont_clusters_strain2

	dict_files['selection_proportion'].write(str(sim_number)+','+str(curr_selec_round)+','+str(final_time)+','+strain_1_name+','+
		strain_2_name+','+str(ini_cont_clusters_strain1)+','+str(ini_cont_clusters_strain2)+','+
		str(ini_cont_cells_strain1)+','+str(ini_cont_cells_strain2)+','+str(ini_total_clusters)+','+str(ini_total_cells)+','+
		str(fin_cont_clusters_strain1)+','+str(fin_cont_clusters_strain2)+','+
		str(fin_cont_cells_strain1)+','+str(fin_cont_cells_strain2)+','+str(fin_total_clusters)+','+str(fin_total_cells)+'\n')


## Using namedtuples
# def write_proportion_populations(initial_population, final_population, dict_files, final_time):
# 	# write proportions file
# 	# sim_number, transfer, experiment_phase, time, strain1, strain2, prop_cluster_pop1, prop_cluster_pop2, prop_cell_pop1, prop_cell_pop2,
# 	# total_clusters, total_cells
# 	curr_selec_round = input_variables['curr_selec_round']
# 	sim_number = input_variables['sim_number']

# 	strain_1_name = input_variables['strain_1_name']
# 	strain_2_name = input_variables['strain_2_name']

# 	# calculate values for the initial population
# 	ini_cont_cells_strain1 = 0
# 	ini_cont_cells_strain2 = 0

# 	ini_cont_clusters_strain1 = 0
# 	ini_cont_clusters_strain2 = 0

# 	for cluster in [list(i) for i in nx.connected_components(initial_population)]:
# 		cluster_size = len(cluster)
# 		temp_clust_strain = initial_population.nodes[cluster[0]]['info'].strain

# 		if temp_clust_strain == strain_1_name:
# 			ini_cont_cells_strain1 += cluster_size
# 			ini_cont_clusters_strain1 += 1
# 		else:
# 			ini_cont_cells_strain2 += cluster_size
# 			ini_cont_clusters_strain2 += 1

# 	ini_total_cells = ini_cont_cells_strain1 + ini_cont_cells_strain2
# 	ini_total_clusters = ini_cont_clusters_strain1 + ini_cont_clusters_strain2

# 	# calculate values for the final population
# 	fin_cont_cells_strain1 = 0
# 	fin_cont_cells_strain2 = 0

# 	fin_cont_clusters_strain1 = 0
# 	fin_cont_clusters_strain2 = 0

# 	for cluster in [list(i) for i in nx.connected_components(final_population)]:
# 		cluster_size = len(cluster)
# 		temp_clust_strain = final_population.nodes[cluster[0]]['info'].strain

# 		if temp_clust_strain == strain_1_name:
# 			fin_cont_cells_strain1 += cluster_size
# 			fin_cont_clusters_strain1 += 1
# 		else:
# 			fin_cont_cells_strain2 += cluster_size
# 			fin_cont_clusters_strain2 += 1

# 	fin_total_cells = fin_cont_cells_strain1 + fin_cont_cells_strain2
# 	fin_total_clusters = fin_cont_clusters_strain1 + fin_cont_clusters_strain2

# 	dict_files['selection_proportion'].write(str(sim_number) + ',' + str(curr_selec_round) + ',' + str(final_time) + ',' + strain_1_name + ',' +
# 		strain_2_name + ',' + str(ini_cont_clusters_strain1) + ',' + str(ini_cont_clusters_strain2) + ',' +
# 		str(ini_cont_cells_strain1) + ',' + str(ini_cont_cells_strain2) + ',' + str(ini_total_clusters) + ',' + str(ini_total_cells) + ',' +
# 		str(fin_cont_clusters_strain1) + ',' + str(fin_cont_clusters_strain2) + ',' +
# 		str(fin_cont_cells_strain1) + ',' + str(fin_cont_cells_strain2) + ',' + str(fin_total_clusters) + ',' + str(fin_total_cells) + '\n')


def write_ending_population(cluster_population):
	#final population file
	# strain, node1, node2, transfer_added, time_grown, cluster_id, cluster_generation,

	# create a copy of all networks 
	list_clusters=[cluster_population.subgraph(c).copy() for c in nx.connected_components(cluster_population)]

	cont_mother_id_bigger=0

	for temp_cluster in list_clusters:
		temp_strain=temp_cluster.nodes[list(temp_cluster.nodes())[0]]['strain']
		temp_cluster_id=temp_cluster.nodes[list(temp_cluster.nodes())[0]]['cluster_id']
		temp_clust_gen=temp_cluster.nodes[list(temp_cluster.nodes())[0]]['cluster_generation']

		## Using namedtuples (make sure to unindent if used)
		# for temp_cluster in list_clusters:
		# 	temp_strain = temp_cluster.nodes[list(temp_cluster.nodes())[0]]['info'].strain
		# 	temp_cluster_id = temp_cluster.nodes[list(temp_cluster.nodes())[0]]['info'].cluster_id
		# 	temp_clust_gen = temp_cluster.nodes[list(temp_cluster.nodes())[0]]['info'].cluster_generation

		for edge in temp_cluster.edges():
			node1, node2 = edge

			#sometimes the ids are in the incorrect order (daughter in node1 or in ndoe2), so the mother id is going to be
			#the smallest id to avoid errors in the programs later on
			if(node1<node2):
				mother_id=node1
				daughter_id=node2
			else:
				mother_id=node2
				daughter_id=node1

			#Time the daughter cell was created
			temp_time_added=temp_cluster.nodes[daughter_id]['time_added']
			temp_transfer_added=temp_cluster.nodes[daughter_id]['transfer_added']

			## Namedtuple
			# temp_time_added = temp_cluster.nodes[daughter_id]['info'].time_added
			# temp_transfer_added = temp_cluster.nodes[daughter_id]['info'].transfer_added

			dict_files['final_population'].write(temp_strain+','+str(mother_id)+','+str(daughter_id)+','+str(temp_transfer_added)+','+
				str(temp_time_added)+','+str(temp_cluster_id)+','+str(temp_clust_gen)+'\n')



def write_nodes_information(clusters):

	curr_selec_round=input_variables['curr_selec_round']
	temp_dir=input_variables['temp_dir']

	# I need to create the file to save this information
	#node_id, strain, number_divisions, cluster_generation, time_added, cluster_id

	node_data = []

	# Iterate through each node in the network
	for node_id in clusters.nodes():
		node_info = clusters.nodes[node_id]

		## Using namedtuples (unindent if used)
		# for node_id in clusters.nodes():
		# 	node_info = clusters.nodes[node_id]['info']
		
		# Extract node information
		strain = node_info.get('strain', None)
		number_divisions = node_info.get('number_divisions', None)
		cluster_generation = node_info.get('cluster_generation', None)
		time_added = node_info.get('time_added', None)
		cluster_id = node_info.get('cluster_id', None)
		temp_transfer_added=node_info.get('transfer_added', None)
		
		## Namedtuple
		# strain = node_info.strain
		# number_divisions = node_info.number_divisions
		# cluster_generation = node_info.cluster_generation
		# time_added = node_info.time_added
		# cluster_id = node_info.cluster_id
		# temp_transfer_added = node_info.transfer_added

		# Append the node information as a dictionary to the list
		node_data.append({
			'node_id': node_id,
			'strain': strain,
			'number_divisions': number_divisions,
			'cluster_generation': cluster_generation,
			'time_added': time_added,
			'cluster_id': cluster_id,
			'transfer_added': temp_transfer_added 
		})

		## Namedtuple
		# node_data.append({
		# 	'node_id': node_id,
		# 	'strain': node_info.strain,
		# 	'number_divisions': node_info.number_divisions,
		# 	'cluster_generation': node_info.cluster_generation,
		# 	'time_added': node_info.time_added,
		# 	'cluster_id': node_info.cluster_id,
		# 	'transfer_added': node_info.transfer_added
		# })

	# Convert the list of dictionaries to a pandas DataFrame
	df = pd.DataFrame(node_data)

	# Save the DataFrame to a CSV file
	df.to_csv(os.path.join(temp_dir, 'node_inf_'+str(curr_selec_round)+'.csv'), index=False)



'''
-final_population_#: csv file with the edges, the time they were added and the id of the cluster
To have the time of when a edge was added I need to keep track of the daughters of each cell
and save the time when they were added (so a list of 2 entries [daughter_cell_id, time])
Cluster id needs to be updated every time a cluster fractures
	strain, node1, node2, transfer_added, time_grown, cluster_id, cluster_generation,

-time_registry: csv file with population at each time point
	sim_number, strain, num_clusters, time, transfer, num_cells, total_clusters, total_cells

-selection_proportions: csv file of proportions of clusters after the growth and settling round
	sim_number, transfer, experiment_phase, end, strain1, strain2, prop_cluster_pop1, prop_cluster_pop2, prop_cell_pop1, prop_cell_pop2,
	total_clusters, total_cells 
	experiment_phase [growth/settling]
	end [True/False] Is it the end of the phase? or the beginning?
'''

#Function to initialize all output file, a dictionary with all the opened files is returned for it
#to be accesed by all functions
def initialize_files(input_variables):
	output_dir=input_variables['output_dir']
	curr_selec_round=input_variables['curr_selec_round']

	dict_files={}

	#Initializing results file if they don't exist
	final_pop_file=os.path.join(temp_dir,"final_population_"+str(curr_selec_round)+".csv")
	if(os.path.exists(final_pop_file)):
		dict_files['final_population']=open(final_pop_file, "a")
	else:
		dict_files['final_population']=open(final_pop_file, "w")
		dict_files['final_population'].write('strain,node1,node2,transfer_added,time_grown,cluster_id,cluster_generation\n')

	time_registry_file=os.path.join(output_dir,"time_registry.csv")
	if(os.path.exists(time_registry_file)):
		dict_files['time_registry']=open(time_registry_file, "a")
	else:
		dict_files['time_registry']=open(time_registry_file, "w")
		dict_files['time_registry'].write('sim_number,strain,num_clusters,time,transfer,number_cells,total_clusters,total_cells\n')

	selection_proportions_file=os.path.join(output_dir,"proportions_growth_registry.csv")
	if(os.path.exists(selection_proportions_file)):
		dict_files['selection_proportion']=open(selection_proportions_file, "a")
	else:
		dict_files['selection_proportion']=open(selection_proportions_file, "w")
		dict_files['selection_proportion'].write('sim_number,transfer,time,strain1,strain2,clusters_pop1_b,clusters_pop2_b,'+
			'cells_pop1_b,cells_pop2_b,total_clusters_b,total_cells_b,clusters_pop1_a,clusters_pop2_a,cells_pop1_a,'+
			'cells_pop2_a,total_clusters_a,total_cells_a\n')

	return(dict_files)


#function to close all files opened in the dictionary
def close_all_files(dict_files):
	for file_obj in dict_files.values():
		file_obj.close()


#test function to know how much memory is being used by the script
def memory_usage():
	process = psutil.Process(os.getpid())
	mem_info = process.memory_info()
	return(mem_info.rss)


def get_graph_size(G):
	graph_size = sys.getsizeof(G)
	for node, data in G.nodes(data=True):
		graph_size += sys.getsizeof(node) #if using the dictionary this should give the ram usage 
		#of the dictionary

		#only when using the named tuples
		# info = data['info']
		# graph_size += sys.getsizeof(info)
	return(graph_size)


def get_size(obj, seen=None):
    """Recursively finds size of objects"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_size(v, seen) for v in obj.values()])
        size += sum([get_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_size(i, seen) for i in obj])
    return size


#### MAIN ####

'''
to do 29mar
Add carrying capacity as an input
modify the function to load the clusters for the initial growth

'''

#Start memory usage
start_mem = memory_usage()

# initializing variables 

strain_1_name=args.strain_names.split(':')[0]
strain_2_name=args.strain_names.split(':')[1]

strain_doubling_times={}
strain_doubling_times[strain_1_name]=load_doubling_time(args.strain1_dt, "minutes")
strain_doubling_times[strain_2_name]=load_doubling_time(args.strain2_dt, "minutes")

temp_dir=args.temp_dir
output_dir=args.output_dir

edge_deg_threshold=int(args.edge_degree)

curr_selec_round=int(args.curr_selec_round)

sim_number=int(args.sim_number)

carrying_capacity=int(args.carrying_cap) #maximum number of cells in the population (before was 10**6)

input_variables={'strain_1_name':strain_1_name,
					'strain_2_name':strain_2_name,
					'strain_doubling_times':strain_doubling_times,
					'temp_dir':temp_dir,
					'output_dir':output_dir,
					'edge_deg_threshold':edge_deg_threshold,
					'curr_selec_round':curr_selec_round,
					'sim_number':sim_number}

# print(input_variables)

# Create or open files
dict_files=initialize_files(input_variables)


#hard-coded variables (I could later move them to be given as an input)
max_minutes=24*60 #24 hours * 60 minutes
# time_step=5 #5 minutes from time lapses
#If time is continues 

# NodeInfo = namedtuple('NodeInfo', ['strain', 'number_divisions', 'cluster_generation', 'time_added', 'cluster_id', 'transfer_added'])


start_time = time.time()

# Initializing or reloading populations
if (curr_selec_round==1):
	#Initialize population
	cluster_population, cells_to_divide, cont_ids, cont_cluster_ids=initialize_population(strain_1_name, strain_2_name, args.pop_concentration, strain_doubling_times)
else:
	#Load population
	cluster_population, cells_to_divide, cont_ids, cont_cluster_ids=load_surviving_population(curr_selec_round, temp_dir, strain_doubling_times)

#saving initial population to track the initial and the final proportion
initial_population=cluster_population.copy()

# Main growth loop
curr_pop_size=len(cluster_population.nodes())
# len([c for c in nx.connected_components(cluster_population)]) #this can be used to count the total amount of clusters


curr_time=cells_to_divide[0][0]

# print(curr_time, curr_pop_size<carrying_capacity, curr_time<max_minutes)

while (curr_pop_size<carrying_capacity) & (curr_time<=max_minutes):

	# print(curr_time, curr_pop_size<carrying_capacity, curr_time<max_minutes)
	# print(curr_time, curr_pop_size, end='\r')

	#divide cells at this time point
	cont_ids, cells_to_divide=divide_cells(cluster_population, cells_to_divide, strain_doubling_times, cont_ids) #it gets a reference to the network variable so it doesn't need to return the variable

	#check highest edge degree for all clusters to trigger fragmentation if any goes beyond the threshold
	cont_cluster_ids=fragmentation_edge_degree_population(cluster_population, edge_deg_threshold, cont_cluster_ids)

	#save information of this timepoint
	write_output_information_timepoint(curr_time, cluster_population, dict_files)

	#delete the times in cells to divide of the current time point
	# cells_to_divide=update_cells_to_divide(cells_to_divide)
	cells_to_divide.pop(0)

	#update current time and population size
	# curr_time=get_smallest_time(cells_to_divide)
	last_time=curr_time
	curr_time=cells_to_divide[0][0] #select the next smallest time
	curr_pop_size=len(cluster_population.nodes())

log_file=open('log.txt', 'a')
if(curr_time>max_minutes):
	# print('Max time reached')
	log_file.write('Growth phase: Max time reached\n')
else:
	# print('Carrying capacity reached')
	log_file.write('Growth phase: Carrying capacity reached\n')
	# if carrying capacity was reached then the information of the timepoint is going to be written again
	write_output_information_timepoint(max_minutes, cluster_population, dict_files)


#Final memory usage
end_mem = memory_usage()

log_file.write("Final Ram usage: "+str((end_mem-start_mem)/1024/1024)+" MB\n")
log_file.write("cluster_population Ram usage: "+str(get_graph_size(cluster_population)/1024/1024)+" MB\n")
log_file.write("cells_to_divide Ram usage: "+str(get_size(cells_to_divide)/1024/1024)+" MB\n")

end_time = time.time()
execution_time = end_time - start_time

log_file.write("Execution time: "+str(execution_time)+" seconds\n\n")

log_file
log_file.close()

if(curr_time>max_minutes):
	final_time=max_minutes
else:
	final_time=last_time

# print([len(c) for c in sorted(nx.connected_components(cluster_population), key=len, reverse=True)])
# print([cluster_population.nodes[list(c)[0]]['cluster_id'] for c in sorted(nx.connected_components(cluster_population), key=len, reverse=True)])

write_proportion_populations(initial_population, cluster_population, dict_files, final_time)

write_ending_population(cluster_population)

write_nodes_information(cluster_population)

#closing all files from the dictionary
close_all_files(dict_files)
