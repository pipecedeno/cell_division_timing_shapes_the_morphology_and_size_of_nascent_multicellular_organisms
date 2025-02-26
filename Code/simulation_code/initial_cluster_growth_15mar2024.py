#!/usr/bin/env python3

'''
Date: 29mar2024
This code is part of the population simulations used to calculate the selection rate of the growth phase and the settling
phase as if it was a competition experiment in the lab. This code is used to initialize a population of clusters before the 
initial growth phase of the experiment.

inputs:
-doubling time strain
-initial concentration 
-strain name (it is necessary to know how to name the output file)
-temp_dir: directory where the files of the edges are going to be saved and the ids of the clusters that survive settling selection
-output_dir: Here the files time_registry and selection_proportions
-edge degree threshold
-cluster generations
-mean_doubling_times: for how many mean doubling times are they allowed to grow after the number of generations is reached.
-other_strain_clusts: this is the other strain node id file, it is just to get the cluster_id and node_id to not have repeating ids
						if the other strain hasn't been processed then it needs to receive 'none' to start the ids as 1


	

'''

import os
import pandas as pd
import argparse
import networkx as nx
# import matplotlib.pyplot as plt
# from random import choice
import random
import numpy as np

# import netlsd
# from scipy.stats import skew, kurtosis
import time

#only used for testing purposes
# from itertools import permutations

parser = argparse.ArgumentParser()
parser.add_argument('-f','--strain_dt',dest="strain_dt",required=True) #csv file of the doubling time distributions
parser.add_argument('-s','--other_strain_clusts',dest="other_strain_clusts",required=True) #other strain node_id file or none
parser.add_argument('-n','--strain_name',dest="strain_name",required=True) #name of the strains (format: strain1:strain2)
parser.add_argument('-t','--temp_dir',dest="temp_dir",required=True) #temp directory
parser.add_argument('-e','--edge_degree',dest="edge_degree",required=True) #edge degree threshold
parser.add_argument('-c','--pop_concentration',dest="pop_concentration",required=True) #starting concentration of the population (integer number)
parser.add_argument('-k','--max_clust_gen',dest="max_clust_gen",required=True) #cluster generation needed to be reached to end the program
args = parser.parse_args()


#### FUNCTIONS #####

#This function returns a dictionary with where each key is the number of division and the value associated to that key
#is the observed doubling times
#function modified to also return the mean doubling time
def load_doubling_time(path_to_file, column):
	temp_file=pd.read_csv(path_to_file, header=0)

	# Group the data by 'division_number' and extract the 'minutes' values into a list
	grouped_data = temp_file.groupby('division_number')['minutes'].apply(list)

	# Create a dictionary where the key is the unique value of 'division_number' and the value is the list of 'minutes'
	doub_t_dist = dict(zip(grouped_data.index, grouped_data.values))

	mean_doubling_time=np.ceil(np.mean(temp_file.minutes))

	return([doub_t_dist, mean_doubling_time])


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


#this is a modified function from the growth phase of the simulation, instead of initializing 2 populations
#it is initializing only one population and instead of receiving 2 strains to start it is only receiving
#on strain, therefore pop_concentration is just an integer instead of a string
def initialize_population(strain_name, pop_concentration, strain_doubling_times, cont_ids, cont_cluster_ids):

	temp_concentration=int(pop_concentration)

	clusters=nx.Graph()
	cells_to_divide=[]


	# create the cells for the first strain
	for i in range(temp_concentration):
		#add node
		clusters.add_node(cont_ids, strain=strain_name, number_divisions=0, cluster_generation=1,
			time_added=0, cluster_id=cont_cluster_ids, transfer_added=0)
		#Just when initializing the clusters the cluster_id is going to be the same as the cell ids
		
		#sample doubling time
		temp_next_doub=sample_doub_t(strain_doubling_times, 0, 0, strain_name)

		#add value to the list in order
		cells_to_divide=add_cells_list(temp_next_doub, cont_ids, cells_to_divide)

		cont_ids+=1
		cont_cluster_ids+=1

	return([clusters, cells_to_divide, cont_ids, cont_cluster_ids])



'''
In the function divide cells all the cells that are dividing in the current time point are going to divided
-cluster_population is the networkx Graph object with the networks of all the clusters in the population
-cells_to_divide is the list from where the current time and the list of cells that are going to be dividing this time
	is going to be extracted from
Note: the function doesn't return the cluster_population variable as it receives a pointer to the object and the real
variable is being modified
It was modified to make the selection_round=0 so that when the cluster is being built in matlab no error happens
'''
def divide_cells(cluster_population, cells_to_divide, dict_doub_t_dist, cont_ids):

	curr_selec_round=0

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


def fragmentation_edge_degree_population(cluster_population, edge_deg_threshold, cont_cluster_ids, max_cluster_generation):
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

			if(cluster_population.nodes[temp_node]['cluster_generation']>max_cluster_generation):
				max_cluster_generation=cluster_population.nodes[temp_node]['cluster_generation']

			#update cluster ids for the propagule
			#use temp cluster to find the nodes of the propagule after fracture
			fractured_temp_cluster=temp_cluster.copy()
			fractured_temp_cluster.remove_edge(*edge_to_remove)
			smallest_component = min(nx.connected_components(fractured_temp_cluster), key=len)
			for temp_node in list(smallest_component):
				cluster_population.nodes[temp_node]['cluster_id']=cont_cluster_ids
			
			cont_cluster_ids+=1

	return([cont_cluster_ids, max_cluster_generation])



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

			dict_files['final_population'].write(temp_strain+','+str(mother_id)+','+str(daughter_id)+','+str(temp_transfer_added)+','+
				str(temp_time_added)+','+str(temp_cluster_id)+','+str(temp_clust_gen)+'\n')



def write_nodes_information(clusters):

	# curr_selec_round=input_variables['curr_selec_round']
	temp_dir=input_variables['temp_dir']
	strain_name=input_variables['strain_name']

	# I need to create the file to save this information
	#node_id, strain, number_divisions, cluster_generation, time_added, cluster_id

	node_data = []

	# Iterate through each node in the network
	for node_id in clusters.nodes():
	    node_info = clusters.nodes[node_id]
	    
	    # Extract node information
	    strain = node_info.get('strain', None)
	    number_divisions = node_info.get('number_divisions', None)
	    cluster_generation = node_info.get('cluster_generation', None)
	    time_added = node_info.get('time_added', None)
	    cluster_id = node_info.get('cluster_id', None)
	    temp_transfer_added=node_info.get('transfer_added', None)
	    
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

	# Convert the list of dictionaries to a pandas DataFrame
	df = pd.DataFrame(node_data)

	# Save the DataFrame to a CSV file
	df.to_csv(os.path.join(temp_dir, 'initial_node_inf_'+strain_name+'.csv'), index=False)



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
	# output_dir=input_variables['output_dir']
	# curr_selec_round=input_variables['curr_selec_round']
	strain_name=input_variables['strain_name']

	dict_files={}

	#Initializing results file if they don't exist
	final_pop_file=os.path.join(temp_dir,"initial_population_"+strain_name+".csv")
	if(os.path.exists(final_pop_file)):
		dict_files['final_population']=open(final_pop_file, "a")
	else:
		dict_files['final_population']=open(final_pop_file, "w")
		dict_files['final_population'].write('strain,node1,node2,transfer_added,time_grown,cluster_id,cluster_generation\n')

	return(dict_files)


#function to close all files opened in the dictionary
def close_all_files(dict_files):
	for file_obj in dict_files.values():
		file_obj.close()




#### MAIN ####


# initializing variables 

strain_name=args.strain_name

strain_doubling_times={}
strain_doubling_times[strain_name], mean_doubling_time=load_doubling_time(args.strain_dt, "minutes")

temp_dir=args.temp_dir

edge_deg_threshold=int(args.edge_degree)

if(args.other_strain_clusts=='none'):
	cont_ids=1
	cont_cluster_ids=1
else:
	nodes_information_other_strain=pd.read_csv(args.other_strain_clusts)
	cont_ids=max(list(nodes_information_other_strain.node_id))+1
	cont_cluster_ids=max(list(nodes_information_other_strain.cluster_id))+1

# curr_selec_round=int(args.curr_selec_round)

# sim_number=int(args.sim_number)

input_variables={'strain_name':strain_name,
					'strain_doubling_times':strain_doubling_times,
					'temp_dir':temp_dir,
					'edge_deg_threshold':edge_deg_threshold}

# print(input_variables)

# Create or open files
dict_files=initialize_files(input_variables)


#hard-coded variables (I could later move them to be given as an input)
cluster_generations_to_grow=int(args.max_clust_gen)
carrying_capacity=10**6 #maximum number of cells in the population
max_minutes=24*60 #24 hours * 60 minutes
# time_step=5 #5 minutes from time lapses


#Initialize population
cluster_population, cells_to_divide, cont_ids, cont_cluster_ids=initialize_population(strain_name, args.pop_concentration, strain_doubling_times, cont_ids, cont_cluster_ids)
max_cluster_generation=1 #This is the value used in initialize_population function

#saving initial population to track the initial and the final proportion
# initial_population=cluster_population.copy()

# Main growth loop
curr_pop_size=len(cluster_population.nodes())
# len([c for c in nx.connected_components(cluster_population)]) #this can be used to count the total amount of clusters


curr_time=cells_to_divide[0][0]

# print(curr_time, curr_pop_size<carrying_capacity, curr_time<max_minutes)

# print(mean_doubling_time*8)

# while (curr_time < (mean_doubling_time*8)):
while (max_cluster_generation <= cluster_generations_to_grow):

	# print(curr_time, curr_pop_size<carrying_capacity, curr_time<max_minutes)
	# print(curr_time, curr_pop_size, end='\r')

	#divide cells at this time point
	cont_ids, cells_to_divide=divide_cells(cluster_population, cells_to_divide, strain_doubling_times, cont_ids) #it gets a reference to the network variable so it doesn't need to return the variable

	#check highest edge degree for all clusters to trigger fragmentation if any goes beyond the threshold
	cont_cluster_ids, max_cluster_generation=fragmentation_edge_degree_population(cluster_population, edge_deg_threshold, cont_cluster_ids, max_cluster_generation)

	#delete the times in cells to divide of the current time point
	# cells_to_divide=update_cells_to_divide(cells_to_divide)
	cells_to_divide.pop(0)

	#update current time and population size
	# curr_time=get_smallest_time(cells_to_divide)
	last_time=curr_time
	curr_time=cells_to_divide[0][0] #select the next smallest time
	# print(curr_time)
	curr_pop_size=len(cluster_population.nodes())

log_file=open('log.txt', 'a')
log_file.write(strain_name+
	'\nFinal time: '+str(last_time/60)+' hours'+
	'\nNumber of clusters grown: '+str(len([c for c in nx.connected_components(cluster_population)]))+
	'\nMean cluster size: '+str(np.mean([len(c) for c in nx.connected_components(cluster_population)]))+
	'\nMean cluster generation: '+str(np.mean([cluster_population.nodes[list(c)[0]]['cluster_generation'] for c in nx.connected_components(cluster_population)]))+'\n')
log_file.close()

# print([len(c) for c in sorted(nx.connected_components(cluster_population), key=len, reverse=True)])
# print([cluster_population.nodes[list(c)[0]]['cluster_id'] for c in sorted(nx.connected_components(cluster_population), key=len, reverse=True)])

write_ending_population(cluster_population)

write_nodes_information(cluster_population)

#closing all files from the dictionary
close_all_files(dict_files)
