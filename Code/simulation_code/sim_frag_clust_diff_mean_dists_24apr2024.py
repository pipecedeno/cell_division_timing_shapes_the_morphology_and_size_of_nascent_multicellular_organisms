#!/usr/bin/env python3

'''
Date: 20oct2023

For this program fragmentation is going to happen when the highest edge degree
of the cluster gets higher than the threshold. The code will grow only one cluster at a time,
so after fragmentation the parent cluster (the one with the highest amount of cells) is discarded
to chech how similarity changes in the case where more differences are going to occur.

Note: this code cannot be executed in parallel, because netlsd already uses all cores

This code is a modification of sim_frag_pop_edge_degree_5oct2023.py.

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
	-fragmentation information (edges removed, size at fracture, proportion of the biggest fragment)
	-networks diameter
	-updates to fragmentation, it saves how many cells divided (updates) between fragmentations
	-Network files are only saved for the first network to reach multiples of 10

Modification history:
Date: 16Aug2023
The function save_degrees_to_csv was replaced by save_degree_distribution_to_csv. The idea of this
change is to save directly the degree distribution table into the "degree_dist" csv file instead
of the degree of each node, the idea of this change is to save storage of the information
saved and to avoid calculating the degree distribution in R after loading all the tables for each
simulation. This new table includes a column of frequency and another column of probability
to be able to plot the degree distribution more easily

Date: 7Sep2023
-Added the output file num_snowflakes_per_update.csv which saves how many clusters there
are by update after each fracture, this file is to know when the population reaches 50 clusters
and compare in which population it is being reached faster
-made the input_variabled dictionary a global variable, it is not longer being passed
as an argument to all functions
-Now the code saves the first network that reaches each generation divisible by %10, so 0, 10, 20...
The idea is to be able to visualize the networks after several generations

Date: 18Oct2023
-Changed the code so that the output files are in a dictionary and their are kept open throughout
the program execution
-Now the file size_at_fracture.csv won't be saved as that information is now contained in a column
on the fragmentation_inf.csv file


Date: 24Apr2024
-Changed the fragmentation part of the code to work similar to sim_frag_clust_edge_degree_after_all_cells_22feb2024.py
so that the fragmentation of the cluster happens at random and it always saves the information of the propagule
-Changed the while loop in the simulate one cluster function so that the function save_final_networks function is no longer
used and as now the while loop is going to stop by checking the size of cells to divide as in the break cluster function
it reinitialize the list if the desired size is reached
-deleted the network similarity functions and now it is not calculated at all 


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



#This function adds the new doubling time to the list of doubling times, because it is expected
#that the new doubling time is bigger than the elements of the list, it is going to iterate it
#in reverse to make this more efficient. The element input is formed of [time until doubling, id 
#of cell], and ordered list is the already ordered list of doublings.
def add_to_ordered_list(element, ordered_list):
	# Find the index to insert the element by iterating in reverse order
	for i, el in reversed(list(enumerate(ordered_list))):
		if element[0] >= el[0]:
			ordered_list.insert(i + 1, element)
			return(ordered_list)
	# If the element's number is smaller than all elements, insert it at the beginning
	ordered_list.insert(0, element)
	return(ordered_list)

#This function is to substract the first value to the whole list, it will return the new updated
#list of values
def subtract_time(number, lst):
	modified_lst = []
	for element in lst:
		time = element[0] - number
		modified_element = [time] + element[1:]
		modified_lst.append(modified_element)
	return(modified_lst)


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
		if(cell_data[1] in component):
			positions_to_remove.append(i)
	
	# Remove the entries from the list of lists in reverse order to avoid index issues
	for pos in reversed(positions_to_remove):
		list_cells.pop(pos)
	
	return(list_cells)


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
		# cells_to_divide=remove_cells_by_component(cells_to_divide, snowflake.nodes())
		cells_to_divide=[]

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

	snowflake.add_node(cont_ids,number_divisions=0, identity='parent') #adding the first cell
	#generations=1, cluster_id='m'
	cells_to_divide=[[sample_doub_t(first_div_params, second_div_params, 0),1,0, 'm']] #sampling doubling time of the first cell
	#cells to divide structure, [time to division, cell_id, num_generations, cluster_id]
	cont_ids+=1

	t_sim=0

	cont_updates=1

	cluster_fract_size=-1

	cont_gen=0

	# while reached_max_generations==False:
	while len(cells_to_divide)>0:

		#Part 1: Update division time and select cell to divide
		#check if the first element is 0, if not do the substraction of the first element to the whole 
		#list, it could be that the first element was 0 if there where several cells that were going
		#to divide at the same time
		if(cells_to_divide[0][0]!=0):
			t_sim+=cells_to_divide[0][0]
			cells_to_divide=subtract_time(cells_to_divide[0][0], cells_to_divide)

		#print(t_sim)
		# print(cells_to_divide)

		mother_id=cells_to_divide[0][1]

		temp_mother_identity=snowflake.nodes[mother_id]['identity']
		
		#Part 2: Create daughter cell and generate next sample times
		#create daughter cell
		daughter_id=cont_ids
		snowflake.add_node(daughter_id, number_divisions=0, identity=temp_mother_identity)
		snowflake.add_edge(mother_id, daughter_id)
		temp_daughter_time=sample_doub_t(first_div_params, second_div_params, 0) #using 0 instead of snowflake.nodes[daughter_id]["number_divisions"]
		cells_to_divide=add_to_ordered_list([temp_daughter_time, daughter_id], cells_to_divide)
		cont_ids+=1

		#mother cell division
		snowflake.nodes[mother_id]["number_divisions"]+=1 #adding that the mother divided one more time (because of the step before)

		temp_mother_div=snowflake.nodes[mother_id]["number_divisions"]
		temp_mother_time=sample_doub_t(first_div_params, second_div_params, temp_mother_div)
		cells_to_divide=add_to_ordered_list([temp_mother_time, mother_id], cells_to_divide)
		#print("mother next doub ",mother_id, t_sim+temp_mother_time)

		#removing first element of the list
		cells_to_divide.pop(0)


		#Part 3: save doubling information
		dict_files['sampled_times'].write(sim_number+','+str(temp_mother_div)+','+str(temp_mother_time)+'\n')
		dict_files['sampled_times'].write(sim_number+','+str(0)+','+str(temp_daughter_time)+'\n')
		#saving difference in doubling time
		dict_files['diff_d_t'].write(sim_number+","+str(cont_gen)+","+str(abs(temp_mother_time-temp_daughter_time))+"\n")

		#Part 4, calculate first order edge degree
		edge_to_remove, edge_degree=calculate_max_edge_degree(snowflake)

		if(edge_degree>=edge_degree_threshold):

			#find if the edge to remove is formed by the highest degree nodes
			if(input_variables['test_frag']):
				is_fractured_edge_between_highest_degree_nodes(snowflake, edge_to_remove, cont_gen)

			#Calling function to break clusters, update list of cells_to_divide, and save intermediate networks
			cells_to_divide=break_cluster(snowflake, cells_to_divide, cont_gen, cont_updates, edge_to_remove, num_generations)

			cont_gen+=1

		cont_updates+=1

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

'''


'''

#initializing variables
#Loading doubling time distributions
# dict_doub_t_dist=load_doubling_time(args.doubling_t, "minutes")

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
"first_div_params":first_div_params, "second_div_params":second_div_params, "cluster_to_keep":cluster_to_keep}
# "dict_doub_t_dist":dict_doub_t_dist


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
# log_file.write("-d: "+args.doubling_t+"\n")
log_file.write("-i: "+args.first_dist_params+"\n")
log_file.write("-j: "+args.second_dist_params+"\n")
log_file.write("-n: "+args.number_sims+"\n")
log_file.write("-o: "+args.output_dir+"\n")
log_file.write("-g: "+args.generations+"\n")
log_file.write("-e: "+args.edge_degree+"\n")
log_file.write("-c: "+args.cluster_to_keep+"\n")
log_file.write("\nExecution times:\n")

for i in range(1, num_iterations+1):

	# print(i)

	input_variables["sim_number"]=str(i) #changing variable type as it is only used as a character for saving results of the simulation

	#counting execution time
	start_time = time.time()

	#grow the cluster of this simulation
	simulate_one_cluster_growth(input_variables)

	end_time=time.time()
	elapsed_time=end_time-start_time


	log_file.write("Simulation "+str(i)+" completed in "+str(elapsed_time)+" seconds.\n")


log_file.write("\nProgram finished execution without any errors.")
log_file.close()

#closing all files from the dictionary
close_all_files(dict_files)
