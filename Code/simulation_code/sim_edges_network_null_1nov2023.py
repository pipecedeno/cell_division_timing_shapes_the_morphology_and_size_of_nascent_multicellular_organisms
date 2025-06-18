#!/usr/bin/env python3
'''
Date: 2023jun2
Version 2: 1Nov2023

This version of the program is made to output a the nodes of the edges that were added each step, and now
the network is going to be simulated until it reaches a desired size. The edges output file is meant to be used
as input of the matlab codes to do the spatial and physics simulations.

###
Notes:
One of the biases that this program has when there are many cells that divide at the same time,
is that the first cell that appeared in the list with that time is the one that is going to divide
first, and also because of how the elements are added to the list, if the mother and the daughter
cell sample the same doubling time, the daughter is going to divide first because that is the cell
that gets added first to the list.
###



'''

import os
import pandas as pd
import argparse
import networkx as nx
from random import choice
import random
import numpy as np
import time


parser = argparse.ArgumentParser()
parser.add_argument("-m", "--max_num_cells",dest="max_num_cells",required=True) #max amount of cells that the simulation is going to have
parser.add_argument('-d','--doubling_t',dest="doubling_t",required=True) #csv file of the doubling time distributions
parser.add_argument('-n','--number_sims',dest="number_sims",required=True) #number of simulations it is going to run
parser.add_argument('-o','--output_dir',dest="output_dir",required=True) #output directory
args = parser.parse_args()



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
def sample_doub_t(dict_doub_t_dist, divisions_cell):
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
	return(choice(dict_doub_t_dist[sampling_dist]))



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
	return ordered_list

#This function is to substract the first value to the whole list, it will return the new updated
#list of values
def subtract_time(number, lst):
	modified_lst = []
	for element in lst:
		time = element[0] - number
		modified_element = [time, element[1]]
		modified_lst.append(modified_element)
	return(modified_lst)



# #Function to save the degree distribution of the network, the degree distribution is saved as a
# #probability distribution to be normalized for network size
# def save_degree_distribution_to_csv(network):
# 	# Get the degree of all nodes
# 	degrees = dict(network.degree())
	
# 	# Create a DataFrame from the degrees dictionary
# 	df_degrees = pd.DataFrame(degrees.items(), columns=['Node', 'Degree'])
	
# 	# Create a DataFrame for the degree distribution
# 	degree_counts=df_degrees['Degree'].value_counts().reset_index()
# 	degree_counts.columns=['Degree', 'Frequency']
	
# 	# Calculate the probability distribution
# 	total_nodes = len(df_degrees)
# 	degree_counts['Probability']=degree_counts['Frequency'] / total_nodes
	
# 	# Add static columns to the dataframe
# 	degree_counts['sim_number']=input_variables['sim_number']

# 	for i in range(len(degree_counts)):
# 		dict_files['degree_dist'].write(",".join([str(j) for j in degree_counts.iloc[i]])+"\n")



# #Function to save properties of the final network
# def save_final_networks(network):
# 	sim_number=input_variables["sim_number"]


# 	temp_diameter_network=str(nx.diameter(network))
# 	dict_files['diameter'].write(sim_number+","+temp_diameter_network+"\n")

# 	save_degree_distribution_to_csv(network)

# 	#calculate parent feature values of netlsd and netsimile
# 	temp_heat=netlsd.heat(network)
# 	#saving heat values of network
# 	dict_files['heat_val'].write(sim_number+","+",".join([str(i) for i in temp_heat])+"\n")

# 	temp_node_features=feature_extraction(network)
# 	temp_feat_vec=graph_signature(temp_node_features)
# 	#saving feature values
# 	dict_files['feature_val'].write(sim_number+","+",".join([str(i) for i in temp_feat_vec])+"\n")



def simulate_one_cluster_growth(input_variables, dict_doub_t_dist):

	#defining general use variables
	output_dir=input_variables['output_dir']
	sim_number=input_variables['sim_number']
	edges_dir=input_variables['edges_dir']

	cont_cells=1

	snowflake=nx.Graph() #creating graph

	snowflake.add_node(cont_cells,number_divisions=0) #adding the first cell
	cells_to_divide=[[sample_doub_t(dict_doub_t_dist, 0),1]] #sampling doubling time of the first cell
	cont_cells+=1

	t_sim=0

	output_file=open(os.path.join(edges_dir, 'edges_network_sim_'+str(sim_number)+'.csv'), 'w')
	output_file.write('node1,node2,time_sim\n')

	while cont_cells<max_num_cells:

		#Part 1: Update division time and select cell to divide
		#check if the first element is 0, if not do the substraction of the first element to the whole 
		#list, it could be that the first element was 0 if there where several cells that were going
		#to divide at the same time
		if(cells_to_divide[0][0]!=0):
			t_sim+=cells_to_divide[0][0]
			cells_to_divide=subtract_time(cells_to_divide[0][0], cells_to_divide)

		#print(t_sim)
		#print(cells_to_divide)

		#select the first element of the list of cells to divide, which should be 0 and do the operation
		#to simulate the division of that cell

		mother_id=cells_to_divide[0][1]
		
		#Part 2: Create daughter cell and generate next sample times
		#create daughter cell
		daughter_id=cont_cells
		snowflake.add_node(daughter_id, number_divisions=0)
		snowflake.add_edge(mother_id, daughter_id)
		temp_daughter_time=sample_doub_t(dict_doub_t_dist, 0) #using 0 instead of snowflake.nodes[daughter_id]["number_divisions"]
		cells_to_divide=add_to_ordered_list([temp_daughter_time, daughter_id], cells_to_divide)
		cont_cells+=1

		#mother cell division
		snowflake.nodes[mother_id]["number_divisions"]+=1 #adding that the mother divided one more time (because of the step before)

		temp_mother_div=snowflake.nodes[mother_id]["number_divisions"]
		temp_mother_time=sample_doub_t(dict_doub_t_dist, temp_mother_div)
		cells_to_divide=add_to_ordered_list([temp_mother_time, mother_id], cells_to_divide)
		#print("mother next doub ",mother_id, t_sim+temp_mother_time)

		#removing first element of the list
		cells_to_divide.pop(0)

		#Part 3: save doubling information
		#saving difference in doubling time
		dict_files['diff_d_t'].write(sim_number+","+str(abs(temp_mother_time-temp_daughter_time))+"\n")

		dict_files['sampled_times'].write(sim_number+','+str(temp_mother_div)+','+str(temp_mother_time)+'\n')
		dict_files['sampled_times'].write(sim_number+','+str(0)+','+str(temp_daughter_time)+'\n')

		#save added edge to the output file
		output_file.write(str(mother_id)+','+str(daughter_id)+','+str(t_sim)+'\n')


	# save_final_networks(snowflake)

	output_file.close()


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
		dict_files['diff_d_t'].write('sim_number,diff_minutes\n')

	sampled_times_file=os.path.join(output_dir,"sampled_times.csv")
	if(os.path.exists(sampled_times_file)):
		dict_files['sampled_times']=open(sampled_times_file, "a")
	else:
		dict_files['sampled_times']=open(sampled_times_file, "w")
		dict_files['sampled_times'].write("sim_number,number_divisions,minutes"+"\n")



	# diam_file=os.path.join(output_dir, "networks_diameter.csv")
	# if(os.path.exists(diam_file)):
	# 	dict_files['diameter']=open(diam_file, "a")
	# else:
	# 	dict_files['diameter']=open(diam_file, "w")
	# 	dict_files['diameter'].write("sim_number,diameter\n")

	# degree_dist_file=os.path.join(output_dir,"degree_distribution.csv")
	# if(os.path.exists(degree_dist_file)):
	# 	dict_files['degree_dist']=open(degree_dist_file, "a")
	# else:
	# 	dict_files['degree_dist']=open(degree_dist_file, "w")
	# 	dict_files['degree_dist'].write("Degree,Frequency,Probability,sim_number\n")

	# heat_values_file=os.path.join(output_dir,"heat_values.csv")
	# if(os.path.exists(heat_values_file)):
	# 	dict_files['heat_val']=open(heat_values_file, "a")
	# else:
	# 	dict_files['heat_val']=open(heat_values_file, "w")
	# 	dict_files['heat_val'].write("sim_number,"+",".join(["val_"+str(i) for i in range(1, 251)])+"\n")

	# feat_values_file=os.path.join(output_dir,"feature_values.csv")
	# if(os.path.exists(feat_values_file)):
	# 	dict_files['feature_val']=open(feat_values_file, "a")
	# else:
	# 	dict_files['feature_val']=open(feat_values_file, "w")
	# 	dict_files['feature_val'].write("sim_number,"+",".join(["val_"+str(i) for i in range(1, 36)])+"\n")


	return(dict_files)


#function to close all files opened in the dictionary
def close_all_files(dict_files):
	for file_obj in dict_files.values():
		file_obj.close()



### MAIN ###

#note: I am getting rid of the division_registry and doub_t_sampled as the files are not 
#useful anymore

#initializing variables
#Loading doubling time distributions
dict_doub_t_dist=load_doubling_time(args.doubling_t, "minutes")


num_iterations=int(args.number_sims)

max_num_cells=int(args.max_num_cells)


output_dir=args.output_dir

#creating output directory if it doesn't exist, and stopping the program execution if the 
#directory exists
if(not os.path.exists(output_dir)):
	os.mkdir(output_dir)
else:
	print("Output directory already exists, stopping program execution to avoid overwriting data")
	exit()
edges_dir=os.path.join(output_dir, 'edges_sim_files')
if(not os.path.exists(edges_dir)):
	os.mkdir(edges_dir)

input_variables={'dict_doub_t_dist':dict_doub_t_dist, 'max_num_cells':max_num_cells, 
'output_dir':output_dir, 'edges_dir':edges_dir}


#Initializing all output files. Note: this variable is accesed as a global variable
dict_files=initialize_files(input_variables)



#creating and writing information of the program in the log file
log_file=open(os.path.join(output_dir,"log.txt"), "a")
log_file.write("sim_edges_network_null_1nov2023.py\nInputs received:\n")
log_file.write("-d: "+args.doubling_t+"\n")
log_file.write("-n: "+args.number_sims+"\n")
log_file.write("-o: "+args.output_dir+"\n")
log_file.write("-m: "+args.max_num_cells+"\n")
log_file.write("\nExecution times:\n")


for i in range(1, num_iterations+1):

	input_variables['sim_number']=str(i)

	#counting execution time
	start_time = time.time()

	#grow the cluster of this simulation
	simulate_one_cluster_growth(input_variables, dict_doub_t_dist)

	end_time=time.time()
	elapsed_time=end_time-start_time


	log_file.write("Simulation "+str(i)+" completed in "+str(elapsed_time)+" seconds.\n")


log_file.write("\nProgram finished execution without any errors.")
log_file.close()

#closing all files from the dictionary
close_all_files(dict_files)

