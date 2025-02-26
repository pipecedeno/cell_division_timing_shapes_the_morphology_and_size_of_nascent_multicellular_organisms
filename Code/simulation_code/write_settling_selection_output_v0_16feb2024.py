#!/usr/bin/env python3


'''

inputs:
-nodes_file: file with the information of each of the cells, used to count how many cells of each strain survived
-output_file: name of the proportions file of settling selection
-strain_names: same variable used in the other scripts (format: strain1:strain2), it is used to avoid any errors if
	one of the populations go extinct
-current transfer
-surviving clusters: file were the ids of the clusters that survived settling selection is saved
-sim_number
-volumes_file: filename were all the volumes of the clusters are saved (this file has 2 columns ids and volume separated by a
	comma without a header). This file is used to know which clusters underwent settling selection as only 10% at random go to 
	settling selection
-processed_vols: file where the volumes of each strain is going to be saved and the flag if they survived settling 
	selection or not

outputs
-proportion_settling_selection_registry.csv
sim_number,transfer,strain1,strain2,cluster_pop1_b,cluster_pop2_b,cells_pop1_b,cells_pop2_b,total_clusters_b,total_cells_b,
cluster_pop1_a,cluster_pop2_a,cells_pop1_a,cells_pop2_a,total_clusters_a,total_cells_a
-volumes_settling_selection.csv
sim_number,strain,cluster_id,transfer,volume,number_cells,survived,initial_pos,final_pos,distance_travelled


Process the pandas by filtering and counting rows as cells and unique cluster ids as clusters

Note: if I decide to save the volumes before and after selection I could do that here, as I could have a file with id and volume
and just process the file to add more columns to the table, being one of those values a flag if they survived settling selectio

'''


import os
import pandas as pd
import argparse
# import matplotlib.pyplot as plt
# from random import choice
import numpy as np



parser = argparse.ArgumentParser()
parser.add_argument('-n','--nodes_file',dest="nodes_file",required=True) #nodes file
parser.add_argument('-o','--output_file',dest="output_file",required=True) #output file
parser.add_argument('-m','--strain_names',dest="strain_names",required=True) #name of the strains (format: strain1:strain2)
parser.add_argument('-r','--curr_selec_round',dest="curr_selec_round",required=True) #current selection round
parser.add_argument('-s','--surviving_file',dest="surviving_file",required=True) #surviving ids file
parser.add_argument('-i','--sim_number',dest="sim_number",required=True) #simulation number
parser.add_argument('-v','--volumes_file',dest="volumes_file",required=True) #volumes file
parser.add_argument('-p','--processed_vols',dest="processed_vols",required=True)
args = parser.parse_args()

#input variables
nodes_filename=args.nodes_file
surviving_filename=args.surviving_file
output_filename=args.output_file
curr_selec_round=args.curr_selec_round
sim_number=args.sim_number
volumes_filename=args.volumes_file
processed_volumes_out=args.processed_vols


#read surviving networks ids
file=open(os.path.join(surviving_filename), 'r')
list_ids=[int(line.rstrip('\n').split(',')[0]) for line in file.readlines()]
file.close()

# get ids of cluster that went to settling selection
file=open(volumes_filename, 'r')
settling_selection_ids=[int(line.rstrip('\n').split(',')[0]) for line in file.readlines()]
file.close()


#load nodes files
nodes_df=pd.read_csv(nodes_filename)

# strain_names=sorted(np.unique(nodes_df.strain))
# strain_1_name=strain_names[0]
# strain_2_name=strain_names[1]
# By receiving this variable it should avoid the error of one of the populations going extinct
strain_1_name=args.strain_names.split(':')[0]
strain_2_name=args.strain_names.split(':')[1]

#filtering out the ids of clusters that were not in the settling selection
nodes_df=nodes_df[nodes_df.cluster_id.isin(settling_selection_ids)]

# Processing volumes to create a csv file with the volume distribution and the flag if they survived settling selection
#sim_number,strain,transfer,volume,number_cells,survived
if(os.path.exists(processed_volumes_out)):
	output=open(processed_volumes_out, "a")
else:
	output=open(processed_volumes_out, "w")
	output.write('sim_number,strain,cluster_id,transfer,diameter,number_cells,survived,initial_pos,final_pos,distance_travelled,final_time,settling_speed\n')


file=open(volumes_filename, 'r')

for line in file.readlines():
	#processing lines from volumes file
	line_split=line.rstrip('\n').split(',')
	temp_id=line_split[0]
	temp_diameter=line_split[1]
	temp_ini_pos=line_split[2]
	temp_final_pos=line_split[3]
	temp_dist_travelled=line_split[4]
	temp_final_time=line_split[5]
	temp_speed=line_split[6]

	temp_nodes=nodes_df[nodes_df.cluster_id==int(temp_id)]

	temp_strain=np.unique(temp_nodes.strain)[0]
	temp_num_cells=len(temp_nodes.node_id)

	flag_survived=int(temp_id) in list_ids

	output.write(sim_number+','+temp_strain+','+temp_id+','+curr_selec_round+','+str(temp_diameter)+','+str(temp_num_cells)+','+str(flag_survived)+','+
		temp_ini_pos+','+temp_final_pos+','+temp_dist_travelled+','+temp_final_time+','+temp_speed+'\n')

output.close()

file.close()


# calculate values before settling selection
ini_cont_cells_strain1=np.sum(nodes_df.strain==strain_1_name)
ini_cont_cells_strain2=np.sum(nodes_df.strain==strain_2_name)
ini_total_cells=ini_cont_cells_strain1+ini_cont_cells_strain2

ini_cont_clusters_strain1=len(np.unique(nodes_df[nodes_df.strain==strain_1_name].cluster_id))
ini_cont_clusters_strain2=len(np.unique(nodes_df[nodes_df.strain==strain_2_name].cluster_id))
ini_total_clusters=ini_cont_clusters_strain1+ini_cont_clusters_strain2

# filtering the dataset to get clusters that survived settling selection
surviving_nodes=nodes_df[nodes_df.cluster_id.isin(list_ids)]

# calculate values after settling selection
fin_cont_cells_strain1=np.sum(surviving_nodes.strain==strain_1_name)
fin_cont_cells_strain2=np.sum(surviving_nodes.strain==strain_2_name)
fin_total_cells=fin_cont_cells_strain1+fin_cont_cells_strain2

if(fin_cont_cells_strain1==0 or fin_cont_cells_strain2==0):
	log_file=open('log.txt', 'a')
	if(fin_cont_cells_strain1==0):
		log_file.write('Population '+strain_1_name+' went extinct!\n')
	else:
		log_file.write('Population '+strain_2_name+' went extinct!\n')
	log_file.close()
	break_file=open('stop_extinction.txt','w')
	break_file.close()

fin_cont_clusters_strain1=len(np.unique(surviving_nodes[surviving_nodes.strain==strain_1_name].cluster_id))
fin_cont_clusters_strain2=len(np.unique(surviving_nodes[surviving_nodes.strain==strain_2_name].cluster_id))
fin_total_clusters=fin_cont_clusters_strain1+fin_cont_clusters_strain2


#initializing output file if it doesn't exist already
if(os.path.exists(output_filename)):
	output=open(output_filename, "a")
else:
	output=open(output_filename, "w")
	output.write('sim_number,transfer,strain1,strain2,'+
		'clusters_pop1_b,clusters_pop2_b,cells_pop1_b,cells_pop2_b,'+
		'total_clusters_b,total_cells_b,clusters_pop1_a,clusters_pop2_a,'+
		'cells_pop1_a,cells_pop2_a,total_clusters_a,total_cells_a\n')

#write output to file
output.write(str(sim_number)+','+str(curr_selec_round)+','+strain_1_name+','+strain_2_name+','+
		str(ini_cont_clusters_strain1)+','+str(ini_cont_clusters_strain2)+','+str(ini_cont_cells_strain1)+','+str(ini_cont_cells_strain2)+','+
		str(ini_total_clusters)+','+str(ini_total_cells)+','+str(fin_cont_clusters_strain1)+','+str(fin_cont_clusters_strain2)+','+
		str(fin_cont_cells_strain1)+','+str(fin_cont_cells_strain2)+','+str(fin_total_clusters)+','+str(fin_total_cells)+'\n')

output.close()


log_file=open('log.txt', 'a')
log_file.write('Settling selection files saved\n')
log_file.close()


