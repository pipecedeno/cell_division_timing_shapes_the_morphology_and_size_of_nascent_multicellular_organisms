#!/usr/bin/python3

'''
Date: 5 May 2025
This is an updated version of the script merge_time_lapse_tables_v2.py where the cell_death and big_clust
columns have been removed as they are no longer necessary in the final version.

Expected file naming format:
(date)_(strain)_(replicate)_(sufix).csv
	-date: format is not important but make it consistant accross files
	-strain: can be the gob id or PA replicate population + time point separated by a "-" (example:
	PA1_t600)
	-replicate: just a number to differentiate time lapses from the same attempt
	-sufix: one of the following: "branch_table" "spots_table" "edges_table"

input:
directories: directory with all the csv file, for each timelapse 3 files are expected and should have the following sufixes
for them to be identified:
	"_branch_table.csv"
	"_spots_table.csv"
	"_edges_table.csv"
	Note: if needed this can be modified.

output:
-o: csv file merging all of the important information.

The columns that the final csv file will have, as of the current version, are the following:
	-From column A to I it has the information of the branch table that is obtained of the branch hierarchy analysis.
	-best_ellipse_angle: it has the angle of the last cell of the branch which needs to have a solidity bigger than 0.95, and it being
	different from 1, so it's not a circle of the spot I manually added. The angle is already in a range from 0 to 90.
	-best_ellipse_id: is the spot id from which the angle is used in best_ellipse_angle.
	-From column L to AN it has the columns of the spots table, except for the label and track_id column. The row of the spots table that is added
	to the row of the final table, is the spot corresponding to the last frame before a bud is detected.
	-date: date that the time lapse was made (format:(year)(three first letters of month)(day))
	-strain
	-strain_date
	-complete_id: id of each timelapse the branch is from, format: (strain)_(timepoint)_(replicate)_(date)
	-strain_timepoint: (strain)_(timepoint). The strain can be PA1 to 5, or gob8 and gob21, for gob8 and gob21 the time point is ancestor
	-cell_id: id to track specific cell, so when a branch divides, the mother cell will have the same cell_id and the new bud will get a new cell_id
	-division_number: number starting from 0 that means how many times a cell has divided, 0 if it hasn't divided yet. [Note: need to consider if this
	number should start from 1, as for the analysis it's more intuitive if it starts from 1]
	-mother_cell_id: it's the id of the mother cell plus the division number, so that the mother and daughter relationship can be tracked outisde from this code,
	and the division number will help find the specific branch.
	-budd_angle: angle that the budd is appearing from mother cell, it is only calculated if best_ellipse_angle is different from NA for both the mother and
	daughter cells. If it wasn't calculated then it will have an NA.
	-mother_delta_t: if a cell divided it will have a mother_delta_t in seconds, which is the amount of time the mother cell took to divide.
	-max_num_div: maximum amount of times that the cell divided, this number is repeated, as the rows of the same cell will have the same value, so if you want unique 
	values should get only one of this values per cell_id.

Modifications history:
Date: 26jan2023
Added the names of the t600 strains into the strain_convertion_dict and time_point_dict so that the information of those time lapses can also be processed.

Date: 1mar2023
Added a dictionary so that it counts how many time lapses are being analyzed by strain and by timepoint, so at the end of the program it will
print that information so that it is easier to build the table of the amount of time lapses analyzed

Date: 5apr2023
Change the name of gob21 and gob8 in the dictionary "strain_convertion_dict" so that they are now refered to as petite and grande, respectively.
The idea of this change is so that the strains are presented as petite ancestor and grande ancestor, and in the plots the fridge names of the strains 
are no longer showed.
Also, added the column id_file so that it is easier to find the file of each specific time lapse as their ids in the other columns are no longer the same
as the file name.

Date: 10apr2023
Changed gob21 and gob8 timepoint to be called 't0' instead of being called ancestor
small_id variable was deleted, as 2 time lapses had the same small_id because in the 22nov2022 t200 and t1000 were imaged together.

Date: 5 May 2025
Removed the dict_cell_death and list_clust_big variables and associated columns as they are no longer needed.

#just the command that I use to run this script, and in which directory, can be used as an example of how to run this code
cd work_dir/timelapses/
merge_time_lapse_tables_v3.py -d timelapses_csv_files/ -o tl_merged_csv_doubling.csv
'''

import pandas as pd
from glob import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output",dest="output",required=True) #name of output csv file
parser.add_argument('-d','--directories',dest="directories",required=True) #directory with all of the csv files
args = parser.parse_args()

#this will generate the new id that is going to be used for the cell_id column
#this function receives the first_part_id which should be in the folowing format:
#(strain)_(replicate)_(day)_(# of cell) [the number of cell will always have 4 numbers]
#and will return the complete id and the counter increased by 1
#example: input: 'gob21_',100 output: 'gob21_0101',101
def generate_new_cell_id(first_part_id, cont_assigned_cells):
	return(first_part_id+str(cont_assigned_cells).zfill(4), cont_assigned_cells+1)


def get_mother_daughter_relationship(merged, edges, spots):
	#dictionary to count how many times they have divided
	dict_num_div={}

	#getting info for the id cell generator
	#first_part_id=merged['small_id'][0]+'_'
	first_part_id=merged['complete_id'][0]+'_'
	cont_assigned_cells=1

	#sorting merged table by 'POSITION_T'
	merged=merged.sort_values('POSITION_T').reset_index(drop=True)

	#create new columns
	merged['cell_id']='NA'
	merged['division_number']='NA'
	merged['mother_cell_id']='NA'
	merged['budd_angle']='NA'
	#merged['cell_id_w_div']='NA'
	#merged['new_budd']=0 #0=False, 1=True
	merged['mother_delta_t']='NA'

	for i in range(len(merged)):
		#if it doesn't has any predecessor then it means that is the first cell and a new id gets assigned
		if(merged.N_PREDECESSORS[i]==0):
			new_id,cont_assigned_cells=generate_new_cell_id(first_part_id, cont_assigned_cells)
			merged.at[i,'cell_id']=new_id
			merged.at[i, 'division_number']=0
			dict_num_div[new_id]=0
			merged.at[i,'mother_cell_id']=merged.cell_id[i]+'_'+str(0)

		#only if it has daughter cells it will do the following:
		#it will give the same cell id for it self, which is the cell with the bigger size, and for the smaller
		#cell, which is the budded cell, will give a new cell id
		if(merged.N_SUCCESSORS[i]==2):
			temp_last_id=merged.LAST[i]

			#how to get mother cell info
			#temp_mother_area=merged.AREA[i]
			temp_mother_doubling=merged['DELTA_T'][i] #converting into hours

			#getting daughter cells id
			list_daughters=edges.loc[edges['SPOT_SOURCE_ID'] == temp_last_id]['SPOT_TARGET_ID'].reset_index(drop=True)

			#getting area of daughter cells
			# area_d_1=int(spots.loc[spots.ID==list_daughters[0]]['AREA'])
			# area_d_2=int(spots.loc[spots.ID==list_daughters[1]]['AREA'])
			area_d_1 = int(spots.loc[spots.ID==list_daughters[0]]['AREA'].iloc[0])
			area_d_2 = int(spots.loc[spots.ID==list_daughters[1]]['AREA'].iloc[0])

			#getting daugther cell index in merged table
			daughter_1_index=merged.index[merged['FIRST']==list_daughters[0]].values[0]
			daughter_2_index=merged.index[merged['FIRST']==list_daughters[1]].values[0]

			mother_id=merged.cell_id[i] #if the cell doesn't has a mother cell then this is going to get NA as value, and it will
			#add that value to the cell, so if the cell is the first cell the mother_cell_id is going to be 0 
			daughter_id,cont_assigned_cells=generate_new_cell_id(first_part_id,cont_assigned_cells)

			#adding to number of divisions dictionary
			dict_num_div[mother_id]+=1
			dict_num_div[daughter_id]=0
			temp_mother_id_div=mother_id+'_'+str(dict_num_div[mother_id])

			#checking which cell is the largest and assigning the id accordingly
			if(area_d_1>area_d_2):
				#assigning cell id of mother cell
				merged.loc[daughter_1_index,'cell_id']=mother_id
				merged.loc[daughter_1_index,'division_number']=dict_num_div[mother_id]
				merged.loc[daughter_1_index, 'mother_cell_id']=temp_mother_id_div
				merged.loc[daughter_1_index, 'mother_delta_t']=temp_mother_doubling

				#getting angle at which mother cell is budding
				mother_angle=merged['best_ellipse_angle'][i]
				daughter_angle=merged['best_ellipse_angle'][daughter_2_index]
				if(mother_angle!='NA' and daughter_angle!='NA'):
					angle_diff=abs(float(daughter_angle)-float(mother_angle))
					if(angle_diff>90):
						angle_diff=180-angle_diff
					merged.loc[daughter_1_index,'budd_angle']=angle_diff

				#assigning new cell id
				merged.loc[daughter_2_index,'cell_id']=daughter_id
				merged.loc[daughter_2_index,'division_number']=0
				merged.loc[daughter_2_index,'mother_cell_id']=temp_mother_id_div
				#merged.loc[daughter_2_index,'new_budd']=1
				merged.loc[daughter_2_index, 'mother_delta_t']=temp_mother_doubling

			else:
				#assigning cell id of mother cell
				merged.loc[daughter_2_index,'cell_id']=mother_id
				merged.loc[daughter_2_index,'division_number']=dict_num_div[mother_id]
				merged.loc[daughter_2_index,'mother_cell_id']=temp_mother_id_div
				merged.loc[daughter_2_index, 'mother_delta_t']=temp_mother_doubling

				#getting angle at which mother cell is doubling
				mother_angle=merged['best_ellipse_angle'][i]
				daughter_angle=merged['best_ellipse_angle'][daughter_1_index]
				if(mother_angle!='NA' and daughter_angle!='NA'):
					angle_diff=abs(float(daughter_angle)-float(mother_angle))
					if(angle_diff>90):
						angle_diff=180-angle_diff
					merged.loc[daughter_1_index,'budd_angle']=angle_diff

				#assigning new cell id
				merged.loc[daughter_1_index,'cell_id']=daughter_id
				merged.loc[daughter_1_index,'division_number']=0
				merged.loc[daughter_1_index,'mother_cell_id']=temp_mother_id_div
				#merged.loc[daughter_1_index,'new_budd']=1
				merged.loc[daughter_1_index, 'mother_delta_t']=temp_mother_doubling



	
	merged['max_num_div']='NA'
	merged['max_num_div']=[dict_num_div[j] for j in merged['cell_id']]


	return(merged)


#This funcion is used to fill the information for the columns best_ellipse_angle and best_ellipse_id
#this value will later be used to obtain the budding angle. 
def get_usable_angle(branches, spots, edges):

	#transform angle to degrees, and add 180 if its a negative value
	spots['angle_transformed']=[i if i>=0 else i+180 for i in spots['ELLIPSE_THETA']*57.2958]

	#creating smaller tables to be easier to use
	small_edges=edges[['SPOT_SOURCE_ID','SPOT_TARGET_ID']]
	small_spots=spots[['ID','angle_transformed','SOLIDITY', 'ELLIPSE_ASPECTRATIO']]

	#creating variables to fill
	branches['best_ellipse_angle']='NA'
	branches['best_ellipse_id']='NA'

	for i in range(len(branches)):

		temp_first_id=temp_df_branches['FIRST'][i]
		temp_last_id=temp_df_branches['LAST'][i]

		#is the branch is only one spot, it is going to be ignored
		if(temp_first_id==temp_last_id):
			continue

		#do while to get the list of edges of a branch with the information of that row.
		list_edges=[]
		list_edges.append(temp_first_id)
		temp_edge=small_edges[small_edges['SPOT_SOURCE_ID']==list_edges[-1]]
		while(temp_edge['SPOT_TARGET_ID'].values.item()!=temp_last_id):
			list_edges.append(temp_edge['SPOT_TARGET_ID'].values.item())
			temp_edge=small_edges[small_edges['SPOT_SOURCE_ID']==list_edges[-1]]
		list_edges.append(temp_edge['SPOT_TARGET_ID'].values.item())

		#iterating through the list of edges to get the best angle to be saved, the best value is a value that has a solidity higher than 0.95, that is not a 
		#spot that I manually added, and that is the latest point possible that satisfy the other 2 conditions.
		assigned=False
		for edge in list_edges:
			temp_row=small_spots[small_spots['ID']==edge]
			#print(temp_row)
			if(temp_row['SOLIDITY'].values>0.95 and temp_row['ELLIPSE_ASPECTRATIO'].values!=1): #different from 1 because if they have a 1 then they are the circles I added to the timelapse
				#19August2024: changed temp_row['SOLIDITY'].values!=1 to temp_row['ELLIPSE_ASPECTRATIO'].values!=1 as here I should have used aspect ratio since the beginning
				temp_value_to_save=temp_row[['ID','angle_transformed']]
				assigned=True

		#if there's a point that satisy the conditions then it's saved.
		if(assigned):
			branches.loc[i,'best_ellipse_id']=temp_value_to_save['ID'].values.item()
			branches.loc[i,'best_ellipse_angle']=temp_value_to_save['angle_transformed'].values.item()
	
	#returns the branches table that was modified
	return(branches)



#getting the ids that are in the directory (id is date_strain_rep-number)
#with glob it retrieves the files that have the branch_table.csv sufix, and then with the split and the replace, the id of each 
#file is extracted
if(args.directories[-1]=='/'):
	directory=args.directories
else:
	directory=args.directories+'/'
list_files=[file.split("/")[-1].replace('_branch_table.csv','') for file in glob(directory+'*branch_table.csv')]


print('Amount of groups found: ',len(list_files))

#to print the files found if necessary
# list_files.sort()
# print(list_files)





#dictionary to change the names from gob2149-53 to PA1-5
strain_convertion_dict={"gob21":"petite","gob2149":"PA1","gob2150":"PA2","gob2151":"PA3","gob2152":"PA4","gob2153":"PA5",
"gob8":"grande","gob440":"PA3","gob386":"PA5","gob463":"PA1","gob385":"PA4","gob383":"PA2",
"gob1405":"PA1","gob1407":"PA2","gob1409":"PA3","gob1413":"PA4","gob2350":"PA5"}

#dictionary to add the timepoint to each strain, for example, so that it links the gob2149-53 as t1000s
time_point_dict={"gob21":"t0","gob2149":"t1000","gob2150":"t1000","gob2151":"t1000","gob2152":"t1000","gob2153":"t1000",
"gob8":"t0","gob440":"t200","gob386":"t200","gob463":"t200","gob385":"t200","gob383":"t200",
"gob1405":"t600","gob1407":"t600","gob1409":"t600","gob1413":"t600","gob2350":"t600"}


#flag to know if the variable (temp_big_clust) has already a table in it or not, if it has a value then the new table is 
#concatenated as this variable will hold the final table that is going to be used
assigned=False

#variables just to count the progress of the program
total_files=len(list_files)
cont_file_done=0

dict_cont_timelapses_per_strain={}

for id_file in list_files:
	# print(id_file)

	#creating variables with information that is going to be added as columns
	temp_date=id_file.split("_")[0]

	#modified section for version 2
	if("-" in id_file.split("_")[1]):
		temp_strain=id_file.split("_")[1].split("-")[0].upper()
		temp_time_point=id_file.split("_")[1].split("-")[1]
		temp_replicate=id_file.split("_")[2]
	else:
		temp_strain=strain_convertion_dict[id_file.split("_")[1]]
		temp_time_point=time_point_dict[id_file.split("_")[1]]
		temp_replicate=id_file.split("_")[2]


	#loading tables that were the output of the time lapse analysis
	temp_df_branches=pd.read_csv(directory+id_file+"_branch_table.csv", header=0, skiprows=[1,2,3])
	temp_df_spots=pd.read_csv(directory+id_file+"_spots_table.csv", header=0, skiprows=[1,2,3])
	temp_df_edges=pd.read_csv(directory+id_file+"_edges_table.csv", header=0, skiprows=[1,2,3])

	#get usable angle value for the branch table
	temp_df_branches=get_usable_angle(temp_df_branches,temp_df_spots,temp_df_edges)

	#merging branches table and spots table, this way the spot information of the branch will be the spot of the frame before a bud was detected and the cell
	#divided
	temp_merged=temp_df_branches.merge(temp_df_spots.drop(['LABEL','TRACK_ID','angle_transformed'], axis=1), left_on='LAST', right_on='ID', how='left')

	#assigning some values to the new columns
	temp_merged['date']=temp_date
	temp_merged['strain']=temp_strain
	temp_merged['strain_date']=temp_strain+'_'+temp_date #added 27oct2022
	temp_merged['complete_id']=temp_strain+'_'+temp_time_point+'_'+temp_replicate+'_'+temp_date #updated no longer id_file, 
	#now is (strain)_(timepoint)_(replicate)_(date)
	temp_date=id_file.split('_')[0]
	#temp_merged['small_id']=temp_strain+'_'+temp_time_point+'_'+temp_replicate #small id format: (strain)_(timepoint)_(replicate) #modified for v2

	temp_merged['timepoint']=temp_time_point

	#new column added for v2
	temp_merged['strain_timepoint']=temp_strain+'_'+temp_time_point
	temp_merged['strain_timepoint_date']=temp_strain+'_'+temp_time_point+'_'+temp_date
	temp_merged['id_file']=id_file

	#Deleted columns for v2
	# temp_merged['name']=strain_convertion_dict[temp_strain]
	# temp_merged['name_timepoint']=strain_convertion_dict[temp_strain]+'_'+temp_time_point

	#getting cell ids
	temp_merged=get_mother_daughter_relationship(temp_merged, temp_df_edges, temp_df_spots)

	#concatenating the dataframes
	if(not assigned):
		final_df=temp_merged
		assigned=True
	else:
		final_df=pd.concat([final_df,temp_merged], ignore_index=True)

	#counting amount of time lapses analyzed per strain timepoint
	if(temp_strain+'_'+temp_time_point in dict_cont_timelapses_per_strain):
		dict_cont_timelapses_per_strain[temp_strain+'_'+temp_time_point]+=1
	else:
		dict_cont_timelapses_per_strain[temp_strain+'_'+temp_time_point]=1

	cont_file_done+=1
	print(str(cont_file_done)+'/'+str(total_files))

#saving the final dataframe as a csv
final_df.to_csv(args.output, sep=',', index=False)

#printing how many time lapses are per strain and per timepoint
cont_total=0
for key in sorted(dict_cont_timelapses_per_strain.keys()):
	print(key+": "+str(dict_cont_timelapses_per_strain[key]))
	cont_total+=dict_cont_timelapses_per_strain[key]

print("Total: "+str(cont_total))

