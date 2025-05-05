#!/usr/bin/python3

'''
Date:23jan2023
The idea of this code is to be used to obtain the difference between the next two doubling times after a cell divisions, so 
for example the difference between the first doubling time of the bud and the second time the mother cell divided. In the output
file the divisions that are reported are of cells that have 2 succesors and the cases where the two daughter cells also have 2 
succesors, and this is because if those cells didn't had 2 sucessors then the difference can be wrongly calculated as there are
some cells in the time lapse that are not tracked entirely or that are lost after they divide, so that measurement will be completely
wrong.
The difference from the version 2 of the code to the version 1 is that this code was modified to save in the strain column from which 
replicate population are instead of saving the gob name, specially because for the t400 I don't know the gob names of the isolates.
Version 3: It added columns with the information of the mother and the daughter as it calculates which spot is the mother after the 
division of a branch based of the which branch has first spot closer to the mother. This should be highly accurate as when I was correcting
the time lapses I added the spots where the bud started or trying to keep the mother in the center.
So now the program loads the 3 tables obtained from the time lapses, the edges, branch and spots table.

Expected file naming format:
(date)_(strain)_(replicate)_(sufix).csv
	-date: format is not important but make it consistant accross files
	-strain: can be the gob id or PA replicate population + time point separated by a "-" (example:
	PA1_t600)
	-replicate: just a number to differentiate time lapses from the same attempt
	-sufix: one of the following: "branch_table" "edges_table"

input:
-d: directory with the csv files of the timelapses, the files that are going to be used are:
	-the output csv file of the branch hierarchy analysis
	-the output cvs file of the edges table

output:
-o: csv file with the following columns:
	-LABEL: from banches table	
	-TRACK_ID: from banches table	
	-N_PREDECESSORS: from banches table
	-N_SUCCESSORS: from banches table
	-minutes: amount of minutes that the branch took to divide
	-FIRST: from banches table
	-LAST: from banches table
	-d1_id: id of the first daughter cell
	-d1_minutes: minutes it took the first daughter cell to divide
	-d2_id: id of the second daughter cell
	-d2_minutes: minutes it took the second daughter cell to divide
	-diff_minutes: absolute difference of minutes it took a cell to divide after the other (abs(d1_minutes-d2_minutes))
	-date: date of the time lapse
	-strain: strain of the timelapse
	-complete_id: id of each timelapse the branch is from, format: (date)_(strain)_(replicate)
	-name: the name of the strains are the replicate population from which the isolates are made, so it can be PA1 to 5, or gob21 or gob8
	-timepoint: it is the amount of transfers of the isolate, so it can be t0, t200, t400, etc.
	-name_timepoint: is the combination of the columns name and timepoint separated by a "_"


Modification history:
Date: 22Mar2023
Changed the variables strain_convertion_dict and time_point_dict to include the information of the t600 strains also, as they are in the 
merge_time_lapse_tables_v2.py code

Date: 5apr2023
Change the name of gob21 and gob8 in the dictionary "strain_convertion_dict" so that they are now refered to as petite and grande, respectively.
The idea of this change is so that the strains are presented as petite ancestor and grande ancestor, and in the plots the fridge names of the strains 
are no longer showed.
Also, added the column id_file so that it is easier to find the file of each specific time lapse as their ids in the other columns are no longer the same
as the file name.

Date: 10apr2023
Changed gob21 and gob8 timepoint to be called 't0' instead of being called ancestor
small_id variable was deleted, as 2 time lapses had the same small_id because in the 22nov2022 t200 and t1000 were imaged together.

Date: 31Jan2024
Added a column to save the mother time and daughter time based to the distance of the closest point to the mother's last spot
Now it will also load the spot table as this one is going to be used to get the XY position of the spots

cd work_dir/timelapses/
./tl_cell_sync_estimate_v2.py -d timelapses_csv_files/ -o test_cell_sync_info.csv
'''

import pandas as pd
from glob import glob
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output",dest="output",required=True) #name of output csv file
parser.add_argument('-d','--directories',dest="directories",required=True) #directory with all of the csv files
args = parser.parse_args()


import math


#Function to calculate the euclidean distance between 2 points
def euclidean_distance(point1, point2):
	
	x1, y1 = point1
	x2, y2 = point2
	
	# x1=float(x1)
	# y1=float(y1)

	# x2=float(x2)
	# y2=float(y2)

	x1 = float(x1.iloc[0])
	y1 = float(y1.iloc[0])
	x2 = float(x2.iloc[0])
	y2 = float(y2.iloc[0])

	distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	return(distance)

#getting the ids that are in the directory (id is date_strain_rep-number)
#with glob it retrieves the files that have the branch_table.csv sufix, and then with the split and the replace, the id of each 
#file is extracted
if(args.directories[-1]=='/'):
	directory=args.directories
else:
	directory=args.directories+'/'
list_files=[file.split("/")[-1].replace('_branch_table.csv','') for file in glob(directory+'*branch_table.csv')]

print('Amount of groups found: ',len(list_files))


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

for id_file in list_files:
	#print(id_file)

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
	print(temp_strain)

	#loading files
	temp_df_branches=pd.read_csv(directory+id_file+"_branch_table.csv", header=0, skiprows=[1,2,3])
	temp_df_edges=pd.read_csv(directory+id_file+"_edges_table.csv", header=0, skiprows=[1,2,3])

	temp_df_spots=pd.read_csv(directory+id_file+"_spots_table.csv", header=0, skiprows=[1,2,3])

	#transforming the column of how many frames a branch is in the video before dividing into minutes
	temp_df_branches['minutes']=temp_df_branches['DELTA_T']/60
	
	#create small tables that are going to be used
	small_edges=temp_df_edges[['SPOT_SOURCE_ID','SPOT_TARGET_ID']]
	small_branches=temp_df_branches[['LABEL','TRACK_ID','N_PREDECESSORS','N_SUCCESSORS','minutes','FIRST','LAST']]

	#doing a left join on small edges so that now it has the information of the branch which first id is the same as SPOT_TARGET_ID
	left_join_div=small_edges.merge(small_branches.drop(['LABEL','TRACK_ID','N_PREDECESSORS','LAST'],axis=1), left_on='SPOT_TARGET_ID', right_on='FIRST', how='left').dropna(axis=0).reset_index(drop=True)

	#creating a dictionary with the id of the last spot of the branch as a key
	dividing_branches={}
	for i in range(len(small_branches)):
		if(small_branches.N_SUCCESSORS[i]==2):
			dividing_branches[small_branches.LAST[i]]=[]

	#looking for the branches that are connected to id of the last branch (looking for the id of the daughter and the id 
	#that the mother cell will keep after dividing)
	#To know save a value the following conditions need to be met:
	#1.- the SPOT_SOURCE_ID needs to be a key of the dictionary, meaning that the branch divides, and
	#2.- the daughter branch needs to also have to daughter, if this condition is met then the estimate of
	#the difference in doubling time will be incorrect because we are using a cell that never divided.
	#If the conditions are met, then a list will be appended to the list of the dictionary, this list will have the
	#SPOT_TARGET_ID and the doubling time in minutes.
	for i in range(len(left_join_div)):
		temp_key=left_join_div.SPOT_SOURCE_ID[i]
		if(temp_key in dividing_branches and left_join_div.N_SUCCESSORS[i]==2):
			dividing_branches[temp_key].append([left_join_div.SPOT_TARGET_ID[i],left_join_div.minutes[i]])

	#here is where the results will be added the table that is going to be the output.
	temp_final_df=small_branches[small_branches['N_SUCCESSORS']==2].reset_index(drop=True)
	temp_final_df['d1_id']='NA'
	temp_final_df['d1_minutes']='NA'
	temp_final_df['d2_id']='NA'
	temp_final_df['d2_minutes']='NA'
	temp_final_df['diff_minutes']='NA'

	temp_final_df['mother_id']='NA'
	temp_final_df['mother_minutes']='NA'
	temp_final_df['daughter_id']='NA'
	temp_final_df['daughter_minutes']='NA'

	#iterating through the table
	for i in range(len(temp_final_df)):
		#this if checks that the branch actually divided (This is not necessary because the only rows that were saved into this table
		#were the ones were this condition was satisfied)
		if(temp_final_df.N_SUCCESSORS[i]==2):
			temp_key=temp_final_df.LAST[i] #getting the key of the last spot
			#print(dividing_branches[temp_key])

			#checking if the branch actually has 2 daughters that divided, if one or non of them divided then the lenght of the 
			#list that corresponds to the id is going to be less than 2, and if it's more than two is won't be processed as it 
			#should be an error or a branch that was not correctly manually corrected
			if(len(dividing_branches[temp_key])==2):
				daughter_1, daughter_2=dividing_branches[temp_key]
				#print(daughter_1,daughter_2)
				temp_final_df.loc[i,'d1_id']=daughter_1[0]
				temp_final_df.loc[i,'d1_minutes']=round(daughter_1[1])
				temp_final_df.loc[i,'d2_id']=daughter_2[0]
				temp_final_df.loc[i,'d2_minutes']=round(daughter_2[1])
				#print(temp_final_df['d1_minutes'][i],temp_final_df['d2_minutes'][i])
				temp_final_df.loc[i,'diff_minutes']=abs(temp_final_df['d1_minutes'][i]-temp_final_df['d2_minutes'][i]) #getting the absolute difference
				#this step is done here and not at the end for all the column because some of the values in the column will have NAs


				# Calculating which cell is the mother and which cell is the daughter
				last_known_mother_info=temp_df_spots[temp_df_spots.ID==temp_final_df.LAST[i]]
				last_mother_XY=[last_known_mother_info.POSITION_X, last_known_mother_info.POSITION_Y]

				cell_1_info=temp_df_spots[temp_df_spots.ID==daughter_1[0]]
				cell_1_XY=[cell_1_info.POSITION_X, cell_1_info.POSITION_Y]

				cell_2_info=temp_df_spots[temp_df_spots.ID==daughter_2[0]]
				cell_2_XY=[cell_2_info.POSITION_X, cell_2_info.POSITION_Y]

				#Calculating euclidian distances from the last known position of the mother to each of the other spots
				dist_mother_to_1=euclidean_distance(last_mother_XY, cell_1_XY)
				dist_mother_to_2=euclidean_distance(last_mother_XY, cell_2_XY)

				if(dist_mother_to_1 < dist_mother_to_2):
					#Case where cell 1 is the mother
					temp_final_df.loc[i, 'mother_id']=daughter_1[0]
					temp_final_df.loc[i, 'mother_minutes']=round(daughter_1[1])

					temp_final_df.loc[i, 'daughter_id']=daughter_2[0]
					temp_final_df.loc[i, 'daughter_minutes']=round(daughter_2[1])

				else:
					#Case where cell 2 is the mother
					temp_final_df.loc[i, 'mother_id']=daughter_2[0]
					temp_final_df.loc[i, 'mother_minutes']=round(daughter_2[1])

					temp_final_df.loc[i, 'daughter_id']=daughter_1[0]
					temp_final_df.loc[i, 'daughter_minutes']=round(daughter_1[1])



	#A row index is added to the list of branches to drop only if diff_minutes is NA, which means that
	#one or both branches didn't divided.
	list_row_to_drop=[]
	for i in range(len(temp_final_df)):
		if(temp_final_df.diff_minutes[i]=='NA'):
			list_row_to_drop.append(i)

	#dropping the rows of the table
	temp_final_df=temp_final_df.drop(list_row_to_drop).reset_index(drop=True)

	#assigning some values to the new columns
	temp_final_df['date']=temp_date
	temp_final_df['strain']=temp_strain
	temp_final_df['timepoint']=temp_time_point
	temp_final_df['strain_timepoint']=temp_strain+"_"+temp_time_point
	temp_final_df['complete_id']=temp_strain+'_'+temp_time_point+'_'+temp_replicate+'_'+temp_date #updated no longer id_file, 
	#now is (strain)_(timepoint)_(replicate)_(date)
	#temp_final_df['small_id']=temp_strain+'_'+temp_replicate+'_'+temp_date #small id format: (strain)_(replicate)_(date) #modified for v2
	temp_final_df['id_file']=id_file
	

	#concatenating the dataframes
	if(not assigned):
		final_df=temp_final_df
		assigned=True
	else:
		final_df=pd.concat([final_df,temp_final_df], ignore_index=True)

	cont_file_done+=1
	print(str(cont_file_done)+'/'+str(total_files))

#saving the final dataframe
final_df.to_csv(args.output, sep=',', index=False)

