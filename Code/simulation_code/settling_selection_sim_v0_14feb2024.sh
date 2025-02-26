#!/bin/bash

# Date: 18Feb2024
# This code will execute the main loop for the growth and settling selection simulations for several
#transfers similar to the evolution experiment
# Programs that are part of the main loop (executed by settling_selection_sim_v0_14feb2024.sh):
# 	-growth_phase_clusters_v0_9feb2024.py (growth phase of the experiment)
# 	-settling_selection_phase_v0_15feb2024.m (Obtaining volume of the clusters and settling selection)
# 	-write_settling_selection_output_v0_16feb2024.py (program to obtain the proportion of clusters that
# 		survived settling selection)
# Note: for the matlab code to be executed it needs to be in the matlab path. The codes that need to
# be executed to set this are:
# addpath('/path/to/your/script');
# savepath;

#modification history:
#Date: 3Apr2024
#the input selection_strength was modified to the number_of_cells that are going to be transfered to
#the next growth phase after the settling selection, so you will need to input what is the final diluation
#that you want from the random sample and the settling selection

#Date: 26Apr2024
# -Added the input carrying capacity for the matlab code, this change was made to make the matlab code automatically
#calculate how many cells should survive the settling selection
# -Added two options for the settling selection "random", "top"

while getopts "f:s:n:o:e:c:r:m:g:p:k:h:i:y:" option
do
case "${option}"
in
f) strain1_dt_file=${OPTARG};;
s) strain2_dt_file=${OPTARG};;
n) strain_names=${OPTARG};;
o) output_dir=${OPTARG};; 
e) edge_degree_threshold=${OPTARG};;
c) pop_concentration=${OPTARG};;
r) num_transfers=${OPTARG};;
m) sim_number=${OPTARG};;
g) selection_strength=${OPTARG};;
p) proportion_sampled_at_random=${OPTARG};;
k) carrying_capacity=${OPTARG};;
h) threshold_clust_gen_initial=${OPTARG};;
i) initial_growth_concentration=${OPTARG};;
y) type_settling_selection=${OPTARG};;
*) echo "Invalid option or missing option argument for -$OPTARG" >&2; exit 1;;
esac
done
#g) selection_strength=${OPTARG};;


# The bash script which executes this code in parallel already checks this so it is no longer 
# necessary to do this step
# #checking that the output directory doesn't exists to avoid overwriting it
# [ -d ${output_dir} ] && echo "output directory already exits, program won't continue" && exit

# It should already have it, but it is a good idea just to make sure that it does
#adding a "/" into the last character of the string if it doesn't has the "/" already
if [[ ! "${output_dir}" = */ ]]; then
	output_dir="${output_dir}/"
fi

path_strain1_dt=$(realpath ${strain1_dt_file})
path_strain2_dt=$(realpath ${strain2_dt_file})

#making output directories
mkdir -p ${output_dir}
mkdir -p ${output_dir}temp_dir

temp_dir=$(realpath ${output_dir}temp_dir)/

# Changing to the output directory as the matlab script receives a file called parameters.txt as 
# input, and because the name of the file is hardcoded we need to be in each directory to avoid any 
# miscommunications between the simulations running in parallel
cd ${temp_dir}

echo "Starting simulation" >> log.txt
echo "" >> log.txt
echo "Growing initial populations" >> log.txt

initial_growth_start=`date +%s`

#Extracting strain names and concentrations
first_strain=$(echo $strain_names | cut -d':' -f1)
second_strain=$(echo $strain_names | cut -d':' -f2)

first_concentration=$(echo $pop_concentration | cut -d':' -f1)
second_concentration=$(echo $pop_concentration | cut -d':' -f2)


#Simulating the initial populations
initial_cluster_growth_15mar2024.py -f ${path_strain1_dt} -s none -n ${first_strain} -t ${temp_dir} -e ${edge_degree_threshold} -c ${initial_growth_concentration} -k ${threshold_clust_gen_initial}
initial_cluster_growth_15mar2024.py -f ${path_strain2_dt} -s ${temp_dir}initial_node_inf_${first_strain}.csv -n ${second_strain} -t ${temp_dir} -e ${edge_degree_threshold} -c ${initial_growth_concentration} -k ${threshold_clust_gen_initial}
#they had first_concentration and second_concentration before


initial_growth_end=`date +%s`
echo "" >> log.txt
echo Initial cluster growth took `expr $initial_growth_end - $initial_growth_start` seconds. >> log.txt
echo "" >> log.txt

final_start=`date +%s`

for i in $(seq 1 $num_transfers); do
	echo "" >> log.txt
	echo "transfer ${i}" >> log.txt

	transfer_start=`date +%s`

	#growth phase of the simulation
	growth_phase_clusters_v0_9feb2024.py -f ${path_strain1_dt} -s ${path_strain2_dt} -n ${strain_names} \
		-o ${output_dir} -t ${temp_dir} -e ${edge_degree_threshold} -c ${pop_concentration} -r ${i} -i ${sim_number} \
		-k ${carrying_capacity}

	#settling simulation in matlab
	# saving the inputs for the matlab script in the parameters.txt file
	# echo "${carrying_capacity}"
	echo "${temp_dir}final_population_${i}.csv ${temp_dir}surviving_ids_${i}.txt ${selection_strength} ${temp_dir}volumes_${i}.txt ${proportion_sampled_at_random} ${carrying_capacity} ${type_settling_selection}" > parameters.txt
	
	matlab -nodisplay -r "settling_selection_phase_v0_15feb2024 ; exit"

	# If there is an error in the matlab script is won't stop the execution of the program, so if 
	# that happens we need to stop the loop to stop errors happening in any of the other scripts, so 
	# this will help with trouble shooting in case of any error
	if [ ! -f "${temp_dir}surviving_ids_${i}.txt" ]; then
		echo "File ${temp_dir}surviving_ids_${i}.txt not found. Stopping the loop." >> log.txt
		break # Exit the loop if file not found
	fi

	#save proportions and information of the settling selection
	write_settling_selection_output_v0_16feb2024.py -n ${temp_dir}node_inf_${i}.csv -o ${output_dir}proportions_settling_selec_registry.csv \
		-r ${i} -s ${temp_dir}surviving_ids_${i}.txt -i ${sim_number} -v ${temp_dir}volumes_${i}.txt -p ${output_dir}volumes_settling_selection.csv \
		-m ${strain_names}

	if [ -f "stop_extinction.txt" ]; then
		echo "Population went extinct, ending simulation." >> log.txt
		break # Exit the loop if file not found
	fi

	transfer_end=`date +%s`
	echo Transfer took `expr $transfer_end - $transfer_start` seconds. >> log.txt

done


final_end=`date +%s`
echo "" >> log.txt
echo Total execution time was `expr $final_end - $final_start` seconds. >> log.txt


#moving log file to the main folder and deleting the temp_dir
mv ${temp_dir}log.txt ${output_dir}
rm -r ${temp_dir}
