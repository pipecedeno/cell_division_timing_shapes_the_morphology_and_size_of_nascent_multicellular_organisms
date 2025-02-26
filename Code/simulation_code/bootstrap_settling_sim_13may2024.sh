#!/bin/bash

# Date: 13May2024
#
# Note: for the matlab code to be executed it needs to be in the matlab path. The codes that need to
# be executed to set this are:
# addpath('/path/to/your/script');
# savepath;


while getopts "f:s:n:o:e:c:m:g:p:k:b:h:i:u:y:" option
do
case "${option}"
in
f) strain1_dt_file=${OPTARG};;
s) strain2_dt_file=${OPTARG};;
n) strain_names=${OPTARG};;
o) output_dir=${OPTARG};; 
e) edge_degree_threshold=${OPTARG};;
c) pop_concentration=${OPTARG};;
m) sim_number=${OPTARG};;
g) selection_strength=${OPTARG};;
p) proportion_sampled_at_random=${OPTARG};;
k) growth_carrying_capacity=${OPTARG};;
b) bootstrap_total_size=${OPTARG};;
h) threshold_clust_gen_initial=${OPTARG};;
i) initial_growth_concentration=${OPTARG};;
u) number_bootstraps=${OPTARG};;
y) type_settling_selection=${OPTARG};;
*) echo "Invalid option or missing option argument for -$OPTARG" >&2; exit 1;;
esac
done



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

output_dir=$(realpath ${output_dir})/
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
initial_cluster_growth_15mar2024.py -f ${path_strain1_dt} -s none -n ${first_strain} -t ${temp_dir} -e ${edge_degree_threshold} \
	-c ${initial_growth_concentration} -k ${threshold_clust_gen_initial}
initial_cluster_growth_15mar2024.py -f ${path_strain2_dt} -s ${temp_dir}initial_node_inf_${first_strain}.csv -n ${second_strain} \
	-t ${temp_dir} -e ${edge_degree_threshold} -c ${initial_growth_concentration} -k ${threshold_clust_gen_initial}
#they had first_concentration and second_concentration before


initial_growth_end=`date +%s`
echo "" >> log.txt
echo Initial cluster growth took `expr $initial_growth_end - $initial_growth_start` seconds. >> log.txt
echo "" >> log.txt

final_start=`date +%s`

#growth phase of the simulation
# echo "growth_phase_clusters_v0_9feb2024.py -f ${path_strain1_dt} -s ${path_strain2_dt} -n ${strain_names} \
# 	-o ${output_dir} -t ${temp_dir} -e ${edge_degree_threshold} -c ${pop_concentration} -r 1 -i ${sim_number} \
# 	-k ${growth_carrying_capacity}"

growth_phase_clusters_v0_9feb2024.py -f ${path_strain1_dt} -s ${path_strain2_dt} -n ${strain_names} \
	-o ${output_dir} -t ${temp_dir} -e ${edge_degree_threshold} -c ${pop_concentration} -r 1 -i ${sim_number} \
	-k ${growth_carrying_capacity}

echo "${temp_dir}final_population_1.csv ${selection_strength} ${proportion_sampled_at_random} ${bootstrap_total_size} ${type_settling_selection} \
${output_dir}proportions_settling_selec_registry.csv ${first_strain} ${second_strain} ${sim_number} ${number_bootstraps}" > parameters.txt
	
matlab -nodisplay -r "settling_selection_bootstrap_v0_10may2024 ; exit"

final_end=`date +%s`
echo "" >> log.txt
echo Total execution time was `expr $final_end - $final_start` seconds. >> log.txt




