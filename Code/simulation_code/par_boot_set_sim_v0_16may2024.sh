#!/bin/bash


# Date: 18Feb2024
# This code is used to execute the growth and settling selection simulations in parallel, the main loop
# of the program is in the code settling_selection_sim_v0_14feb2024.sh which executes each of the
# individual simulations.
# For notes about the program check each the programs executed and the notes in the lab notebook in
# OneNote.
# Programs that are part of the main loop (executed by settling_selection_sim_v0_14feb2024.sh):
# 	-growth_phase_clusters_v0_9feb2024.py (growth phase of the experiment)
# 	-settling_selection_phase_v0_15feb2024.m (Obtaining volume of the clusters and settling selection)
# 	-write_settling_selection_output_v0_16feb2024.py (program to obtain the proportion of clusters that
# 		survived settling selection)

# Modifications:
# Date: 26Apr2024
# Added two options for the settling selection "random", "top"


while getopts "f:s:n:o:e:c:m:t:g:p:k:b:h:i:u:y:" option
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
t) num_cores=${OPTARG};;
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

# echo ${output_dir}

#checking that the output directory doesn't exists to avoid overwriting it
[ -d ${output_dir} ] && echo "output directory already exits, program won't continue" && exit

# if [[ "$type_settling_selection" != "top" && "$type_settling_selection" != "random" ]]; then
#     echo "Error in parameters, -y only accepts the options 'top' or 'random'"
#     exit
# fi


#adding a "/" into the last character of the string if it doesn't has the "/" already
if [[ ! "${output_dir}" = */ ]]; then
	output_dir="${output_dir}/"
fi


# Making the output directory
mkdir -p ${output_dir}

# Obtaining the realpath of the directories as in the other codes the directory is going to be changed
# to acces the temporary folders
output_dir=$(realpath ${output_dir})/

path_strain1_dt=$(realpath ${strain1_dt_file})
path_strain2_dt=$(realpath ${strain2_dt_file})

# Test code to show what is going to be executed
# echo "$(seq 1 ${sim_number})" | parallel -P ${num_cores} echo "bootstrap_settling_sim_13may2024.sh -f ${path_strain1_dt} -s ${path_strain2_dt} \
# -n ${strain_names} -o ${output_dir}sim_{} -e ${edge_degree_threshold} -c ${pop_concentration} -m {} -g ${selection_strength} \
# -p ${proportion_sampled_at_random} -k ${growth_carrying_capacity} -b ${bootstrap_total_size} -h ${threshold_clust_gen_initial} \
# -i ${initial_growth_concentration} -u ${number_bootstraps}"

time_start=`date +%s`

echo "Command used:" >> ${output_dir}log_boot_set_sim.txt
echo "par_boot_set_sim_v0_16may2024.sh" >> ${output_dir}log_boot_set_sim.txt
echo "  -f ${strain1_dt_file}" >> ${output_dir}log_boot_set_sim.txt
echo "  -s ${strain2_dt_file}" >> ${output_dir}log_boot_set_sim.txt
echo "  -n ${strain_names}" >> ${output_dir}log_boot_set_sim.txt
echo "  -o ${output_dir}" >> ${output_dir}log_boot_set_sim.txt
echo "  -e ${edge_degree_threshold}" >> ${output_dir}log_boot_set_sim.txt
echo "  -c ${pop_concentration}" >> ${output_dir}log_boot_set_sim.txt
echo "  -r ${num_transfers}" >> ${output_dir}log_boot_set_sim.txt
echo "  -m ${sim_number}" >> ${output_dir}log_boot_set_sim.txt
echo "  -t ${num_cores}" >> ${output_dir}log_boot_set_sim.txt
echo "  -g ${selection_strength}" >> ${output_dir}log_boot_set_sim.txt
echo "  -p ${proportion_sampled_at_random}" >> ${output_dir}log_boot_set_sim.txt
echo "  -k ${carrying_capacity}" >> ${output_dir}log_boot_set_sim.txt
echo "  -h ${threshold_clust_gen_initial}" >> ${output_dir}log_boot_set_sim.txt
echo "  -i ${initial_growth_concentration}" >> ${output_dir}log_boot_set_sim.txt
echo "  -y ${type_settling_selection}" >> ${output_dir}log_boot_set_sim.txt
echo "" >> ${output_dir}log_boot_set_sim.txt

echo "$(seq 1 ${sim_number})" | parallel -P ${num_cores} bootstrap_settling_sim_13may2024.sh -f ${path_strain1_dt} -s ${path_strain2_dt} \
-n ${strain_names} -o ${output_dir}sim_{} -e ${edge_degree_threshold} -c ${pop_concentration} -m {} -g ${selection_strength} \
-p ${proportion_sampled_at_random} -k ${growth_carrying_capacity} -b ${bootstrap_total_size} -h ${threshold_clust_gen_initial} \
-i ${initial_growth_concentration} -u ${number_bootstraps} -y ${type_settling_selection}

echo "Simulations completed!"

echo "Concatenating result files..."

# this one is commented out because is not interesting anymore for this simulations
# find ${output_dir} -name 'proportions_growth_registry.csv' > files_to_concat.txt
# concat_files_output_18feb2024.py -f files_to_concat.txt -o ${output_dir}proportions_growth_registry_all_sim.csv

find ${output_dir} -name 'proportions_settling_selec_registry.csv' > files_to_concat.txt
concat_files_output_18feb2024.py -f files_to_concat.txt -o ${output_dir}proportions_settling_selec_registry_all_sim.csv

# find ${output_dir} -name 'time_registry.csv' > files_to_concat.txt
# concat_files_output_18feb2024.py -f files_to_concat.txt -o ${output_dir}time_registry_all_sim.csv

find ${output_dir} -name 'test_max_initial_pos.csv' > files_to_concat.txt
concat_files_output_18feb2024.py -f files_to_concat.txt -o ${output_dir}test_max_initial_pos_all_sim.csv

rm files_to_concat.txt

# to delete the intermediate files if needed
# rm -r ${output_dir}/sim_*

echo "Program completed"


time_end=`date +%s`
echo "" >> ${output_dir}log_boot_set_sim.txt
echo Total execution time was `expr $time_end - $time_start` seconds. >> ${output_dir}log_boot_set_sim.txt
