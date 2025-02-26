#!/bin/bash

while getopts "d:n:o:g:e:.pyt:" option
do
case "${option}"
in
d) doubling_time_file=${OPTARG};;
n) number_simulations=${OPTARG};;
o) output_dir=${OPTARG};; 
g) num_generations=${OPTARG};;
e) edge_degree_threshold=${OPTARG};;
t) num_cores=${OPTARG};;
*) echo "Invalid option or missing option argument for -$OPTARG" >&2; exit 1;;
esac
done


#checking that the output directory doesn't exists to avoid overwriting it
[ -d ${output_dir} ] && echo "output directory already exits, program won't continue" && exit

#adding a "/" into the last character of the string if it doesn't has the "/" already
if [[ ! "${output_dir}" = */ ]]; then
	output_dir="${output_dir}/"
fi


# Making the output directory
mkdir -p ${output_dir}
mkdir -p ${output_dir}temp_dir

start=`date +%s`

echo "Command used:" >> ${output_dir}log.txt
echo "par_sim_frag_all_clusters_exponential_11mar2024.sh -d ${doubling_time_file} -n ${number_simulations} -o ${output_dir} -g ${num_generations} -e ${edge_degree_threshold} -t ${num_cores}" >> ${output_dir}log.txt
echo "" >> ${output_dir}log.txt

echo "$(seq 1 ${number_simulations})" | parallel -P ${num_cores} sim_frag_clust_edge_degree_all_clusters_exponential_9mar2024.py -d ${doubling_time_file} \
-n {} -o ${output_dir}temp_dir/ -g ${num_generations} -e ${edge_degree_threshold}

find ${output_dir}temp_dir -name 'degree_distribution.csv' > ${output_dir}files_to_concat.txt
concat_files_output_18feb2024.py -f ${output_dir}files_to_concat.txt -o ${output_dir}degree_distribution.csv

find ${output_dir}temp_dir -name 'diff_doub_t.csv' > ${output_dir}files_to_concat.txt
concat_files_output_18feb2024.py -f ${output_dir}files_to_concat.txt -o ${output_dir}diff_doub_t.csv

find ${output_dir}temp_dir -name 'fragmentation_inf.csv' > ${output_dir}files_to_concat.txt
concat_files_output_18feb2024.py -f ${output_dir}files_to_concat.txt -o ${output_dir}fragmentation_inf.csv

find ${output_dir}temp_dir -name 'networks_diameter.csv' > ${output_dir}files_to_concat.txt
concat_files_output_18feb2024.py -f ${output_dir}files_to_concat.txt -o ${output_dir}networks_diameter.csv

#deleting the files_to_concat.txt file
rm ${output_dir}files_to_concat.txt

#deleting temp_dir
rm -rf ${output_dir}temp_dir

echo "" >> ${output_dir}log.txt
echo "Code completed!" >> ${output_dir}log.txt

end=`date +%s`
echo Total execution time was `expr $end - $start` seconds. >> ${output_dir}log.txt