
## Time lapse dataset creation codes

In order to create the codes for the datasets the following codes need to be executed.

[Add explanation of what the codes do]

./merge_time_lapse_tables_v2.py -d ~/Desktop/timelapses_article/timelapse_csv_files/ -o ~/Desktop/timelapses_article/timelapse_joined_files_2025may5.csv

The file ```timelapse_joined_files_2025may5.csv``` still needs to be processed to add some columns so that it can be used in the R code synchrony_mother_daughter_31jan2024.csv

./tl_cell_sync_estimate_v3.py -d ~/Desktop/timelapses_article/timelapse_csv_files/ -o ~/Desktop/timelapses_article/cell_sync_mother_daughter_2025may5.csv
