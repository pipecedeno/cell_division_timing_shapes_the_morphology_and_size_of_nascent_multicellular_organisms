% Thomas C. Day, 2022 edit
% This script generates N clusters of snowflake yeast, modeled as connected
% cells, where each cell is a prolate ellipsoid. The stop condition of the
% experiment can either be a number of generations or a total overlap stop
% condition. Cell aspect ratio can also be varied. Last, there are options
% to include force relaxation to cell positions and orientations.

%{
% Luis Felipe 2023
Date: 18Feb2024

This code is based on settling_selection_phase_v0_15feb2024.m, some parts of the code may
be redundant or unnecesary as not all of its functionalities are used.
The idea of this code is to make the networks from the growth simulations into clusters and calculate
their volume in order to be able to perform settling selection, but in here there is a bootstrap performed
to increase the size of the population.

% inputs:
Network file: File with the edges of the networks and when each edge was created to build the 
    cluster in order
Output file: It is going to save the id of the clusters that survived settling selection and their 
    volume. The volume may not be necessary as the volume of all clusters is saved in another file
    Each value is separated by a comma and there is no header in the file 
Selection strength: What percentage of the volume is going to survive the settling selection
Volume file: what is the name of the file in which the volumes of all clusters are going to be saved
    -cluster_id
    -diameter
    -initial_position
    -final_position
    -distance_travelled
    -final_time
    -speed
Percentage sampled at random: Percentage of the number of clusters which is going to be sampled at 
    random from the whole population
Carrying Capacity: Carrying capacity of the grown population, this is to know how many cells should 
    survive
Type of settling selection: random if all clusters start at a random position, top if all clusters
    start at the top of the tube
Proportions_output: File name of the number of cells and clusters that survived settling selection
    Same columns as the settling selection proportions registry file

Network file
Selection Strength
Percentage sampled at random
Carrying capacity
Type of settling selection
Proportions output
Strain_1
Strain_2
sim_number
number_bootstraps



Parts of the code
Load variables
Load table
Get unique values in the cluster_id table

Initialize vectors using the length of the cluster ids list
-list of ids
-list of volumes

Make the main loop to iterate through the vector of unique cluster ids
-subset table
-order it according to the transfer_added and time_grown columns
-Grow cluster
-Save volume and id

-Do settling selection
-Get ids and save them into the output file
    settling_results_(curr_transfer).csv

Note: in order for this code to be executed in the terminal you need to add the path of the script
to the matlab path, here are the commands needed to do that:
# addpath('/path/to/your/script');
# savepath;
# You may need to edit the file pathdef.m in /Applications/MATLAB_R2023b.app/toolbox/local


Date: 30Jan2026
-Corrected how the diameter of the clusters is calculated, now the gyration radius of the clusters is calculated and 
then converted in the radius of a sphere to calculate the diameter of the clusters assuming that they have spherical shapes


%}

%% Inputs:
% folder = pwd;
% cd(folder);


%Loading input values

% % Open the file
fileID = fopen('parameters.txt','r');

% Read the line as a string
line = fgetl(fileID);

% Close the file
fclose(fileID);

% Split the line into string and float parts
parts = strsplit(line, ' ');
network_file = parts{1};
selectionStrength = str2double(parts{2});
percentage_sampled_at_random=str2double(parts{3});
carrying_capacity=str2num(parts{4});
type_settling_selection=parts{5};
proportions_filename=parts{6};
strain_1_name=parts{7};
strain_2_name=parts{8};
sim_number=str2num(parts{9});
number_bootstraps=str2num(parts{10});

num_surviving_cells=carrying_capacity*percentage_sampled_at_random*selectionStrength;
cells_to_sample=carrying_capacity*percentage_sampled_at_random;


% Load table of network edges
clusters_edges=readtable(network_file);

unique_cluster_ids=unique(clusters_edges.cluster_id);
num_cluster_ids=length(unique_cluster_ids);

% fprintf(['before', num2str(length(unique_cluster_ids)), '\n'])

%% This sampling is not necessary to do
% Sample 10% of the clusters at random 
% sampleSize = ceil(length(unique_cluster_ids) * percentage_sampled_at_random);
% sampleIndices = randperm(length(unique_cluster_ids), sampleSize);
% unique_cluster_ids = unique_cluster_ids(sampleIndices);
% num_cluster_ids=length(unique_cluster_ids);

% fprintf(['after', num2str(length(unique_cluster_ids)), '\n'])


% N               = 50; % number of clusters to generate
% numGens         = 5; % number of generations of cell division in a group
diam            = 4.58; % smallest diameter of the cell
err_diam        = 0; % variation in cell size, taken from coefficient of variation data
AR              = 1.2;%:.1:2.8; % aspect ratio value %# DEFAULT was 1.2
err_AR          = 0; % standard deviation in aspect ratio from Shane's measurements
pole_theta      = deg2rad(10); % buds nearest the pole will be chosen from between 0 and 10 degrees in polar angle
THETA           = deg2rad(45); % polar angle average from SEM data
thetaVariance   = deg2rad(0); % variation in polar angle from SEM data
distance_thresh = 1.1672; % minimum distance (um) separating bud scars
new_bud_prob    = .8; % probability that the first cell will bud near the pole
check_overlap   = 0; % do we check the overlap?
overlap_thresh  = 4.0e1; % threshold of total overlaps %#2.5e2

% Force input parameters:
include_forces  = 0; % do we include forces for rearrangements? 1 = Y, 0 = N
NEIGHBOR_THRESH = 3*diam; % radius at which you consider cells "neighbors"
STERIC_MAG      = 2; % relative strength of steric interactions(from cells colliding)
CHITIN_MAG      = 1; % relative strength of the bonds holding cells together
BOND_TORQ_MAG   = 1; % relative torsional strength of chitin bonds
mobility_pos    = .1;
mobility_rot    = .1;
dt              = .1; % timestep for mechanical relaxation
T               = 25; % number of timesteps to allow mechanical relaxation

% Do we visualize figures?
figure_viz = 0; % 1 = Y, 0 = N



% Creating input variables dictionary
global input_variables;
input_variables=containers.Map();
% input_variables('clust_size_output')=clust_size_output;
input_variables('aspect_ratio')=AR;
input_variables('neighbor_thresh')=NEIGHBOR_THRESH;
input_variables('THETA')=THETA;
input_variables('diam')=diam;
input_variables('err_diam')=err_diam;
input_variables('err_AR')=err_AR;
input_variables('pole_theta')=pole_theta;
input_variables('thetaVariance')=thetaVariance;
input_variables('distance_thresh')=distance_thresh;
input_variables('new_bud_prob')=new_bud_prob;
input_variables('check_overlap')=check_overlap;
input_variables('overlap_thresh')=overlap_thresh;
input_variables('NEIGHBOR_THRESH')=NEIGHBOR_THRESH;

input_variables('include_forces')=include_forces;
input_variables('NEIGHBOR_THRESH')=NEIGHBOR_THRESH;
input_variables('STERIC_MAG')=STERIC_MAG;
input_variables('CHITIN_MAG')=CHITIN_MAG;
input_variables('BOND_TORQ_MAG')=BOND_TORQ_MAG;
input_variables('mobility_pos')=mobility_pos;
input_variables('mobility_rot')=mobility_rot;
input_variables('dt')=dt;
input_variables('figure_viz')=figure_viz;
input_variables('T')=T;
input_variables('strain_1_name')=strain_1_name;
input_variables('strain_2_name')=strain_2_name;
input_variables('sim_number')=sim_number;
input_variables('proportions_filename')=proportions_filename;
input_variables('num_surviving_cells')=num_surviving_cells;
input_variables('cells_to_sample')=cells_to_sample;
input_variables('type_settling_selection')=type_settling_selection;



% Initialization:
fprintf('Welcome to the cell simulator v2.0\n');
fprintf('Written by Thomas C. Day, 2020\n');

%% MAIN CODE: -------------------------------------------------------------
% -------------------------------------------------------------------------


% creating a vector of each cluster strain
list_input_strains = cell(1, length(unique_cluster_ids));

for i = 1:length(unique_cluster_ids)
    temp_clust_id = unique_cluster_ids(i);
    clust_edge = clusters_edges(clusters_edges.cluster_id == temp_clust_id, :);
    list_input_strains{i} = clust_edge.strain{1};
end


% Calculate the volume of all the clusters (the idea is to perform this calculation only once)
list_diameters=zeros(1, length(unique_cluster_ids));
list_volumes=zeros(1, length(unique_cluster_ids));
cluster_sizes=zeros(1, length(unique_cluster_ids));

for i = 1:length(unique_cluster_ids)

    temp_clust_id=unique_cluster_ids(i);

    % Subsetting information with the cluster id
    clust_edge=clusters_edges(clusters_edges.cluster_id==temp_clust_id, :);

    % Sorting in both rows to avoid any errors (sorting only by time_grown was not enough as
    %that value gets reset every transfer)
    sorted_clust_edges = sortrows(clust_edge, {'transfer_added', 'time_grown'});
    % sorted_clust_edges=sortrows(clust_edge, 'time_grown');

    node1_vals      = sorted_clust_edges.node1;
    node2_vals      = sorted_clust_edges.node2;


    % aspRat    = AR(a);
    aspRat = input_variables('aspect_ratio');
    theta     = input_variables('THETA');
    

    % Grow cluster of the snowflake
    [cell_list] = ELYES_SIM_fixed_daughter(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, NEIGHBOR_THRESH, check_overlap, overlap_thresh, node1_vals, node2_vals);

    %Calculating the gyration radius
    temp_gyr_rad=gyration_radius(cell_list);

    % Transforming the gyration radius to a spherical radius, and obtaining the diameter
    temp_diameter=2*(sqrt(5*(temp_gyr_rad^2)/3));

    list_diameters(i)=temp_diameter;
    temp_V=(1/6) * pi * (temp_diameter^3);
    list_volumes(i)=temp_V;
    cluster_sizes(i)=length(cell_list);
    
    if(length(cell_list)~=length(node1_vals)+1)
        fprintf('Different amount of cells\n');
    end

end


% fprintf('Size of list_volumes: %d\n', length(list_volumes));
% fprintf('Size of list_diameters: %d\n', length(list_diameters));
% fprintf('Size of cluster_sizes: %d\n', length(cluster_sizes));


% Performing the settling selections in the bootstrap
for num_boot=1:number_bootstraps
    [survivors, list_sampled_strains, sampled_volumes, sampled_cluster_sizes]=perform_boostrap_settling_selection(num_boot, clusters_edges, list_input_strains, unique_cluster_ids, list_diameters, list_volumes, cluster_sizes);
    fprintf('%d/%d\n', num_boot, number_bootstraps);
end


function [survivors, list_sampled_strains, sampled_volumes, sampled_cluster_sizes]=perform_boostrap_settling_selection(num_boot, clusters_edges, list_input_strains, unique_cluster_ids, list_diameters, list_volumes, cluster_sizes)

    global input_variables;
    strain_1_name=input_variables('strain_1_name');
    strain_2_name=input_variables('strain_2_name');
    sim_number=input_variables('sim_number');
    proportions_filename=input_variables('proportions_filename');
    num_surviving_cells=input_variables('num_surviving_cells');
    cells_to_sample=input_variables('cells_to_sample');
    type_settling_selection=input_variables('type_settling_selection');

    % diam=input_variables('diam');
    % err_diam=input_variables('err_diam');
    % err_AR=input_variables('err_AR');
    % pole_theta=input_variables('pole_theta');
    % thetaVariance=input_variables('thetaVariance');
    % distance_thresh=input_variables('distance_thresh');
    % new_bud_prob=input_variables('new_bud_prob');
    % check_overlap=input_variables('check_overlap');
    % overlap_thresh=input_variables('overlap_thresh');
    % NEIGHBOR_THRESH=input_variables('NEIGHBOR_THRESH');

    sampled_cluster_ids=[];
    sampled_diameters=[];
    sampled_volumes=[];
    sampled_cluster_sizes=[];
    list_sampled_strains={};

    sampled_num_cells=0;

    cont_clusters=1;

    while sampled_num_cells<cells_to_sample

        % fprintf([num2str(unique_cluster_ids(i)), '\n']);
        % fprintf([num2str(i), '/', num2str(length(unique_cluster_ids)), '\n']);

        randomIndex = randi(length(unique_cluster_ids));
        temp_clust_id=unique_cluster_ids(randomIndex);
        % fprintf('%d\n', randomIndex);

        sampled_cluster_ids(cont_clusters)=temp_clust_id;
        sampled_diameters(cont_clusters)=list_diameters(randomIndex);
        sampled_volumes(cont_clusters)=list_volumes(randomIndex);
        sampled_cluster_sizes(cont_clusters)=cluster_sizes(randomIndex);
        list_sampled_strains{cont_clusters}=list_input_strains(randomIndex);

        % Update number of cells sampled
        sampled_num_cells=sum(sampled_cluster_sizes);
        % fprintf('%d\n', length(cell_list));
        % fprintf('%d/%.2f\n', sampled_num_cells, cells_to_sample);

        cont_clusters=cont_clusters+1;
    end

    log_file = fopen('log.txt', 'a');
    fprintf(log_file, '\nBootstrap Information\nNumber of cells bootstrapped: %d/%.2f\n', sampled_num_cells, cells_to_sample);
    fprintf(log_file, 'Number of clusters selected: %d\n', length(sampled_cluster_sizes));
    fclose(log_file);


    % run settling selection function
    if contains(type_settling_selection, "top")
        [survivors, initial_pos, final_pos, distance_travelled, final_time, speedsize]=settling_selection_simulated_top(sampled_diameters, sampled_cluster_sizes, num_surviving_cells);
    else
        [survivors, initial_pos, final_pos, distance_travelled, final_time, speedsize]=settling_selection_simulated_random(sampled_diameters, sampled_cluster_sizes, num_surviving_cells);
    end


    % Save number of clusters and cells that survived

    %Create header if file doesn't exist

    if isfile(proportions_filename)
        proportions_output = fopen(proportions_filename, 'a');
    else
        proportions_output = fopen(proportions_filename, 'w');
        fprintf(proportions_output, 'sim_number,bootstrap_num,strain1,strain2,clusters_pop1_b,clusters_pop2_b,cells_pop1_b,cells_pop2_b,total_clusters_b,total_cells_b,clusters_pop1_a,clusters_pop2_a,cells_pop1_a,cells_pop2_a,total_clusters_a,total_cells_a\n');
    end

    % Convert cell arrays in list_sampled_strains to strings
    list_sampled_strains_str = cellfun(@char, list_sampled_strains, 'UniformOutput', false);

    strain_1_inf=containers.Map();
    strain_1_inf('num_clusters_before')=sum(strcmp(list_sampled_strains_str, strain_1_name));
    strain_1_inf('num_clusters_after')=sum(strcmp(list_sampled_strains_str, strain_1_name) & survivors);
    strain_1_inf('num_cells_before')=sum(sampled_cluster_sizes(strcmp(list_sampled_strains_str, strain_1_name)));
    strain_1_inf('num_cells_after')=sum(sampled_cluster_sizes(strcmp(list_sampled_strains_str, strain_1_name) & survivors));


    strain_2_inf=containers.Map();
    strain_2_inf('num_clusters_before')=sum(strcmp(list_sampled_strains_str, strain_2_name));
    strain_2_inf('num_clusters_after')=sum(strcmp(list_sampled_strains_str, strain_2_name) & survivors);
    strain_2_inf('num_cells_before')=sum(sampled_cluster_sizes(strcmp(list_sampled_strains_str, strain_2_name)));
    strain_2_inf('num_cells_after')=sum(sampled_cluster_sizes(strcmp(list_sampled_strains_str, strain_2_name) & survivors));


    total_cells_before=strain_2_inf('num_cells_before')+strain_2_inf('num_cells_before');
    total_clusters_before=strain_2_inf('num_clusters_before')+strain_2_inf('num_clusters_before');
    total_cells_after=strain_2_inf('num_cells_after')+strain_2_inf('num_cells_after');
    total_clusters_after=strain_2_inf('num_clusters_after')+strain_2_inf('num_clusters_after');

    fprintf(proportions_output, '%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n', sim_number, num_boot, strain_1_name, strain_2_name, ...
        strain_1_inf('num_clusters_before'), strain_2_inf('num_clusters_before'), strain_1_inf('num_cells_before'), ...
        strain_2_inf('num_cells_before'), total_clusters_before, total_cells_before, strain_1_inf('num_clusters_after'), ...
        strain_2_inf('num_clusters_after'), strain_1_inf('num_cells_after'), strain_2_inf('num_cells_after'), ...
        total_clusters_after, total_cells_after);

    fclose(proportions_output);


    % Save file of maximum distance of the bottom of the cells that survived per strain

    if isfile('../test_max_initial_pos.csv')
        test_initial_pos = fopen('../test_max_initial_pos.csv', 'a');
    else
        test_initial_pos = fopen('../test_max_initial_pos.csv', 'w');
        fprintf(test_initial_pos, 'sim_number,bootstrap_num,strain1,strain2,max_position_1,max_position_2\n');
    end

    max_position_strain_1=max(initial_pos(strcmp(list_sampled_strains_str, strain_1_name) & survivors));
    max_position_strain_2=max(initial_pos(strcmp(list_sampled_strains_str, strain_2_name) & survivors));

    fprintf(test_initial_pos, '%d,%d,%s,%s,%f,%f\n', sim_number, num_boot, strain_1_name, strain_2_name, ...
        max_position_strain_1, max_position_strain_2);

    fclose(test_initial_pos);

end


%% OLD function, this function is ineficient as it is growing the clusters every time they are sampled, which shouldn't be
%% necesary as the size of the clusters is always going to be somewhat similar
% function [survivors, list_sampled_strains, list_volumes, cluster_sizes]=perform_boostrap_settling_selection(num_boot, clusters_edges, list_input_strains, unique_cluster_ids)

%     global input_variables;
%     strain_1_name=input_variables('strain_1_name');
%     strain_2_name=input_variables('strain_2_name');
%     sim_number=input_variables('sim_number');
%     proportions_filename=input_variables('proportions_filename');
%     num_surviving_cells=input_variables('num_surviving_cells');
%     cells_to_sample=input_variables('cells_to_sample');
%     diam=input_variables('diam');
%     err_diam=input_variables('err_diam');
%     err_AR=input_variables('err_AR');
%     pole_theta=input_variables('pole_theta');
%     thetaVariance=input_variables('thetaVariance');
%     distance_thresh=input_variables('distance_thresh');
%     new_bud_prob=input_variables('new_bud_prob');
%     check_overlap=input_variables('check_overlap');
%     overlap_thresh=input_variables('overlap_thresh');
%     NEIGHBOR_THRESH=input_variables('NEIGHBOR_THRESH');
%     type_settling_selection=input_variables('type_settling_selection');

%     list_cluster_ids=[];
%     list_diameters=[];
%     list_volumes=[];
%     cluster_sizes=[];
%     list_sampled_strains={};

%     sampled_num_cells=0;

%     cont_clusters=1;

%     while sampled_num_cells<cells_to_sample

%         % fprintf([num2str(unique_cluster_ids(i)), '\n']);
%         % fprintf([num2str(i), '/', num2str(length(unique_cluster_ids)), '\n']);

%         randomIndex = randi(numel(unique_cluster_ids));
%         temp_clust_id=unique_cluster_ids(randomIndex);

%         % Subsetting information with the cluster id
%         clust_edge=clusters_edges(clusters_edges.cluster_id==temp_clust_id, :);

%         % Sorting in both rows to avoid any errors (sorting only by time_grown was not enough as
%         %that value gets reset every transfer)
%         sorted_clust_edges = sortrows(clust_edge, {'transfer_added', 'time_grown'});
%         % sorted_clust_edges=sortrows(clust_edge, 'time_grown');

%         node1_vals      = sorted_clust_edges.node1;
%         node2_vals      = sorted_clust_edges.node2;


%         % aspRat    = AR(a);
%         aspRat = input_variables('aspect_ratio');
%         theta     = input_variables('THETA');
        

%         % Grow cluster of the snowflake
%         [cell_list] = ELYES_SIM_fixed_daughter(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, NEIGHBOR_THRESH, check_overlap, overlap_thresh, node1_vals, node2_vals);

%         %Calculating cluster diameter using the gyration radius
%         temp_diameter=gyration_diameter(cell_list);

%         list_cluster_ids(cont_clusters)=temp_clust_id;
%         list_diameters(cont_clusters)=temp_diameter;
%         temp_V=(1/6) * pi * (temp_diameter^3);
%         list_volumes(cont_clusters)=temp_V;
%         cluster_sizes(cont_clusters)=length(cell_list);
%         list_sampled_strains{end+1}=list_input_strains(randomIndex);
        
%         if(length(cell_list)~=length(node1_vals)+1)
%             fprintf('Different amount of cells\n');
%         end

%         cont_clusters=cont_clusters+1;

%         % Update number of cells sampled
%         sampled_num_cells=sum(cluster_sizes);
%         % fprintf('%d\n', length(cell_list));
%         % fprintf('%d/%.2f\n', sampled_num_cells, cells_to_sample);
%     end

%     log_file = fopen('log.txt', 'a');
%     fprintf(log_file, 'Bootstrap Information\nNumber of cells bootstrapped: %d/%.2f\n', sampled_num_cells, cells_to_sample);
%     fprintf(log_file, 'Number of clusters selected: %d\n', length(cluster_sizes));
%     fclose(log_file);


%     % run settling selection function
%     if contains(type_settling_selection, "top")
%         [survivors, initial_pos, final_pos, distance_travelled, final_time, speedsize]=settling_selection_simulated_top(list_diameters, cluster_sizes, num_surviving_cells);
%     else
%         [survivors, initial_pos, final_pos, distance_travelled, final_time, speedsize]=settling_selection_simulated_random(list_diameters, cluster_sizes, num_surviving_cells);
%     end


%     % Save number of clusters and cells that survived

%     %Create header if file doesn't exist

%     if isfile(proportions_filename)
%         proportions_output = fopen(proportions_filename, 'a');
%     else
%         proportions_output = fopen(proportions_filename, 'w');
%         fprintf(proportions_output, 'sim_number,bootstrap_num,strain1,strain2,clusters_pop1_b,clusters_pop2_b,cells_pop1_b,cells_pop2_b,total_clusters_b,total_cells_b,clusters_pop1_a,clusters_pop2_a,cells_pop1_a,cells_pop2_a,total_clusters_a,total_cells_a\n');
%     end

%     % Convert cell arrays in list_sampled_strains to strings
%     list_sampled_strains_str = cellfun(@char, list_sampled_strains, 'UniformOutput', false);

%     strain_1_inf=containers.Map();
%     strain_1_inf('num_clusters_before')=sum(strcmp(list_sampled_strains_str, strain_1_name));
%     strain_1_inf('num_clusters_after')=sum(strcmp(list_sampled_strains_str, strain_1_name) & survivors);
%     strain_1_inf('num_cells_before')=sum(cluster_sizes(strcmp(list_sampled_strains_str, strain_1_name)));
%     strain_1_inf('num_cells_after')=sum(cluster_sizes(strcmp(list_sampled_strains_str, strain_1_name) & survivors));


%     strain_2_inf=containers.Map();
%     strain_2_inf('num_clusters_before')=sum(strcmp(list_sampled_strains_str, strain_2_name));
%     strain_2_inf('num_clusters_after')=sum(strcmp(list_sampled_strains_str, strain_2_name) & survivors);
%     strain_2_inf('num_cells_before')=sum(cluster_sizes(strcmp(list_sampled_strains_str, strain_2_name)));
%     strain_2_inf('num_cells_after')=sum(cluster_sizes(strcmp(list_sampled_strains_str, strain_2_name) & survivors));


%     total_cells_before=strain_2_inf('num_cells_before')+strain_2_inf('num_cells_before');
%     total_clusters_before=strain_2_inf('num_clusters_before')+strain_2_inf('num_clusters_before');
%     total_cells_after=strain_2_inf('num_cells_after')+strain_2_inf('num_cells_after');
%     total_clusters_after=strain_2_inf('num_clusters_after')+strain_2_inf('num_clusters_after');

%     fprintf(proportions_output, '%d,%d,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n', sim_number, num_boot, strain_1_name, strain_2_name, ...
%         strain_1_inf('num_clusters_before'), strain_2_inf('num_clusters_before'), strain_1_inf('num_cells_before'), ...
%         strain_2_inf('num_cells_before'), total_clusters_before, total_cells_before, strain_1_inf('num_clusters_after'), ...
%         strain_2_inf('num_clusters_after'), strain_1_inf('num_cells_after'), strain_2_inf('num_cells_after'), ...
%         total_clusters_after, total_cells_after);

%     fclose(proportions_output);

% end



%% SETTLING SELECTION FUNCTION %%
%{
The idea of the simulation is to sample a random starting position for where they are going to start
in the "eppendorf tube" for the settling selection. Only the volume that is defined by the 
selectionStrength is going to survive to the next simulation.
How the settling slection is done is that the clusters start at a random position and they need
to go below a distance threshold from the bottom of the eppendorf tube (defined as 0) for them to 
survive, the loop is going to stop until more than the desired volume has passed the threshold.
The clusters position is going to be updated in small time steps
%}


function [survivors, initial_pos, final_pos, distance_travelled, final_time, speedsize] = settling_selection_simulated_random(diameters, cluster_sizes, num_surviving_cells)
    % Define physical constants and eppendorf tube size
    massDensWater = 997;       % kg/m^3
    dynViscoWater = .000891;   % Pa.s
    massDensParticle = 1112.6; % kg/m^3
    gravity = 9.81;            % m/s^2
    eppSize = 2.5E-2;          % 2mL eppendorf size in meters

    % Calculate the speed for each cluster
    speedsize = (1/18) * ((massDensParticle - massDensWater) / dynViscoWater) * gravity * (diameters.^2);

    % Sample random starting positions for each cluster
    positions = rand(size(diameters)) * eppSize;
    initial_pos=positions;

    % Define a delta t value (modifiable for tuning)
    delta_t = 1e-2; % seconds
    cont_steps=0; %count how many times the clusters were moved 

    below_threshold=true;
    % Main loop for settling selection
    while below_threshold
        % Move clusters towards the bottom of the tube
        positions = positions - speedsize * delta_t;
        positions(positions < 0) = 0; % Ensure no position is less than 0

        %counting how many cells have reached the bottom of the tube
        surviving_cells=sum(cluster_sizes(positions==0));

        % Stop if the surviving volume is greater than the selection strength
        if surviving_cells >= num_surviving_cells
            below_threshold=false;
        end

        cont_steps=cont_steps+1;
    end

    %final positions 
    final_pos=positions;

    % Return the vector of true/false for clusters that survived
    % survivors = positions < survival_threshold;
    survivors=positions==0;

    % fprintf(['Steps settling selection: ', num2str(cont_steps), '\n']);
    % fprintf(['Initial surviving proportion: ', num2str(initialSurvivingVolume), '\n']);
    % fprintf(['Final surviving proportion: ', num2str(survivingVolume), '\n']);
    % fprintf(['Number of surviving clusters: ', num2str(sum(survivors)), '\n']);
    % fprintf(['Final time: ', num2str(cont_steps*delta_t), 'sec  ', num2str(cont_steps*delta_t/60), 'min\n']);

    %this information is written in the log file to adjust some parameters if necessary
    log_file = fopen('log.txt', 'a');
    fprintf(log_file, 'Random settling selection being used\n');
    fprintf(log_file, 'Steps settling selection: %d\n', cont_steps);
    fprintf(log_file, 'Final surviving number of cells: %d\n', surviving_cells);
    fprintf(log_file, 'Number of surviving clusters: %d\n', sum(survivors));
    fprintf(log_file, 'Final time: %f sec, %f min\n', cont_steps*delta_t, cont_steps*delta_t/60);
    
    fclose(log_file);

    %Total distance travelled
    %note: it is initial_pos-final_pos because the distance is how far
    %away they are from the bottom, so this way the distance is going to be
    %positive
    distance_travelled=initial_pos-final_pos;

    %final time 
    final_time=cont_steps*delta_t;
end


function [survivors, initial_pos, final_pos, distance_travelled, final_time, speedsize] = settling_selection_simulated_top(diameters, cluster_sizes, num_surviving_cells)
    % Define physical constants and eppendorf tube size
    massDensWater = 997;       % kg/m^3
    dynViscoWater = .000891;   % Pa.s
    massDensParticle = 1112.6; % kg/m^3
    gravity = 9.81;            % m/s^2
    eppSize = 2.5E-2;          % 2mL eppendorf size in meters

    % Calculate the speed for each cluster
    speedsize = (1/18) * ((massDensParticle - massDensWater) / dynViscoWater) * gravity * (diameters.^2);

    % Sample random starting positions for each cluster
    % positions = rand(size(diameters)) * eppSize;
    positions=ones(1, length(diameters)) * eppSize;
    initial_pos=positions;

    % Define a delta t value (modifiable for tuning)
    delta_t = 1e-2; % seconds
    cont_steps=0; %count how many times the clusters were moved 

    below_threshold=true;
    % Main loop for settling selection
    while below_threshold
        % Move clusters towards the bottom of the tube
        positions = positions - speedsize * delta_t;
        positions(positions < 0) = 0; % Ensure no position is less than 0

        %counting how many cells have reached the bottom of the tube
        surviving_cells=sum(cluster_sizes(positions==0));

        % Stop if the surviving volume is greater than the selection strength
        if surviving_cells >= num_surviving_cells
            below_threshold=false;
        end

        cont_steps=cont_steps+1;
    end

    %final positions 
    final_pos=positions;

    % Return the vector of true/false for clusters that survived
    % survivors = positions < survival_threshold;
    survivors=positions==0;

    % fprintf(['Steps settling selection: ', num2str(cont_steps), '\n']);
    % fprintf(['Initial surviving proportion: ', num2str(initialSurvivingVolume), '\n']);
    % fprintf(['Final surviving proportion: ', num2str(survivingVolume), '\n']);
    % fprintf(['Number of surviving clusters: ', num2str(sum(survivors)), '\n']);
    % fprintf(['Final time: ', num2str(cont_steps*delta_t), 'sec  ', num2str(cont_steps*delta_t/60), 'min\n']);

    %this information is written in the log file to adjust some parameters if necessary
    log_file = fopen('log.txt', 'a');
    fprintf(log_file, 'Top settling selection being used\n');
    fprintf(log_file, 'Steps settling selection: %d\n', cont_steps);
    fprintf(log_file, 'Final surviving number of cells: %f\n', surviving_cells);
    fprintf(log_file, 'Number of surviving clusters: %d\n', sum(survivors));
    fprintf(log_file, 'Final time: %f sec, %f min\n', cont_steps*delta_t, cont_steps*delta_t/60);
    
    fclose(log_file);

    %Total distance travelled
    %note: it is initial_pos-final_pos because the distance is how far
    %away they are from the bottom, so this way the distance is going to be
    %positive
    distance_travelled=initial_pos-final_pos;

    %final time 
    final_time=cont_steps*delta_t;
end


% Function to obtain the gyration radius of a cluster, it calculates the radius, converts it to meters
function [diameter]=gyration_radius(cell_list)
    % Number of points (cells)
    num_points = length(cell_list);
    
    % Extract the center of each cell into a Nx3 matrix
    centers = zeros(num_points, 3);
    for i = 1:num_points
        centers(i, :) = cell_list(i).Center.';
    end
    
    % Find cluster center
    cluster_center = mean(centers, 1);
    
    % Distance from every point to the center
    distances = sqrt(sum((centers - cluster_center).^2, 2));
    
    % Gyration radius
    Rg = sqrt(sum(distances.^2) / num_points)*10^-6;
    diameter=Rg;
end


%% FUNCTIONS


function [cell_list] = ELYES_SIM_fixed_daughter(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, neighbor_thresh, check_overlap, overlap_thresh, node1_vals, node2_vals)

    global input_variables;
    % clust_size_output=input_variables('clust_size_output');
    T_include_forces=input_variables('T');



    include_forces=input_variables('include_forces');
    NEIGHBOR_THRESH=input_variables('NEIGHBOR_THRESH');
    STERIC_MAG=input_variables('STERIC_MAG');
    CHITIN_MAG=input_variables('CHITIN_MAG');
    BOND_TORQ_MAG=input_variables('BOND_TORQ_MAG');
    mobility_pos=input_variables('mobility_pos');
    mobility_rot=input_variables('mobility_rot');
    dt=input_variables('dt');
    figure_viz=input_variables('figure_viz');

    cell_list       = [];  % this will eventually be the tabulated list
    
    % Create root cell:
    rootCell.Center       = [0; 0; 0];
    rootCell.Radii        = 1/2 * [aspRat * diam; diam; diam];
    rootCell.Rmatrix      = eye(3);
    rootCell.Generation   = 0;
    rootCell.IDnumber     = 0; 
    rootCell.networkID    = node1_vals(1); %# this line was modified as not all clusters are going to start with their first node being 1
    rootCell.Daughters    = [];
    rootCell.Mother       = [];
    rootCell.BudXYZ       = [];
    rootCell.Overlaps     = [];
    rootCell.Neighborhood = [];
    cell_list             = [cell_list, rootCell];
    counting_ix           = 0; % this will count how many cells are in the cluster

    cont_failed_to_add_edge=0;
    cont_missing_mother_cell=0;


    % Iterate over generations:
    % for g = 1:numGens
    % Iterate over the edges of the loaded file
    num_rows_edges=length(node1_vals);
    for cont_rows = 1:num_rows_edges


        %# I need to find the position of the mother cell in the list cell_list by checking if cell_list(n).IDnumber
        %# is equal to node1_vals(cont_rows), if it isn't found then we need to skip the row as that cell was not able to
        %# be added in a past iteration
        %# After the for loop if the cell id of the mother was found, then n is going to be the position
        mother_cell_id=node1_vals(cont_rows);
        daughter_cell_id=node2_vals(cont_rows);

        nPossible   = length(cell_list);
        cell_id_in_list=0;
        for n = 1:nPossible
            %# Here I had before IDnumber
            if(cell_list(n).networkID == mother_cell_id)
                cell_id_in_list=1;
                break;
            end
        end

        % fprintf(['cell_id_in_list: ', num2str(cell_id_in_list),'\n']);

        %%% This part of the code no longer should be executed
        %# If the cellid of the mother wasn't found in the list, then the code will skip to the next edge to add 
        if cell_id_in_list==0
            % fprintf(['cell not found in list\n'])
            cont_missing_mother_cell=cont_missing_mother_cell+1;
            continue;
        end


        % Finding the position of the daughter even if it is close to an existing bud scar
        variedTheta = theta + thetaVariance * randn(1);
        [newBud, newBudRel, newAxis] = getDaughterPos_fixed_daughter(cell_list(n), pole_theta, variedTheta, distance_thresh, new_bud_prob); % finds the new bud xyz position and relative position to the old cell
        % this if can be deleted as it always returns newBud as the position of the new bud
        if newBud ~= 0 % if the budding chance was successful
            counting_ix = counting_ix + 1; % the number of cells goes up by one

            % record new bud in mother cell information array:
            cell_list(n).BudXYZ = [cell_list(n).BudXYZ, newBudRel];
            cell_list(n).Daughters = [cell_list(n).Daughters, counting_ix];

            % record new cell information for the cell list:
            newCell.IDnumber     = counting_ix;
            newCell.Generation   = 0; %# how to count generations? this are cells generation it was g before but I set it as 0 for the moment
            newCell.Mother       = cell_list(n).IDnumber;
            newCell.networkID    = daughter_cell_id;
            newCell.Daughters    = [];
            newCell.BudXYZ       = [];
            newCell.Overlaps     = [];
            newCell.Neighborhood = [];
            newCell.Radii        = 1/2 * [aspRat * diam; diam; diam];

            % Find new cell center:
            [Sm, Rm, Tm]   = GET_SURFACE_MATRICES(cell_list(n), .1);
            A              = newCell.Radii(1);
            p              = A * newAxis + newBudRel;
            q              = Tm * Rm * [p; 1];
            newCell.Center = q(1:3);

            % Choose a cell orientation based upon surface normal axis:
            avec   = newCell.Center - newBud(1:3);
            avec   = avec./norm(avec);
            bvec_x = rand;
            bvec_y = rand;
            bvec_z = - (avec(1)*bvec_x + avec(2)*bvec_y)/avec(3);
            bvec   = [bvec_x; bvec_y; bvec_z];
            bvec   = bvec./norm(bvec);
            cvec   = cross(avec, bvec);
            newCell.Rmatrix = [avec, bvec, cvec];

            % add new cell to cell list:
            cell_list = [cell_list, newCell];
        end
        
        % fprintf(['num_cells: ',num2str(length(cell_list)),'\n']);
        
    end

end


function [newPos, relativePos, axis] = getDaughterPos_fixed_daughter(CELL, POLE_TH, TH, DISTANCE_THRESH, NEW_BUD_PROB_THRESH)

    % Surface definition for possible mother cell:
    [Sm, Rm, Tm] = GET_SURFACE_MATRICES(CELL, .1);

    % Determine where on cell body to bud the next scar:
    if isempty(CELL.Daughters) % reproduce near the pole with prob 0.7
        frac = rand;
        if frac < NEW_BUD_PROB_THRESH
            % bud at pole
            th = POLE_TH * rand;
            ph = 2 * pi * rand;
            % ix_too_close = [];
            r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
            x = r * cos(th);
            y = r * sin(th) * cos(ph);
            z = r * sin(th) * sin(ph);
        else
            % bud on the polar angle
            ph = 2 * pi * rand;
            th = TH;
            % ix_too_close = [];
            r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
            x = r * cos(th);
            y = r * sin(th) * cos(ph);
            z = r * sin(th) * sin(ph);
        end

    else % reproduce based upon the polar angle theta, sometimes back bud (if 3 or more bud scars)
        
        nDaughters = size(CELL.BudXYZ,2);
        if nDaughters < 4
            ph = 2 * pi * rand;
            th = TH;
        else
            frac = rand;
            if frac < 1 % bud at the distal pole
                ph = 2 * pi * rand;
                th = TH;
            else % bud at proximal pole
                ph = 2 * pi * rand;
                th = deg2rad(180) - TH;
            end
        end

        % Check if new location is too close to existing scars:
        r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
        x = r * cos(th);
        y = r * sin(th) * cos(ph);
        z = r * sin(th) * sin(ph);
        t = [x; y; z];
        s = CELL.BudXYZ - t;
        d = sqrt(s(1,:).^2 + s(2,:).^2 + s(3,:).^2);
        ix_too_close = find(d < DISTANCE_THRESH);
    end

    % if isempty(ix_too_close) % new bud is successful
    %     relativePos = [x; y; z];
    %     newPos      = Tm * Rm * [relativePos; 1];
    %     normal      = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
    %     axis        = normal./norm(normal);
    % else % this bud is too close, loses its chance
    %     relativePos = 0;
    %     newPos      = 0;
    %     axis        = 0;
    % end

    relativePos = [x; y; z];
    newPos      = Tm * Rm * [relativePos; 1];
    normal      = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
    axis        = normal./norm(normal);

end


function [flag] = CHECK_OVERLAPS2(cell_list, overlap_thresh)

    % Obtain distances between all particles:
    [D, ~, ~, ~] = get_particle_distances(cell_list);
    
    % Obtain combined equatorial radii of two particles:
    radii = [cell_list.Radii];
    eq_radii = radii(2,:);
    R = eq_radii + eq_radii';
    
    %{
    % Obtain the energy associated with each interaction:
    d = D - R;
    w = d(d < 0);
    u = w.^2;
    U = sum(u,'all'); % sum up all interaction energies
    
    flag = U > overlap_thresh; % if U is too big, flag is thrown
    %}
    
    % New method:
    overlapAmt = D./R;
    overlap = overlapAmt < 1;
    overlap = overlap - eye(size(overlap));
    ix = find(overlap);
    w = 1 - overlapAmt(ix);
    u = w.^(5/2);
    U = sum(u);
    flag = U > overlap_thresh;
    
end

function [flag] = CHECK_OVERLAPS(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude, overlap_thresh)

    [~, ~, ~, ~, overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude);
    Overlaps = sum(overlaps, 2);
    if ~isempty(Overlaps)
        flag = Overlaps(2) > overlap_thresh;
    else
        flag = 0;
    end
    

end

function parsave(filename, cell_list)
    save(filename, 'cell_list');
end

function [S, R, T] = GET_SURFACE_MATRICES(cell_of_interest, scaling)
    % A surface matrix is a 4x4 matrix that will transform a point on the
    % unit sphere to a point on an arbitrary ellipsoid. There are three
    % sub-matrices that compose a surface matrix: a scaling matrix (S), a
    % rotation matrix (R) and a translation matrix (T).

    % Scaling matrix
    S = [cell_of_interest.Radii(1) + scaling, 0, 0, 0;...
            0, cell_of_interest.Radii(2) + scaling, 0, 0; ...
            0, 0, cell_of_interest.Radii(3) + scaling, 0; ...
            0,0,0,1];
    %  cell wall is usually ~75 nm thick, this accounts for that
    
    % rotation matrix
    R = cell_of_interest.Rmatrix;
    R = [R; 0,0,0]; R = horzcat(R,[0;0;0;1]);
    
    % translation matrix
    T = eye(4); 
    T(1,end) = cell_of_interest.Center(1);
    T(2,end) = cell_of_interest.Center(2);
    T(3,end) = cell_of_interest.Center(3);

end


function [cell_list_in, cell_list_out, Overlaps] = INCLUDE_FORCES(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude, mobility_pos, mobility_rot, dt, T, figure_viz)
% This function takes as input a grown snowflake and calculates the forces
% and torques on each cell. From these forces, it the calculates any
% rearrangements.

cell_list_in = cell_list;

for t = 1:T
    fprintf(['Time = ',num2str(t),'\n']);
    [Forces, Torques, flog1, flog2, Overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude);
    [cell_list] = UPDATE_POSITIONS(cell_list, Forces, Torques, mobility_pos, mobility_rot, dt);
    
    if figure_viz == 1
        if mod(t,1) == 0
            figure(1); clf;
            hold on; box on; set(gca,'linewidth',2);
            for n = 1:length(cell_list)
                c_o_i = cell_list(n);
                [x,y,z] = VISUALIZE_ELLIPSOID(c_o_i, 30);
                surf(x,y,z,'facealpha',1,'edgecolor','none');
                %{
%                 Fnet = 3*Forces(:,n) + c_o_i.Center;
%                 Tnet = 3*Torques(:,n) + c_o_i.Center;
%                 plot3([c_o_i.Center(1), Fnet(1)], [c_o_i.Center(2),Fnet(2)], [c_o_i.Center(3), Fnet(3)],'r-','linewidth',3);
%                 plot3([c_o_i.Center(1), Tnet(1)], [c_o_i.Center(2),Tnet(2)], [c_o_i.Center(3), Tnet(3)],'b-','linewidth',3);
                forces = flog1{n};
                ix = find(vecnorm(forces) ~= 0);
                forces(:,ix) = forces(:,ix)./vecnorm(forces(:,ix));
                pts = flog2{n};
                plot3(pts(1,1), pts(2,1), pts(3,1), 'rx','markersize',12,'linewidth',2); % plot birth scar location
                plot3([pts(1,1), pts(1,1)+forces(1,1)], [pts(2,1), pts(2,1)+forces(2,1)], [pts(3,1), pts(3,1)+forces(3,1)],'r-','linewidth',2); % plot force from birth scar
                nDaughters = length(cell_list(n).Daughters);
                if nDaughters == 0
                    nDaughters = 1;
                end
                for jj = 1:nDaughters
                    plot3(pts(1,jj+1), pts(2,jj+1), pts(3,jj+1),'b.','markersize',10,'linewidth',2);
                    plot3([pts(1,jj+1), pts(1,jj+1)+forces(1,jj+1)], [pts(2,jj+1), pts(2,jj+1)+forces(2,jj+1)], [pts(3,jj+1), pts(3,jj+1)+forces(3,jj+1)],'b-','linewidth',2); % plot force from daughters
                end
                plot3(pts(1,nDaughters+2:end), pts(2,nDaughters+2:end), pts(3,nDaughters+2:end),'go','markersize',14,'linewidth',2);
                plot3([pts(1,nDaughters+2:end), pts(1,nDaughters+2:end) + forces(1,nDaughters+2:end)],[pts(2,nDaughters+2:end), pts(2,nDaughters+2:end) + forces(2,nDaughters+2:end)], [pts(3,nDaughters+2:end), pts(3,nDaughters+2:end) + forces(3,nDaughters+2:end)],'g-','linewidth',2);
                %}
            end
            view(3); axis equal;
            lightangle(15,15);
            lighting gouraud;
            material dull;
            centers = [cell_list.Center];
            xlim([min(centers(1,:)) - 5, max(centers(1,:)) + 5]);
            ylim([min(centers(2,:)) - 5, max(centers(2,:)) + 5]);
            zlim([min(centers(3,:)) - 5, max(centers(3,:)) + 5]);
            title([num2str(t)]);
            drawnow;
            % print(['test_sims_t=',num2str(t,'%03.f')],'-dpng','-r500');
            % uncomment this line to print the file
        end
    end
end
    
cell_list_out = cell_list;

end


function [cell_list_out] = UPDATE_POSITIONS(cell_list, Forces, Torques, mobility_pos, mobility_rot, dt)
    % Update cell positions and orientations based on the forces and
    % torques input, and mobility factors.

    % Overdamped motion
    cell_list_out = cell_list;
    dx = mobility_pos * Forces * dt; % magnitude and direction of displacement
    zero_magnitude = zeros(size(vecnorm(Torques)));
    da = - mobility_rot * vecnorm(Torques) * dt; % magnitude of rotational change
    Tvec = Torques./vecnorm(Torques); % axis of rotation
    
    % Update cell positions:
    for n = 1:length(cell_list)
        cell_list_out(n).Center = cell_list(n).Center + dx(:,n);
    end
    
    % Update cell orientations:
    for n = 1:length(cell_list)
        if da(n) ~= 0 % only rotate if there is a torque
            Rrot = rotation_around_arb_axis(Tvec(:,n), da(n));
            cell_list_out(n).Rmatrix = Rrot * cell_list(n).Rmatrix;
        end
    end

end


function [Forces, Torques, forces_log1, forces_log2, final_overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude)
    % Calculate the forces and torques acting on all cells in the group.
    % This function can get quite heavy.

    % Obtain distances between all particles:
    [D, ~, ~, ~] = get_particle_distances(cell_list);
    
    % Find all pairs of cells within the neighborhood threshold:
    N = D < neighbor_thresh; % finds all neighbor cells
    N = N - eye(size(N)); % we don't care about counting a cell as its own neighbor
    
    % Obtain forces/torques for each cell:
    Forces = zeros(3, length(cell_list));
    Torques = zeros(3, length(cell_list));
    final_overlaps = [];
    for n = 1:length(cell_list)
        neighbors = find(N(:,n) == 1); % only need to check the neighboring cells for interactions
        
        % Steric interactions:
        [F_pts, O_pts, Overlaps, Directions] = get_overlaps(cell_list, n, neighbors);
        final_overlaps = [final_overlaps, Overlaps];
        
        % Chitin interactions:
        [m_pts, d_pts] = CHITIN_INTERACTIONS(cell_list, n);
        
        % Calculate forces:
        % Force from chitin bond with mother cell
        if isempty(m_pts)
            F_mother = [0;0;0];
        else
            F_mother = chitin_magnitude * m_pts;
        end
        % Forces from chitin bonds with daughters
        if isempty(d_pts)
            F_daughter = [0;0;0];
        else
            F_daughter = chitin_magnitude * d_pts;
        end
        % Forces from steric interactions
        if ~isempty(Overlaps)
            F_steric = steric_magnitude * Overlaps(2,:) .* Directions;
        else
            F_steric = [0;0;0];
        end
        % All forces
        forces = [F_mother, F_daughter, F_steric];
        Forces(:,n) = sum(forces, 2);
        
        % Log forces and force pts for visualization:
        forces_log1{n} = forces;
        [S,R,T] = GET_SURFACE_MATRICES(cell_list(n),0);
        M = T*R*S;
        mother_log = M*[-1;0;0;1];
        if ~isempty(cell_list(n).Daughters)
            daughter_log = T*R*[cell_list(n).BudXYZ; ones(1,length(cell_list(n).Daughters))];
            if isempty(F_pts)
                F_pts = mother_log(1:3);            
            end
            forces_log2{n} = [mother_log(1:3), daughter_log(1:3,:), F_pts];
        else
            if isempty(F_pts)
                F_pts = mother_log(1:3);
            end
            forces_log2{n} = [mother_log(1:3), mother_log(1:3), F_pts];
        end
        
        % Calculate torques:
        % First calculate the r-vector:
        [S, R, T]  = GET_SURFACE_MATRICES(cell_list(n), 0); 
        M = T*R*S;
        if isempty(cell_list(n).Mother)
            r_mother = [0;0;0;1];
        else
            r_mother   = M * [-1; 0; 0; 1] - M * [0;0;0;1];
        end
        if ~isempty(cell_list(n).Daughters)
            r_daughter = T * R * [cell_list(n).BudXYZ; ones(1,size(cell_list(n).BudXYZ,2))] - M*[0;0;0;1];
        else
            r_daughter = [0;0;0;1];
        end
        r_steric   = [F_pts; ones(1,size(F_pts,2))] - M*[0;0;0;1];
        if isempty(r_steric)
            r_steric = [0;0;0;1];
        end
        r_angle = -r_mother(1:3);
        if ~isempty(cell_list(n).Mother)
            % n 
            cell_mother = cell_list(n).Mother + 1;
            ix = find(cell_list(cell_mother).Daughters == n - 1);
            budxyz = cell_list(cell_mother).BudXYZ(:,ix);
            [s, r, t] = GET_SURFACE_MATRICES(cell_list(cell_mother),0);
            surf_norm = 2 * [budxyz(1)/s(1,1)^2; budxyz(2)/s(2,2)^2; budxyz(3)/s(3,3)^2];
            preferred_axis = r * [surf_norm; 1];
            current_axis = R * S * [1;0;0;1];
            preferred_axis = preferred_axis(1:3)/norm(preferred_axis(1:3));
            current_axis = current_axis(1:3)/norm(current_axis(1:3));
            v = bond_torque_magnitude * (preferred_axis - current_axis);
        else
            v = [0;0;0];
        end
        forces = [forces, v];
        rdisp = [r_mother(1:3,:), r_daughter(1:3,:), r_steric(1:3,:), r_angle];
        torques = zeros(size(forces));
        for i = 1:size(forces,2)
            torques(:,i) = cross(forces(:,i), rdisp(:,i));
        end
        Torques(:,n) = sum(torques, 2);
    end

end

function [Force_pts, Overlap_pts, Overlaps, Directions] = get_overlaps(cell_list, ix_o_i, nbors)
    
    % Obtain overlapping volume and center point of overlap:
    cell1       = cell_list(ix_o_i);
    Force_pts   = [];
    Overlap_pts = [];
    Overlaps    = [];
    Directions  = [];
    for j = 1:length(nbors)
        cell2 = cell_list(nbors(j));
        
        % Check if the two cells overlap at all: --------------------------  
        % First numerically approximate the surface of each cell
        [x1,y1,z1] = VISUALIZE_ELLIPSOID(cell1, 20); % numerical surface
        [x2,y2,z2] = VISUALIZE_ELLIPSOID(cell2, 20);
        x1 = x1(:); y1 = y1(:); z1 = z1(:);
        x2 = x2(:); y2 = y2(:); z2 = z2(:);
        r1 = [x1,y1,z1,ones(size(x1))]';
        r2 = [x2,y2,z2,ones(size(x2))]';
        
        % Second find the analytic surface of each ellipsoid
        [S1,R1,T1] = GET_SURFACE_MATRICES(cell1, 0); % find analytic surface matrix
        [S2,R2,T2] = GET_SURFACE_MATRICES(cell2, 0);
        M1 = T1 * R1 * S1; % surface matrices
        M2 = T2 * R2 * S2;
        
        % Third transform the approximate surface of E2 into the analytic
        % coords of E1
        r2_E1 = M1 \ r2;        % all the points defining E2, in E1 coords
        D = vecnorm(r2_E1(1:3,:)); % distance from origin in E1 space
        ix = find(D < 1); % any points less than 1 are within the surface of E1
        
        % Main calculation:
        if length(ix) < 10 % they do not overlap
            overlap = []; % add nothing to the overlaps
            COM     = [];
            Fpt_1   = [];
            direction = [];
        elseif (10 <= length(ix)) && (length(ix) < 100) % they only overlap by a little bit
            r1_E2   = M2 \ r1;
            D_help  = vecnorm(r1_E2(1:3,:));
            ix_E2   = find(D_help < 1);
            inters  = [r1(1:3,ix_E2), r2(1:3,ix)]';
            %# Sometimes convhull gives and error saying there are not enough unique points to 
            % calculate the volume, and by testing it seems that the error is caused because the
            % inters variable repeats the same value for each of it columns, as if the intesection
            % of the ellipses was just a single point, this code is going to test if that is the 
            % case, and assing the value as if there was no overlap in that case, and it they are
            % different then it is going to do the same computation as it had in the original function

            hasMultipleUniqueValues = false(1, size(inters, 2));

            % Check each column for multiple unique values
            for i = 1:size(inters, 2)
                uniqueValues = unique(inters(:, i));
                hasMultipleUniqueValues(i) = length(uniqueValues) > 3;
            end

            if all(hasMultipleUniqueValues)
                [~,vol] = convhull(inters(:,1),inters(:,2),inters(:,3));
                overlap = [nbors(j); vol];
                COM     = mean(inters)';
                Fpt_1   = get_force_point(COM, cell1);% Find force point on the surface of the cell of interest
                direc   = M1 * [0;0;0;1];
                direction = direc(1:3) - Fpt_1;
                direction = direction/norm(direction);
            else
                overlap = []; % add nothing to the overlaps
                COM     = [];
                Fpt_1   = [];
                direction = [];
            end
        else % they overlap by a significant amount
            v1       = [x1,y1,z1]; % vertices defining cell1
            v2       = [x2,y2,z2]; % vertices defining cell2
            r1_E2    = M2 \ r1;
            D_help   = vecnorm(r1_E2(1:3,:));
            ix_E2    = find(D_help < 1);
            inters   = [r1(1:3,ix_E2), r2(1:3,ix)]';
            COM      = mean(inters)'; % a point common to both cells
            [inters] = INTERSECTION(v1, v2, COM');
            [~,vol]  = convhull(inters(:,1), inters(:,2), inters(:,3)); % returns the volume of overlap
            overlap  = [nbors(j); vol]; % add to overlaps array
            COM      = mean(inters)'; % center of the overlap volume
            Fpt_1    = get_force_point(COM, cell_list(ix_o_i));
            direc    = M1 * [0;0;0;1];
            direction = direc(1:3) - Fpt_1;
            direction = direction/norm(direction);
            % Find force point on the surface of the cell of interest:
        end
        % Store info to cell array:
        Force_pts = [Force_pts, Fpt_1];
        Overlap_pts = [Overlap_pts, COM];
        Overlaps = [Overlaps, overlap];
        Directions = [Directions, direction];
%         cell_list(ix_o_i).Force_pts   = [cell_list(ix_o_i).Force_pts, Fpt_1];
%         cell_list(ix_o_i).Overlap_pts = [cell_list(ix_o_i).Overlap_pts, COM];
%         cell_list(ix_o_i).Overlaps    = [cell_list(ix_o_i).Overlaps, overlap];
    end
    
%     % Trim to unique overlap readings:
%     if ~isempty(cell_list(ix_o_i).Overlaps)
%         [~,ix_unique,~] = unique(cell_list(ix_o_i).Overlaps(1,:));
%         cell_list(ix_o_i).Overlaps    = cell_list(ix_o_i).Overlaps(:,ix_unique);
%         cell_list(ix_o_i).Overlap_pts = cell_list(ix_o_i).Overlap_pts(:,ix_unique);
%         cell_list(ix_o_i).Force_pts   = cell_list(ix_o_i).Force_pts(:,ix_unique);
%     end

end

function [m_pts, d_pts] = CHITIN_INTERACTIONS(cell_list, n)
% Obtain all the chitin interactions acting on cell n:

    % Mother interaction:
    mother_cell = cell_list(n).Mother;
    if ~isempty(mother_cell) % if this cell has a mother
        cell_1 = cell_list(mother_cell + 1);
        cell_2 = cell_list(n);
        [Sm, Rm, Tm] = GET_SURFACE_MATRICES(cell_1, 0);
        [Sd, Rd, Td] = GET_SURFACE_MATRICES(cell_2, 0);
        ix = find(cell_1.Daughters == n - 1);
        pt_1 = Tm * Rm * [cell_1.BudXYZ(:,ix); 1]; % point on mother cell surface
        pt_2 = Td * Rd * Sd * [-1; 0; 0; 1]; % point on daughter cell surface
        m_pts = pt_1(1:3) - pt_2(1:3); % displacement vector
    else
        m_pts = [0;0;0];
    end
    
    % Daughter interactions:
    daughter_cells = cell_list(n).Daughters;
    d_pts = [];
    if ~isempty(daughter_cells) % if this cell has any daughters
        for d = 1:length(daughter_cells)
            cell_1 = cell_list(n);
            cell_2 = cell_list(daughter_cells(d)+1);
            [Sm, Rm, Tm] = GET_SURFACE_MATRICES(cell_1, 0);
            [Sd, Rd, Td] = GET_SURFACE_MATRICES(cell_2, 0);
            pt_1 = Tm * Rm * [cell_1.BudXYZ(:,d); 1];
            pt_2 = Td * Rd * Sd * [-1; 0; 0; 1];
            d_pts = [d_pts, pt_2(1:3) - pt_1(1:3)];
        end
    else
        d_pts = [0;0;0];
    end
end

function [Fpt] = get_force_point(COM, cell_o_i)
    % Given a point that lies within the ellipsoid, extrapolate out to the
    % surface of the ellipsoid.
    [S,R,T] = GET_SURFACE_MATRICES(cell_o_i,0);
    M       = T * R * S;
    COM_E   = M \ [COM; 1];
    Fpt_E   = [COM_E(1:3)/norm(COM_E(1:3)); 1];
    Fpt_g   = M * Fpt_E;
    Fpt     = Fpt_g(1:3);

end

function [x,y,z] = VISUALIZE_ELLIPSOID(cell_of_interest, resolution)
    % Make plottable information for cells

    % Extract info from cell_of_interest:
    radii   = cell_of_interest.Radii;
    centers = cell_of_interest.Center;
    R       = cell_of_interest.Rmatrix;
    
    % Generate data for "unrotated" ellipsoid
    [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3), resolution);
    
    % Rotate data with orientation matrix R and center T
    a = kron(R(:,1), xc);
    b = kron(R(:,2), yc);
    c = kron(R(:,3), zc);
    data = a+b+c; n = size(data,2);
    
    % Store for output:
    x = data(1:n,:) + centers(1); 
    y = data(n+1:2*n,:) + centers(2); 
    z = data(2*n+1:end,:) + centers(3);

end

function [D, X, Y, Z] = get_particle_distances(cell_list)
    % Determine distances between each pair of particles:
    r_in = [cell_list.Center]';
    X = r_in(:,1) - r_in(:,1)';
    Y = r_in(:,2) - r_in(:,2)';
    Z = r_in(:,3) - r_in(:,3)';
    D = sqrt(X.^2 + Y.^2 + Z.^2);
end

function [vertices_intersection] = INTERSECTION(vertices1, vertices2, center)
    % This is a reaally fun function. If you are reading this, please enjoy
    % how clever this is. If you want to nerd out about it, or just learn
    % more, you can shoot me an email!

    % Let vertices1 denote the vertices of a polyhedron, and let vertices2
    % denote the vertices of a second polyhedron. The center is a point which
    % is common to both polyhedra. Then, computing the intersection of the two 
    % polyhedra, vertices_intersect = vertices1 intersect vertices2, is the
    % same problem as computing the union of the duals. Whaaaaattt????
    % v_dual_intersect = v_dual1 union v_dual2
    % The dual of polyhedron 1 is calculated in the following way:
    %   1. all vertices are shifted to place the center at the origin.
    %   2. the constraint equations such that A*x <= b are calculated. (this is
    %       just an alternative method of writing down the polyhedron. we can
    %       write this as a1*x + a2*y + a3*z <= d
    %   3. the values of the constraint equations are swapped with the vertex
    %       coordinates. in other words, the dual vertices are (a1, a2, a3) and
    %       the planes defining the polyhedron are x1*x + x2*y + x3*z <= d.
    %   4. Take the union of all the dual vertices D = [A1, A2].
    %   5. Now find the convex hull of the polyhedron D, which returns vertices
    %       verticesD.
    %   6. Compute the dual polyhedron of verticesD.
    %   7. Shift by center back to the original place
    % The dual polyhedron is our intersection volume.
    
    % Determine size/orientation of vertices:
    swap_later = 0;
    if size(center,1) == 3
        swap_later = 1;
        vertices1 = vertices1';
        vertices2 = vertices2';
        center    = center';
    end
    
    % 1. Shift all vertices to the origin given by center:
    vertices1a = vertices1 - center;
    vertices2a = vertices2 - center;
    
    % 2. Calculate the constraint equations:
    [A1,b1] = vert2lcon(vertices1a);
    [A2,b2] = vert2lcon(vertices2a);
    
    % Set b to very small nonzero constant if it is zero:
    b1_zero = b1==0;
    b2_zero = b2==0;
    for m = 1:length(b1)
        if b1_zero(m)
            b1(m) = 1e-3*rand;
        end
    end
    for m = 1:length(b2)
        if b2_zero(m)
            b2(m) = 1e-3*rand;
        end
    end
    
    % Normalize the constraint equations:
    A1 = -A1./b1; % must normalize the constraint equations
    A2 = -A2./b2; % normalize
    
    % 3. The dual vertices are now given by A1 and A2:
    vdual1 = A1;
    vdual2 = A2;
    
    % 4. Take the union of all dual vertices:
    DualA = [vdual1; vdual2];
    
    % 5. Find the convex hull of the dual:
    [k, ~] = convhull(DualA);
    v = unique(k); % only care about the unique vertices
    vertices_dual = DualA(v,:);
    
    % 6. Compute the dual of the polyhedron formed by vertices_dual:
    [Ax, bx] = vert2lcon(vertices_dual);
    Ax = -Ax./bx; % must renormalize the constraint equations again
    
    % 7. Shift by the center back to original place:
    vertices_intersection = Ax + center;
    
    % ReSwap the orientation of the lists if necessary:
    if swap_later == 1
        vertices_intersection = vertices_intersection';
    end
    
    % Plot for visuals:
    %{
    figure; hold on; box on; grid on; set(gca,'linewidth',2);
    [k1,~] = convhull(vertices1);
    [k1d,~] = convhull(vdual1);
    [k2,~] = convhull(vertices2);
    [k2d,~] = convhull(vdual2);
    [k3,~] = convhull(vertices_intersection);
    [k3d,~] = convhull(vertices_dual);
    plot3(center(1), center(2), center(3), '.','markersize',20,'color',[.8,.5,0]);
    trisurf(k1, vertices1(:,1), vertices1(:,2), vertices1(:,3),'facecolor',[.7,.7,.7],'facealpha',.3,'edgecolor','k','linewidth',3);
    trisurf(k2, vertices2(:,1), vertices2(:,2), vertices2(:,3),'facecolor',[.7,.1,.7],'facealpha',.2,'edgecolor',[.4,0,.4],'linewidth',3);
    trisurf(k3, vertices_intersection(:,1), vertices_intersection(:,2), vertices_intersection(:,3), 'facecolor',[.7,0,0],'facealpha',.5,'edgecolor',[.4,0,0],'linewidth',4);
    view(3);
    axis equal;
    title('Real space');
    
    figure; hold on; box on; grid on; set(gca,'linewidth',2);
    trisurf(k1d, vdual1(:,1), vdual1(:,2), vdual1(:,3), 'facecolor',[.7,.7,.7],'facealpha',.5,'edgecolor','b','linewidth',3);
    trisurf(k2d, vdual2(:,1), vdual2(:,2), vdual2(:,3), 'facecolor',[.7,.1,.7],'facealpha',.5,'edgecolor','k','linewidth',2);
    trisurf(k3d, vertices_dual(:,1), vertices_dual(:,2), vertices_dual(:,3),'facecolor', [.7,0,0],'facealpha',.5,'edgecolor',[.4,0,0],'linewidth',4);
    view(3);
    axis equal;
    title('Dual space');
    %}

end

function [R] = rotation_around_arb_axis(axis, angle)
    % This function returns a rotation matrix R for a rotation around the
    % arbitrary axis by some angle. See the Wikipedia on Rotation matrices for
    % more information about how to derive this.
    
    % INPUTS:
    % axis: must be a normalized 3-vector
    % angle: the amount to rotate
    
    % OUTPUTS:
    % R: rotation matrix
    
    Rxx = cos(angle) + axis(1)^2 * (1-cos(angle));
    Rxy = axis(1)*axis(2)*(1-cos(angle)) - axis(3)*sin(angle);
    Rxz = axis(1)*axis(3)*(1-cos(angle)) + axis(2)*sin(angle);
    Ryx = axis(1)*axis(2)*(1-cos(angle)) + axis(3)*sin(angle);
    Ryy = cos(angle) + axis(2)^2 * (1-cos(angle));
    Ryz = axis(2)*axis(3)*(1-cos(angle)) - axis(1)*sin(angle);
    Rzx = axis(1)*axis(3)*(1-cos(angle)) - axis(2)*sin(angle);
    Rzy = axis(2)*axis(3)*(1-cos(angle)) + axis(1)*sin(angle);
    Rzz = cos(angle) + axis(3)^2 * (1-cos(angle));
    
    R = [Rxx, Rxy, Rxz;
         Ryx, Ryy, Ryz;
         Rzx, Rzy, Rzz];

end
%% End Thomas C. Day functions

function [A,b,Aeq,beq]=vert2lcon(V,tol)
    %An extension of Michael Kleder's vert2con function, used for finding the 
    %linear constraints defining a polyhedron in R^n given its vertices. This 
    %wrapper extends the capabilities of vert2con to also handle cases where the 
    %polyhedron is not solid in R^n, i.e., where the polyhedron is defined by 
    %both equality and inequality constraints.
    % 
    %SYNTAX:
    %
    %  [A,b,Aeq,beq]=vert2lcon(V,TOL)
    %
    %The rows of the N x n matrix V are a series of N vertices of a polyhedron
    %in R^n. TOL is a rank-estimation tolerance (Default = 1e-10).
    %
    %Any point x inside the polyhedron will/must satisfy
    %  
    %   A*x  <= b
    %   Aeq*x = beq
    %
    %up to machine precision issues.
    %
    %
    %EXAMPLE: 
    %
    %Consider V=eye(3) corresponding to the 3D region defined 
    %by x+y+z=1, x>=0, y>=0, z>=0.
    %
    % 
    %   >>[A,b,Aeq,beq]=vert2lcon(eye(3))
    %
    %
    %     A =
    % 
    %         0.4082   -0.8165    0.4082
    %         0.4082    0.4082   -0.8165
    %        -0.8165    0.4082    0.4082
    % 
    % 
    %     b =
    % 
    %         0.4082
    %         0.4082
    %         0.4082
    % 
    % 
    %     Aeq =
    % 
    %         0.5774    0.5774    0.5774
    % 
    % 
    %     beq =
    % 
    %         0.5774
      %%initial stuff
      
        if nargin<2, tol=1e-10; end
        [M,N]=size(V);
        
        if M==1
          A=[];b=[];
          Aeq=eye(N); beq=V(:);
          return
        end
        
        
        
        
        p=V(1,:).';
        X=bsxfun(@minus,V.',p);
        
        
        %In the following, we need Q to be full column rank 
        %and we prefer E compact.
        
        if M>N  %X is wide
            
         [Q, R, E] = qr(X,0);  %economy-QR ensures that E is compact.
                               %Q automatically full column rank since X wide
                               
        else%X is tall, hence non-solid polytope
            
         [Q, R, P]=qr(X);  %non-economy-QR so that Q is full-column rank.
         
         [~,E]=max(P);  %No way to get E compact. This is the alternative. 
            clear P
        end
        
        
       diagr = abs(diag(R));
        
       if nnz(diagr)    
           
            %Rank estimation
            r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
        
        
            iE=1:length(E);
            iE(E)=iE;
           
           
            Rsub=R(1:r,iE).';
            if r>1
              [A,b]=vert2con(Rsub,tol);
             
            elseif r==1
                
               A=[1;-1];
               b=[max(Rsub);-min(Rsub)];
            end
            A=A*Q(:,1:r).';
            b=bsxfun(@plus,b,A*p);
            
            if r<N
             Aeq=Q(:,r+1:end).';      
             beq=Aeq*p;
            else
               Aeq=[];
               beq=[];
            end
       else %Rank=0. All points are identical
          
           A=[]; b=[];
           Aeq=eye(N);
           beq=p;
           
       end
       
       
    %            ibeq=abs(beq);
    %             ibeq(~beq)=1;
    %            
    %            Aeq=bsxfun(@rdivide,Aeq,ibeq);
    %            beq=beq./ibeq;
end
           
           
function [A,b] = vert2con(V,tol)
    % VERT2CON - convert a set of points to the set of inequality constraints
    %            which most tightly contain the points; i.e., create
    %            constraints to bound the convex hull of the given points
    %
    % [A,b] = vert2con(V)
    %
    % V = a set of points, each ROW of which is one point
    % A,b = a set of constraints such that A*x <= b defines
    %       the region of space enclosing the convex hull of
    %       the given points
    %
    % For n dimensions:
    % V = p x n matrix (p vertices, n dimensions)
    % A = m x n matrix (m constraints, n dimensions)
    % b = m x 1 vector (m constraints)
    %
    % NOTES: (1) In higher dimensions, duplicate constraints can
    %            appear. This program detects duplicates at up to 6
    %            digits of precision, then returns the unique constraints.
    %        (2) See companion function CON2VERT.
    %        (3) ver 1.0: initial version, June 2005.
    %        (4) ver 1.1: enhanced redundancy checks, July 2005
    %        (5) Written by Michael Kleder, 
    %
    %Modified by Matt Jacobson - March 29,2011
    % 
    k = convhulln(V);
    c = mean(V(unique(k),:));
    V = bsxfun(@minus,V,c);
    A  = nan(size(k,1),size(V,2));
    dim=size(V,2);
    ee=ones(size(k,2),1);
    rc=0;
    for ix = 1:size(k,1)
        F = V(k(ix,:),:);
        if lindep(F,tol) == dim
            rc=rc+1;
            A(rc,:)=F\ee;
        end
    end
    A=A(1:rc,:);
    b=ones(size(A,1),1);
    b=b+A*c';
    % eliminate duplicate constraints:
    [A,b]=rownormalize(A,b);
    [discard,I]=unique( round([A,b]*1e6),'rows');
    A=A(I,:); % NOTE: rounding is NOT done for actual returned results
    b=b(I);
    return
end
  
function [A,b]=rownormalize(A,b)
    %Modifies A,b data pair so that norm of rows of A is either 0 or 1
    
    if isempty(A), return; end
    
    normsA=sqrt(sum(A.^2,2));
    idx=normsA>0;
    A(idx,:)=bsxfun(@rdivide,A(idx,:),normsA(idx));
    b(idx)=b(idx)./normsA(idx);    
end
                  
 
function [r,idx,Xsub]=lindep(X,tol)
    %Extract a linearly independent set of columns of a given matrix X
    %
    %    [r,idx,Xsub]=lindep(X)
    %
    %in:
    %
    %  X: The given input matrix
    %  tol: A rank estimation tolerance. Default=1e-10
    %
    %out:
    %
    % r: rank estimate
    % idx:  Indices (into X) of linearly independent columns
    % Xsub: Extracted linearly independent columns of X
    if ~nnz(X) %X has no non-zeros and hence no independent columns
       Xsub=[]; idx=[];
       return
    end
    if nargin<2, tol=1e-10; end
    
           
     [Q, R, E] = qr(X,0); 
     
     diagr = abs(diag(R));
     %Rank estimation
     r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
     if nargout>1
      idx=sort(E(1:r));
        idx=idx(:);
     end
     
     
     if nargout>2
      Xsub=X(:,idx);                      
     end     
end

function [V,nr,nre]=lcon2vert(A,b,Aeq,beq,TOL,checkbounds)
%An extension of Michael Kleder's con2vert function, used for finding the 
%vertices of a bounded polyhedron in R^n, given its representation as a set
%of linear constraints. This wrapper extends the capabilities of con2vert to
%also handle cases where the  polyhedron has zero volume in R^n, i.e., where the
%polyhedron is defined by both equality and inequality constraints.
% 
%SYNTAX:
%
%  [V,nr,nre]=lcon2vert(A,b,Aeq,beq,TOL)
%
%The rows of the N x n matrix V are a series of N vertices of the polyhedron
%in R^n, defined by the linear constraints
%  
%   A*x  <= b
%   Aeq*x = beq
%
%By default, Aeq=beq=[], implying no equality constraints. The output "nr"
%lists non-redundant inequality constraints, and "nre" lists non-redundant 
%equality constraints.
%
%The optional TOL argument is a tolerance used for both rank-estimation and 
%for testing feasibility of the equality constraints. Default=1e-10. 
%The default can also be obtained by passing TOL=[];
%
%NOTE: It is important that the region specified by the inequality system A*x<=b
%have non-zero volume in R^n. For example, A=b=[1;-1] is not legal input data, 
%because the only solution to A*x<=b is x=1, which has zero volume in R^1. The
%proper way to express a zero-volume region is with the addition of
%equality constraint data, as for example Aeq=1, beq=1, A=1,b=100.
%
%EXAMPLE: 
%
%The 3D region defined by x+y+z=1, x>=0, y>=0, z>=0
%is described by the following constraint data.
% 
%
%     A =
% 
%         0.4082   -0.8165    0.4082
%         0.4082    0.4082   -0.8165
%        -0.8165    0.4082    0.4082
% 
% 
%     b =
% 
%         0.4082
%         0.4082
%         0.4082
% 
% 
%     Aeq =
% 
%         0.5774    0.5774    0.5774
% 
% 
%     beq =
% 
%         0.5774
%
%
%  >> V=lcon2vert(A,b,Aeq,beq)
%
%         V =
% 
%             1.0000    0.0000    0.0000
%             0.0000    0.0000    1.0000
%            -0.0000    1.0000    0.0000
%
%
  %%initial argument parsing
  
  nre=[];
  nr=[];
  if nargin<5 || isempty(TOL), TOL=1e-10; end
  if nargin<6, checkbounds=true; end
  
  switch nargin 
      
      case 0
          
           error 'At least 1 input argument required'
       
      case 1
        
         b=[]; Aeq=[]; beq=[]; 
        
          
      case 2
          
          Aeq=[]; beq=[];
          
      case 3
          
          beq=[];
          error 'Since argument Aeq specified, beq must also be specified'
            
  end
  
  
  b=b(:); beq=beq(:);
  
  if xor(isempty(A), isempty(b)) 
     error 'Since argument A specified, b must also be specified'
  end
      
  if xor(isempty(Aeq), isempty(beq)) 
        error 'Since argument Aeq specified, beq must also be specified'
  end
  
  
  nn=max(size(A,2)*~isempty(A),size(Aeq,2)*~isempty(Aeq));
  
  if ~isempty(A) && ~isempty(Aeq) && ( size(A,2)~=nn || size(Aeq,2)~=nn)
      
      error 'A and Aeq must have the same number of columns if both non-empty'
      
  end
  
  
  inequalityConstrained=~isempty(A);  
  equalityConstrained=~isempty(Aeq);
 [A,b]=rownormalize(A,b);
 [Aeq,beq]=rownormalize(Aeq,beq);
 
  if equalityConstrained && nargout>2
 
        
        [discard,nre]=lindep([Aeq,beq].',TOL);  %#ok<ASGLU>
          
        if ~isempty(nre) %reduce the equality constraints
            
            Aeq=Aeq(nre,:);
            beq=beq(nre);
            
        else    
            equalityConstrained=false;
        end
        
   end
      
  
   %%Find 1 solution to equality constraints within tolerance
  
            
   if equalityConstrained
        
        
       Neq=null(Aeq);   
       x0=pinv(Aeq)*beq;
       if norm(Aeq*x0-beq)>TOL*norm(beq)  %infeasible
          nre=[]; nr=[]; %All constraints redundant for empty polytopes
          V=[]; 
          return;
          
       
       elseif isempty(Neq)
           
           if inequalityConstrained && ~all(A*x0<=b)
            
              nre=[]; nr=[]; %All constraints redundant for empty polytopes
              V=[]; 
              return;              
               
           else %inequality constraints all satisfied, including vacuously
               
               V=x0(:).'; 
               nre=(1:nn).'; %Equality constraints determine everything. 
               nr=[];%All inequality constraints are therefore redundant.             
               return
           
           end
           
           
       end
 
       %rkAeq= nn - size(Neq,2);
       
       
  end  
   
    %%
  if inequalityConstrained && equalityConstrained
     
   AAA=A*Neq;
   bbb=b-A*x0;
    
  elseif inequalityConstrained
      
    AAA=A;
    bbb=b;
   
  elseif equalityConstrained && ~inequalityConstrained
      
       error('Non-bounding constraints detected. (Consider box constraints on variables.)')
      
    
  end
  
  nnn=size(AAA,2);
  
  if nnn==1 %Special case
      
     idxu=sign(AAA)==1;
     idxl=sign(AAA)==-1;
     idx0=sign(AAA)==0;
     
     Q=bbb./AAA;
     U=Q; 
       U(~idxu)=inf;
     L=Q;
       L(~idxl)=-inf;
     
     [ub,uloc]=min(U);
     [lb,lloc]=max(L);
     
     if ~all(bbb(idx0)>=0) || ub<lb %infeasible
         
         V=[]; nr=[]; nre=[];
         return
         
     elseif ~isfinite(ub) || ~isfinite(lb)
         
         error('Non-bounding constraints detected. (Consider box constraints on variables.)')
         
     end
      
     Zt=[lb;ub];
     
     if nargout>1
        nr=unique([lloc,uloc]); nr=nr(:);
     end
     
      
  else    
      
          if nargout>1
           [Zt,nr]=con2vert(AAA,bbb,TOL,checkbounds);
          else
            Zt=con2vert(AAA,bbb,TOL,checkbounds); 
          end
  
  end
  
  if equalityConstrained && ~isempty(Zt)
     
      V=bsxfun(@plus,Zt*Neq.',x0(:).'); 
      
  else
      
      V=Zt;
      
  end
 
  if isempty(V)
     nr=[]; nre=[]; 
  end
  
end

function [V,nr] = con2vert(A,b,TOL,checkbounds)
% CON2VERT - convert a convex set of constraint inequalities into the set
%            of vertices at the intersections of those inequalities;i.e.,
%            solve the "vertex enumeration" problem. Additionally,
%            identify redundant entries in the list of inequalities.
% 
% V = con2vert(A,b)
% [V,nr] = con2vert(A,b)
% 
% Converts the polytope (convex polygon, polyhedron, etc.) defined by the
% system of inequalities A*x <= b into a list of vertices V. Each ROW
% of V is a vertex. For n variables:
% A = m x n matrix, where m >= n (m constraints, n variables)
% b = m x 1 vector (m constraints)
% V = p x n matrix (p vertices, n variables)
% nr = list of the rows in A which are NOT redundant constraints
% 
% NOTES: (1) This program employs a primal-dual polytope method.
%        (2) In dimensions higher than 2, redundant vertices can
%            appear using this method. This program detects redundancies
%            at up to 6 digits of precision, then returns the
%            unique vertices.
%        (3) Non-bounding constraints give erroneous results; therefore,
%            the program detects non-bounding constraints and returns
%            an error. You may wish to implement large "box" constraints
%            on your variables if you need to induce bounding. For example,
%            if x is a person's height in feet, the box constraint
%            -1 <= x <= 1000 would be a reasonable choice to induce
%            boundedness, since no possible solution for x would be
%            prohibited by the bounding box.
%        (4) This program requires that the feasible region have some
%            finite extent in all dimensions. For example, the feasible
%            region cannot be a line segment in 2-D space, or a plane
%            in 3-D space.
%        (5) At least two dimensions are required.
%        (6) See companion function VERT2CON.
%        (7) ver 1.0: initial version, June 2005
%        (8) ver 1.1: enhanced redundancy checks, July 2005
%        (9) Written by Michael Kleder
%
%Modified by Matt Jacobson - March 30, 2011
% 
%%%3/4/2012 Improved boundedness test - unfortunately slower than Michael Kleder's
if checkbounds
   
[aa,bb,aaeq,bbeq]=vert2lcon(A,TOL);

if any(bb<=0) || ~isempty(bbeq)
    error('Non-bounding constraints detected. (Consider box constraints on variables.)')
end

clear aa bb aaeq bbeq

end

dim=size(A,2);

%%%Matt J initialization
if strictinpoly(b,TOL)  
   
   c=zeros(dim,1);

else

        
        slackfun=@(c)b-A*c;
        %Initializer0
        c = pinv(A)*b; %02/17/2012 -replaced with pinv()
        s=slackfun(c);
        if ~approxinpoly(s,TOL) %Initializer1
            c=Initializer1(TOL,A,b,c);
            s=slackfun(c);
        end
        if  ~approxinpoly(s,TOL)  %Attempt refinement
            %disp 'It is unusually difficult to find an interior point of your polytope. This may take some time... '
            %disp ' '   
            c=Initializer2(TOL,A,b,c);
            %[c,fval]=Initializer1(TOL,A,b,c,10000);
            s=slackfun(c);
        end
        if ~approxinpoly(s,TOL)
                %error('Unable to locate a point near the interior of the feasible region.')
                V=[];
                nr=[];
                return
        end
       if ~strictinpoly(s,TOL) %Added 02/17/2012 to handle initializers too close to polytope surface
            %disp 'Recursing...'
            idx=(  abs(s)<=max(s)*TOL );
            Amod=A; bmod=b; 
             Amod(idx,:)=[]; 
             bmod(idx)=[];
            Aeq=A(idx,:); %pick the nearest face to c
            beq=b(idx);
            faceVertices=lcon2vert(Amod,bmod,Aeq,beq,TOL,1);
            if isempty(faceVertices)
               disp 'Something''s wrong. Couldn''t find face vertices. Possibly polyhedron is unbounded.'
               keyboard
            end
            c=faceVertices(1,:).';  %Take any vertex - find local recession cone vector
            s=slackfun(c);
            idx=(  abs(s)<=max(s)*TOL );
            Asub=A(idx,:); bsub=b(idx,:);
            [aa,bb,aaeq,bbeq]=vert2lcon(Asub);
            aa=[aa;aaeq;-aaeq];
            bb=[bb;bbeq;-bbeq];
            clear aaeq bbeq
            [bmin,idx]=min(bb);
             if bmin>=-TOL
               disp 'Something''s wrong. We should have found a recession vector (bb<0).'
               keyboard
             end      
            Aeq2=null(aa(idx,:)).';
            beq2=Aeq2*c;  %find intersection of polytope with line through facet centroid.
            linetips = lcon2vert(A,b,Aeq2,beq2,TOL,1);
            if size(linetips,1)<2
               disp 'Failed to identify line segment through interior.'
               disp 'Possibly {x: Aeq*x=beq} has weak intersection with interior({x: Ax<=b}).'
               keyboard
            end
            lineCentroid=mean(linetips);%Relies on boundedness
            clear aa bb
            c=lineCentroid(:);
            s=slackfun(c);
        end
        b = s;
end
%%%end Matt J initialization


D=bsxfun(@rdivide,A,b); 


k = convhulln(D);
nr = unique(k(:));



G  = zeros(size(k,1),dim);
ee=ones(size(k,2),1);
discard=false( 1, size(k,1) );

for ix = 1:size(k,1) %02/17/2012 - modified
    
    F = D(k(ix,:),:);
    if lindep(F,TOL)<dim; 
        discard(ix)=1;
        continue; 
    end
    G(ix,:)=F\ee;
    
end

G(discard,:)=[];

V = bsxfun(@plus, G, c.'); 

[discard,I]=unique( round(V*1e6),'rows');
V=V(I,:);

return
 end

function [c,fval]=Initializer1(TOL, A,b,c,maxIter)
    
    thresh=-10*max(eps(b));
    
    if nargin>4
     [c,fval]=fminsearch(@(x) max([thresh;A*x-b]), c,optimset('MaxIter',maxIter));
    else
     [c,fval]=fminsearch(@(x) max([thresh;A*x-b]), c); 
    end
    
    return   
end

function c=Initializer2(TOL,A,b,c)
     %norm(  (I-A*pinv(A))*(s-b) )  subj. to s>=0 
        
        maxIter=100000;
     
        [mm,nn]=size(A);
        
        
        
        
         Ap=pinv(A);        
         Aaug=speye(mm)-A*Ap;
         Aaugt=Aaug.';
        
        M=Aaugt*Aaug;
        C=sum(abs(M),2);
         C(C<=0)=min(C(C>0));
        
        slack=b-A*c;
        slack(slack<0)=0;
        
         
            %     relto=norm(b);
            %     relto =relto + (relto==0); 
            %     
            %      relres=norm(A*c-b)/relto;
         
        IterThresh=maxIter; 
        s=slack; 
        ii=0;
        %for ii=1:maxIter
        while ii<=2*maxIter %HARDCODE
            
           ii=ii+1; 
           if ii>IterThresh, 
               %warning 'This is taking a lot of iterations'
               IterThresh=IterThresh+maxIter;
           end          
              
         s=s-Aaugt*(Aaug*(s-b))./C;   
         s(s<0)=0;
          
           c=Ap*(b-s);
           %slack=b-A*c;
           %relres=norm(slack)/relto;
           %if all(0<slack,1)||relres<1e-6||ii==maxIter, break;  end
           
        end
       
    return 
end
        
function tf=approxinpoly(s,TOL)
 
    smax=max(s);
    
    if smax<=0
      tf=false; return 
    end
    
    tf=all(s>=-smax*TOL);
end
   
function tf=strictinpoly(s,TOL)
  
    smax=max(s);
    
    if smax<=0
      tf=false; return 
    end
    
    tf=all(s>=smax*TOL);
         
end