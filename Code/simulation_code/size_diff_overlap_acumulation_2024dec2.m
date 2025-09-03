% Thomas C. Day, 2022 edit
% This script generates N clusters of snowflake yeast, modeled as connected
% cells, where each cell is a prolate ellipsoid. The stop condition of the
% experiment can either be a number of generations or a total overlap stop
% condition. Cell aspect ratio can also be varied. Last, there are options
% to include force relaxation to cell positions and orientations.

%{
% Luis Felipe 2023
% this code is based on snowflake_sim_growth_from_python.m
The modifications for this code are in ELYES_SIM and the main loop and that a lot of the plot 
functions are missing.

This code allows for the clusters to grow until they reach a specific size as the overlap threshold is
setup to be really high. In addition it also has the option to run the force relaxation function after 
each cell is added.

To execute this code use the bash script execute_size_overlap_through_directory.sh

Note: 
to run the results of overlaps_cluster_size_22nov2023 you need to make include_forces_each_step=0
to run the results of overlaps_force_each_step_28nov2023 you need to make include_forces_each_step=1

% September 26 2024
If after 10 attempts it hasn't found how a place to position the new bud then it will just add it in a random position
Added to also calculate the volume of the cluster, and diameter of gyration


% 4nov2024
I deleted the include forces functions and the get cell overlap functions, and now the code is running much faster, so at some
point this function were being run and that is why the code was being executed slower than expected, so now it is garanteed that
the function CHECK_OVERLAPS is the only one being used to calculate the cluster overlaps

%}

%% Inputs:
folder = pwd;
cd(folder);


% strain="petite";
strain="grande";

%directory paths
folder_path='~/work_dir/observed_synchrony/paper_results_2dec2024/physics_sim/overlap_accumulation_2dec2024/';


% Network directory
directory_path = folder_path+"test_"+strain+"_1000m_100n/edges_sim_files/";



% Open output file
% output_folder=folder_path;
output_folder='~/Desktop/results_edge_degree_15/supp_fig3_physics_sim';


clust_filename = strain+"_clust_size_50sim_500n_1.2aspr.csv";

cluster_size_file_name=fullfile(output_folder, clust_filename);

write_to_output_file=0; % 1=YES, 0=NO # For size at fracture calculations

%# Calculate individual cells overlap? This is for the main loop
% calc_ind_overlap=1; % 1 = Y, 0 = N


%# save individual cell overlaps to file? this is in the include forces calculation
save_ind_overlaps=0; % 1 = Y, 0 = N
overl_file = strain+"_cell_overlap_1c_2r_1.2ar.csv";
overlaps_filename = fullfile(output_folder, overl_file);

tot_overl_file=strain+"_total_overlap_1c_2r_1.2ar.csv";
tot_overlaps_filename=fullfile(output_folder, tot_overl_file);


%# save total overlap volume for each network size
save_overlap_each_size=1; % 1 = Y, 0 = N
after_how_many_divisions=5;
overlap_size_file=strain+"_overlap_clust_size_1.2ar_10attempts.csv";
overlap_size_filename=fullfile(output_folder, overlap_size_file);

%# save each cell overlap for each network size
save_cell_overlap_each_size=0; % 1 = Y, 0 = N
cell_overlap_size_file=strain+"_cell_overlap_clust_size_1.2ar.csv";
cell_overlap_size_filename=fullfile(output_folder, cell_overlap_size_file);


% Cell budding angle information
use_cell_angles = 1;

cell_angle_path = '/Users/pipe/emergence_of_coordinated_cell_division_during_the_evolution_of_multicellularity/Data/cell_budding_angles_grande_ancestor_30july2025.csv';
cell_angle_data = readtable(cell_angle_path, 'Delimiter', ',', 'ReadVariableNames', true);
cell_angles_degrees = cell_angle_data.angle_to_bud_scar;
cell_angles_radians = deg2rad(cell_angles_degrees);  % Convert to radians


max_network_size=500;
add_cell_after_all_attempts=0; % 1 it always add the cell, 0 only adds it if there was space in the first 10 attempts

N               = 30; % number of clusters to generate
numGens         = 5; % number of generations of cell division in a group
diam            = 4.58; % smallest diameter of the cell
err_diam        = 0; % variation in cell size, taken from coefficient of variation data
AR              = 1.2;%:.1:2.8; % aspect ratio value %# DEFAULT was 1.2
err_AR          = 0; % standard deviation in aspect ratio from Shane's measurements
pole_theta      = deg2rad(10); % buds nearest the pole will be chosen from between 0 and 10 degrees in polar angle
THETA           = deg2rad(45); % polar angle average from SEM data
thetaVariance   = deg2rad(0); % variation in polar angle from SEM data
distance_thresh = 1.1672; % minimum distance (um) separating bud scars
new_bud_prob    = .8; % probability that the first cell will bud near the pole
check_overlap   = 1; % do we check the overlap?
overlap_thresh  = 4.0e9; % threshold of total overlaps %#2.5e2 
%really high value for the overlap threshold so that it never reaches it

% Force input parameters:
include_forces  = 0; % do we include forces for rearrangements? 1 = Y, 0 = N
NEIGHBOR_THRESH = 3*diam; % radius at which you consider cells "neighbors"
STERIC_MAG      = 2; % relative strength of steric interactions(from cells colliding)
CHITIN_MAG      = 1; % relative strength of the bonds holding cells together
BOND_TORQ_MAG   = 1; % relative torsional strength of chitin bonds
mobility_pos    = .1;
mobility_rot    = .1;
dt              = .1; % timestep for mechanical relaxation
T               = 10; % number of timesteps to allow mechanical relaxation

% Do we visualize figures?
figure_viz = 0; % 1 = Y, 0 = N

%# Include forces after each cell added
include_forces_each_step=1; % 1 = Y, 0 = N
T_each_step=1; % This should always be one

%# Creating file for saving size at fracture distributions
clust_size_output=0; %#initializing variable to avoid errors even if it is not used
if write_to_output_file==1
    clust_size_output = create_file(cluster_size_file_name);
    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if clust_size_output == -1
        return;
    end
    % Add file header
    fprintf(clust_size_output, 'file_num,sim_number,size_fracture,volume,gyration_diam\n');
end

%# Creating file for saving cell overlap information
overlaps_output=0; %#initializing variable to avoid errors even if it is not used
if save_ind_overlaps==1
    overlaps_output = create_file(overlaps_filename);
    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if overlaps_output == -1
        return;
    end
    % Add file header
    fprintf(overlaps_output, 'file_num,sim_number,node_id,degree,overlap_vol,num_small_overlaps,num_big_overlaps,relax_step\n');
end

%# creating file for the total overlap information
total_overlaps_output=0; %#initializing variable to avoid errors even if it is not used
if save_ind_overlaps==1
    total_overlaps_output = create_file(tot_overlaps_filename);
    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if total_overlaps_output == -1
        return;
    end
    % Add file header
    fprintf(total_overlaps_output, 'file_num,sim_number,relax_step,volume,cluster_size,gyration_diam\n');
end

%# creating file for overlap for each cluster size
overlap_size_output=0; %#initializing variable to avoid errors even if it is not used
if save_overlap_each_size==1
    overlap_size_output = create_file(overlap_size_filename);
    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if overlap_size_output == -1
        return;
    end
    % Add file header
    fprintf(overlap_size_output, 'file_num,sim_number,cluster_size,overlap_vol,not_added_nodes\n');
end

%# creating file for overlap for each cluster size
cell_overlap_size_output=0; %#initializing variable to avoid errors even if it is not used
if save_cell_overlap_each_size==1
    cell_overlap_size_output = create_file(cell_overlap_size_filename);
    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if cell_overlap_size_output == -1
        return;
    end
    % Add file header
    fprintf(cell_overlap_size_output, 'file_num,sim_number,node_id,degree,cluster_size,overlap_vol,num_small_overlaps,num_big_overlaps\n');
end


% Creating input variables dictionary
global input_variables;
input_variables=containers.Map();
input_variables('clust_size_output')=clust_size_output;
input_variables('overlaps_output')=overlaps_output;
input_variables('aspect_ratio')=AR;
input_variables('save_ind_overlaps')=save_ind_overlaps;
input_variables('write_to_output_file')=write_to_output_file;
input_variables('total_overlaps_output')=total_overlaps_output;
input_variables('neighbor_thresh')=NEIGHBOR_THRESH;
input_variables('max_network_size')=max_network_size;
input_variables('file_overlap_size')=overlap_size_output;
input_variables('save_overlap_each_size')=save_overlap_each_size;
input_variables('include_forces_each_step')=include_forces_each_step;
input_variables('T_each_step')=T_each_step;
input_variables('save_cell_overlap_each_size')=save_cell_overlap_each_size;
input_variables('cell_overlap_size_output')=cell_overlap_size_output;
input_variables('after_how_many_divisions')=after_how_many_divisions;
input_variables('add_cell_after_all_attempts')=add_cell_after_all_attempts;
input_variables('use_cell_angles')=use_cell_angles;
input_variables('cell_angles_radians')=cell_angles_radians;

%# inputs necesary to execute include_forces function in the ELYES_sim function
input_variables('NEIGHBOR_THRESH')=NEIGHBOR_THRESH;
input_variables('STERIC_MAG')=STERIC_MAG;
input_variables('CHITIN_MAG')=CHITIN_MAG;
input_variables('BOND_TORQ_MAG')=BOND_TORQ_MAG;
input_variables('mobility_pos')=mobility_pos;
input_variables('mobility_rot')=mobility_rot;
input_variables('dt')=dt;
input_variables('figure_viz')=figure_viz;

% Initialization:
fprintf('Welcome to the cell simulator v2.0\n');
fprintf('Written by Thomas C. Day, 2020\n');

%% MAIN CODE: -------------------------------------------------------------
% -------------------------------------------------------------------------
% Parallel iteration over varied params:
% save('params.mat','N','numGens','diam','pole_theta','THETA','distance_thresh','new_bud_prob','overlap_thresh');
% Ncells = zeros(length(AR),1);
% Scells = Ncells;
% for a = 1:length(AR)

file_list = dir(fullfile(directory_path, 'edges_network_sim_*.csv'));
tot_files=length(file_list);

for i = 1:length(file_list)

    % Get the current file's name
    file_name = file_list(i).name;
    
    % Use regular expressions to extract the number from the file name
    match = regexp(file_name, 'edges_network_sim_(\d+)\.csv', 'tokens');

    full_file_path = fullfile(directory_path, file_name);

    edges_data = readtable(full_file_path, 'Delimiter', ',', 'ReadVariableNames', true);
    node1_vals      = edges_data.node1;
    node2_vals      = edges_data.node2;

    if ~isempty(match)
        % Extracted number as a string
        number_str = match{1}{1};
        
        % Convert the number string to a numeric value
        file_number = str2double(number_str);

        input_variables('file_number')=file_number;

        cells_sim = cell(N,1);
        % aspRat    = AR(a);
        aspRat = AR;
        theta     = THETA;
        filename = ['sim','_AspRat=',num2str(aspRat)];
        filename = strrep(filename, '.','-');
        COM = zeros(3,N);
        RAD = zeros(1,N);
    	ncells = zeros(1,N);
    	
        for ii = 1:N
            
            input_variables('sim_number')=ii;

            % Generate one cluster:
            % fprintf(['Cluster #: ',num2str(ii), ' / ', num2str(N), '\n']);
            [cell_list] = ELYES_SIM(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, numGens, NEIGHBOR_THRESH, check_overlap, overlap_thresh, node1_vals, node2_vals);
            ncells(ii) = length(cell_list);
    		
            % Allow for force relaxations from steric interactions:
            if include_forces == 1
                warning('off');
                [~, cell_list, overlaps] = INCLUDE_FORCES(cell_list, NEIGHBOR_THRESH, STERIC_MAG, CHITIN_MAG, BOND_TORQ_MAG, mobility_pos, mobility_rot, dt, T, figure_viz);
                warning('on');
            end

            % Calculating individual cell overlaps
            % if calc_ind_overlap == 1
            %     [cell_overlaps, big_overlaps_list, small_overlaps_list] = calculate_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);

            %     if save_ind_overlaps==1
            %         save_ind_overlaps_to_file(overlaps_output, file_number, ii, cell_list, cell_overlaps, big_overlaps_list, small_overlaps_list, 0);
            %     end
            % end
     
            % Save data to file:
    		%{
            if mod(ii,10) == 0
                parsave([filename,'.mat'], cells_sim);
            end
    		%}
    		
            
    		
        end
    	% fprintf(['Mean number of cells per group: ',num2str(mean(ncells)),'\n']);
    	% fprintf(['Stdev number of cells per group: ',num2str(std(ncells)),'\n']);
    	% Ncells(a) = mean(ncells);
    	% Scells(a) = std(ncells);
    end

    fprintf(['File #: ',num2str(i), ' / ', num2str(tot_files), '\n']);
end

% Show figure of one cluster:
if figure_viz == 1
    figure;
    hold on; box on; set(gca,'linewidth',2);
    for k = 1:length(cell_list)
        [x,y,z] = VISUALIZE_ELLIPSOID(cell_list(k), 30);
        surf(x,y,z,'facealpha',1,'edgecolor','none');
    end
    view(3); axis equal;
    lighting gouraud;
    lightangle(0,30);
    material dull;
end



%# Close output files
if write_to_output_file==1
    fclose(clust_size_output);
end

if save_ind_overlaps==1
    fclose(overlaps_output);
end

if save_ind_overlaps==1
    fclose(total_overlaps_output);
end

%% Added functions by Luis

function fileID = create_file(fileName)
    % Check if the file exists
    if isfile(fileName)
        % File exists, so display a message and stop execution
        fprintf('The file "%s" already exists. Simulation results are already available.\n', fileName);
        fileID = -1; % Return -1 to indicate that the file already exists
        return;
    end

    % If the code reaches this point, the file does not exist, so create it for writing
    fileID = fopen(fileName, 'w');

    % Check if the file was successfully opened
    if fileID == -1
        error('Error opening file for writing');
    end
end


%# Function used to calculate each cells overlap volume with their neighbors
% This function iterates through each cell and find it neighbors to calculate the overlaps
% for that cell and the vector final_overlaps is going to be of the same size as cell_list 
% so the indexes are going to map for each cell overlap
% This function was modified from get_forces_torques
function [final_overlaps, big_overlaps_list, small_overlaps_list] = calculate_ind_cell_overlaps(cell_list, neighbor_thresh)
    % Calculate the forces and torques acting on all cells in the group.
    % This function can get quite heavy.

    % Obtain distances between all particles:
    [D, ~, ~, ~] = get_particle_distances(cell_list);
    
    % Find all pairs of cells within the neighborhood threshold:
    N = D < neighbor_thresh; % finds all neighbor cells
    N = N - eye(size(N)); % we don't care about counting a cell as its own neighbor
    

    final_overlaps = [];
    big_overlaps_list=[];
    small_overlaps_list=[];
    for n = 1:length(cell_list)
        neighbors = find(N(:,n) == 1); % only need to check the neighboring cells for interactions
        
        % Steric interactions:
        [temp_cell_overlap, num_big_overlaps, num_small_overlaps] = get_cell_overlap(cell_list, n, neighbors);
        final_overlaps = [final_overlaps, temp_cell_overlap];
        big_overlaps_list = [big_overlaps_list, num_big_overlaps];
        small_overlaps_list = [small_overlaps_list, num_small_overlaps];
    end
end



function save_ind_overlaps_to_file(cell_list, cell_overlaps, big_overlaps_list, small_overlaps_list, relax_step)

    global input_variables;
    file_number=input_variables('file_number');
    sim_number=input_variables('sim_number');
    overlaps_output=input_variables('overlaps_output');


    for i = 1:length(cell_list)
        %cell id (network)
        temp_cell_id=cell_list(i).networkID;

        %The plus one is to take into account their connection to their mother
        temp_degree=length(cell_list(i).Daughters)+1;

        fprintf(overlaps_output, '%d,%d,%d,%d,%d,%d,%d,%d\n', file_number, sim_number, temp_cell_id, temp_degree, cell_overlaps(i), big_overlaps_list(i), small_overlaps_list(i), relax_step);
    end
end


function save_ind_overlaps_to_file_each_step(cell_list, cell_overlaps, big_overlaps_list, small_overlaps_list, curr_network_size, cont_failed_to_add_edge)

    global input_variables;
    file_number=input_variables('file_number');
    sim_number=input_variables('sim_number');
    cell_overlap_size_output=input_variables('cell_overlap_size_output');


    for i = 1:length(cell_list)
        %cell id (network)
        temp_cell_id=cell_list(i).networkID;

        %The plus one is to take into account their connection to their mother
        temp_degree=length(cell_list(i).Daughters)+1;

        fprintf(cell_overlap_size_output, '%d,%d,%d,%d,%d,%d,%d,%d\n', file_number, sim_number, temp_cell_id, temp_degree, curr_network_size, cell_overlaps(i), big_overlaps_list(i), small_overlaps_list(i));
    end
end

function save_tot_overlap(total_overlap, relax_step, cluster_size)
    global input_variables;
    file_number=input_variables('file_number');
    sim_number=input_variables('sim_number');
    total_overlaps_output=input_variables('total_overlaps_output');

    fprintf(total_overlaps_output, '%d,%d,%d,%d,%d\n', file_number, sim_number, relax_step, total_overlap, cluster_size);

end


% Function to obtain the gyration diameter of a cluster, it calculates the radius, converts it to meters and then
%multiplies by two to get the diameter
function [diameter]=gyration_diameter(cell_list)
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
    diameter=Rg*2;
end



%% FUNCTIONS
function [cell_list] = ELYES_SIM(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, numGens, neighbor_thresh, check_overlap, overlap_thresh, node1_vals, node2_vals)

    global input_variables;
    write_to_output_file=input_variables('write_to_output_file');
    sim_number=input_variables('sim_number');
    file_number=input_variables('file_number');
    clust_size_output=input_variables('clust_size_output');
    max_network_size=input_variables('max_network_size');
    file_overlap_size=input_variables('file_overlap_size');
    save_overlap_each_size=input_variables('save_overlap_each_size');
    include_forces_each_step=input_variables('include_forces_each_step');
    save_cell_overlap_each_size=input_variables('save_cell_overlap_each_size');
    after_how_many_divisions=input_variables('after_how_many_divisions');
    add_cell_after_all_attempts=input_variables('add_cell_after_all_attempts');
    use_cell_angles=input_variables('use_cell_angles');
    cell_angles_radians=input_variables('cell_angles_radians');

    NEIGHBOR_THRESH=input_variables('NEIGHBOR_THRESH');
    STERIC_MAG=input_variables('STERIC_MAG');
    CHITIN_MAG=input_variables('CHITIN_MAG');
    BOND_TORQ_MAG=input_variables('BOND_TORQ_MAG');
    mobility_pos=input_variables('mobility_pos');
    mobility_rot=input_variables('mobility_rot');
    dt=input_variables('dt');
    T_each_step=input_variables('T_each_step');
    figure_viz=input_variables('figure_viz');

    cell_list=[];  % this will eventually be the tabulated list
    
    % Create root cell:
    rootCell.Center       = [0; 0; 0];
    rootCell.Radii        = 1/2 * [aspRat * diam; diam; diam];
    rootCell.Rmatrix      = eye(3);
    rootCell.Generation   = 0;
    rootCell.IDnumber     = 0;
    rootCell.networkID    = 1; %#
    rootCell.Daughters    = [];
    rootCell.Mother       = [];
    rootCell.BudXYZ       = [];
    rootCell.Overlaps     = [];
    rootCell.Neighborhood = [];
    cell_list             = [cell_list, rootCell];
    counting_ix           = 0; % this will count how many cells are in the cluster

    cont_failed_to_add_edge=0;
    cont_missing_mother_cell=0;

    last_network_size=1;

    % Iterate over generations:
    % for g = 1:numGens
    % Iterate over the edges of the loaded file
    num_rows_edges=length(node1_vals);
    for cont_rows = 1:num_rows_edges

        curr_network_size=length(cell_list);
        
        % First check the total amount of overlap in the cluster:edges_sim_files
        %fprintf(['Gen = ',num2str(g),'\n']);
        if check_overlap == 1
            [flag, total_overlap] = CHECK_OVERLAPS2(cell_list, overlap_thresh);


            % saving the overlap for each size and breaking the loop if the size reaches the defined maximum
            % the value is only saved if a cell was added to the cluster
            if save_overlap_each_size==1 & mod(curr_network_size, after_how_many_divisions)==0
                
                if curr_network_size>last_network_size

                    %# total overlap
                    fprintf(file_overlap_size, '%d,%d,%d,%d,%d\n', file_number,sim_number,curr_network_size,total_overlap, cont_failed_to_add_edge);

                    %# cell overlaps
                    if save_cell_overlap_each_size==1
                        [cell_overlaps, big_overlaps_list, small_overlaps_list] = calculate_ind_cell_overlaps(cell_list, neighbor_thresh);
                        save_ind_overlaps_to_file_each_step(cell_list, cell_overlaps, big_overlaps_list, small_overlaps_list, curr_network_size);
                    end
                    last_network_size=curr_network_size;
                end
                if last_network_size>=max_network_size
                    break
                end
            end
        else
            flag = 0;
        end
        
        % disp(flag);

        if flag ~= 1
            

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

            if cell_id_in_list==1

                %# Doing at least 10 attempts to find where the cell is going to be located
                cont_attempts=0;
                newBud=0; %# Default value of new bud, it is assigned here as this would be the first time it is assigned

                if use_cell_angles==0
                    % Sampling buds without using cell budding distributions
                    force_selection=0;
                    while cont_attempts<10 & newBud == 0 %if one of this conditions is met then it needs to stop
                        %# This may need to be outside the loop, now we are only sampling 1 value because it is processing one cell at a time
                        variedTheta = theta + thetaVariance * randn(1); %theta value corresponding to the cell
                        %# To test if variedtheta gets different values or not when
                        %# thetaVariance is 0
                        % fprintf(['variedTheta: ', num2str(variedTheta), '\n']); %# it prints always the same value

                        [newBud, newBudRel, newAxis] = getDaughterPos(cell_list(n), pole_theta, variedTheta, distance_thresh, new_bud_prob, force_selection); % finds the new bud xyz position and relative position to the old cell
                        cont_attempts=cont_attempts+1;
                    end

                    % If there wasn't space for the cell then it is just sampled randomly if add_cell_after_all_attempts is 1
                    force_selection=1;
                    if newBud==0 & add_cell_after_all_attempts==1
                        variedTheta = theta + thetaVariance * randn(1);
                        [newBud, newBudRel, newAxis] = getDaughterPos(cell_list(n), pole_theta, variedTheta, distance_thresh, new_bud_prob, force_selection); % finds the new bud xyz position and relative position to the old cell
                    end
                else
                    % Sample buds using cell budding distributions
                    force_selection=0;
                    while cont_attempts<10 & newBud == 0 %if one of this conditions is met then it needs to stop
                        [newBud, newBudRel, newAxis] = getDaughterPos_measured_angles(cell_list(n), distance_thresh, cell_angles_radians, force_selection); % finds the new bud xyz position and relative position to the old cell
                    end

                    % If there wasn't space for the cell then it is just sampled randomly if add_cell_after_all_attempts is 1
                    force_selection=1;
                    if newBud==0 & add_cell_after_all_attempts==1
                        [newBud, newBudRel, newAxis] = getDaughterPos_measured_angles(cell_list(n), distance_thresh, cell_angles_radians, force_selection); % finds the new bud xyz position and relative position to the old cell
                    end
                end


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
            end
        else
            cell_list = cell_list;
            break;
        end
        
        % fprintf(['cont_rows: ',num2str(cont_rows),'\n']);

    end

    %# write information to output file
    if write_to_output_file==1
        % Calculating cluster volume
        Centers = [cell_list.Center];
        Centers = Centers';
        [T,V] = convhull(Centers);

        temp_diameter=gyration_diameter(cell_list);

        fprintf(clust_size_output, '%d,%d,%d,%d,%d\n', file_number, cont_sim, length(cell_list), V, temp_diameter);
    end
    
end

function [flag, total_overlap] = CHECK_OVERLAPS2(cell_list, overlap_thresh)

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
    total_overlap=U;
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


function [newPos, relativePos, axis] = getDaughterPos(CELL, POLE_TH, TH, DISTANCE_THRESH, NEW_BUD_PROB_THRESH, force_selection)
    % Surface definition for possible mother cell:
    [Sm, Rm, Tm] = GET_SURFACE_MATRICES(CELL, .1);
    
    % Determine where on cell body to bud the next scar:
    if isempty(CELL.Daughters) % reproduce near the pole with prob 0.7
        frac = rand;
        if frac < NEW_BUD_PROB_THRESH
            % bud at pole
            th = POLE_TH * rand;
            ph = 2 * pi * rand;
            r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2 ).^(-1));
            x = r * cos(th);
            y = r * sin(th) * cos(ph);
            z = r * sin(th) * sin(ph);
        else
            % bud on the polar angle
            ph = 2 * pi * rand;
            th = TH;
            r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2 ).^(-1));
            x = r * cos(th);
            y = r * sin(th) * cos(ph);
            z = r * sin(th) * sin(ph);
        end
        
        
        ix_too_close = [];  % No need to check distances as it is the first cell added
        
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
        
        % Calculate position
        r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2 ).^(-1));
        x = r * cos(th);
        y = r * sin(th) * cos(ph);
        z = r * sin(th) * sin(ph);
        
        % Apply force_selection logic for subsequent daughters
        if force_selection == 1
            ix_too_close = [];  % Force selection - ignore distance checks
        else
            % Check if new location is too close to existing scars:
            t = [x; y; z];
            s = CELL.BudXYZ - t;
            d = sqrt(s(1,:).^2 + s(2,:).^2 + s(3,:).^2);
            ix_too_close = find(d < DISTANCE_THRESH);
        end
    end
    
    if isempty(ix_too_close) % new bud is successful
        relativePos = [x; y; z];
        newPos = Tm * Rm * [relativePos; 1];
        normal = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
        axis = normal./norm(normal);
    else % this bud is too close, loses its chance
        relativePos = 0;
        newPos = 0;
        axis = 0;
    end
end

% Simplified getDaughterPos function - uses only measured angles
function [newPos, relativePos, axis] = getDaughterPos_measured_angles(CELL, DISTANCE_THRESH, cell_angles_radians, force_selection)
    % Surface definition for possible mother cell:
    [Sm, Rm, Tm] = GET_SURFACE_MATRICES(CELL, .1);
    
    % Sample angle from measured distribution
    th = cell_angles_radians(randi(length(cell_angles_radians)));
    ph = 2 * pi * rand;  % Random azimuthal angle for 3D positioning
    
    % Calculate position on ellipsoid surface
    r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2 ).^(-1));
    x = r * cos(th);
    y = r * sin(th) * cos(ph);
    z = r * sin(th) * sin(ph);
    
    % Check distance to existing buds only if force_selection is 0
    if force_selection == 1
        ix_too_close = [];  % Force selection - ignore distance checks
    else
        % Normal behavior - check if new location is too close to existing scars
        if ~isempty(CELL.BudXYZ)
            t = [x; y; z];
            s = CELL.BudXYZ - t;
            d = sqrt(s(1,:).^2 + s(2,:).^2 + s(3,:).^2);
            ix_too_close = find(d < DISTANCE_THRESH);
        else
            ix_too_close = [];  % No existing buds, so no distance check needed
        end
    end
    
    if isempty(ix_too_close) % new bud is successful
        relativePos = [x; y; z];
        newPos = Tm * Rm * [relativePos; 1];
        normal = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
        axis = normal./norm(normal);
    else % this bud is too close, loses its chance
        relativePos = 0;
        newPos = 0;
        axis = 0;
    end
end





function [total_overlap] = return_CHECK_OVERLAPS2(cell_list)

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
    total_overlap=U;
    % flag = U > overlap_thresh;
    % fprintf([num2str(t),',',num2str(U), '\n']);
    
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
