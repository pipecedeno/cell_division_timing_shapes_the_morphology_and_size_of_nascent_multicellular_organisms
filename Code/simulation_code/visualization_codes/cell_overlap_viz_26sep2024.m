 % Thomas C. Day, 2022 edit
% This script generates N clusters of snowflake yeast, modeled as connected
% cells, where each cell is a prolate ellipsoid. The stop condition of the
% experiment can either be a number of generations or a total overlap stop
% condition. Cell aspect ratio can also be varied. Last, there are options
% to include force relaxation to cell positions and orientations.

%{
Luis Felipe 2023
Note: Most of my comments use %# to differentiate them from Tom's comments
this code is based on snowflake_generation_TCD_2023_09_05.m
This version of the code works for calculating cluster sizes using the edges loaded from python simulations
the part of the include_forces is broken for this code as the ID variable was modified to use the network 
ids and those ids are used in other parts of the code

Nov 16 2023
The code now also includes the option to calculate the individual cell overlaps (it still doesn't include the
option to save them into a file), it can make the 3D plot of the true overlaps between cells and it can save
a rotating video of that 3D plot

Input variables added:
    Needed:
        -file_path: path to a file of the list of edges created in the network simulation

    Flags:
        -write_to_output_file: this file is to save the cluster size at fracture and how many edges where added or not
            If the file is going to be saved the you can change the name of the output file by modifying output_filename
        -calc_ind_overlap: Should each cell overlap be calculated? (They aren't saved for the moment)
        -flag_overlap_plot: Should the plot of the cells overlap be made? (If several clusters are being simulated only
            the last one is shown)
        -flag_save_overlap_video: Should a rotating video of the overlap plot be made?
            the video file name is 'rotating_overlaps.avi' by default
            The angle of the rotation and more parameters could be changed in save_video_3d_plot
            If you saved a video you need to clear the workspace before saving another video (sometimes it still causes an error
            even when doing that, so you could clear the workspace again or just restart matlab)
            The video is saved of the last cluster, or after force is calculated
        -save_ind_overlaps: Should individual overlaps be saved?
            change overlaps_filename value to change the name of the saved file


Note: I modified the functions of get_overlap and its modified versions so that it doesn't mark the error
of not enough values to use convhull, so now I can try to modify include_forces to see if it works

Nov 18
Now include forces works
I had to fix some parameters that I moved from ELYES_SIM when creating the daughters and how some
variables were initialized. And the I just modified get_forces_torques to its original state as
before I made some modifications on it whe trying to fix the code

Nov 30 
Adding option to save edges that where not added

%}

%% Inputs:
folder = pwd;
cd(folder);


% Evolved test file
% file_path = '/Users/pipe/work_dir/observed_synchrony/evolution_results/physics_sim_26sep2024/clust_size_sync_vs_async_26sep2024/test_pa2-t200_1000m_50n/edges_sim_files/edges_network_sim_2.csv';
% Petite test file
file_path = '/Users/pipe/work_dir/observed_synchrony/evolution_results/physics_sim_26sep2024/clust_size_sync_vs_async_26sep2024/test_petite_1000m_50n/edges_sim_files/edges_network_sim_3.csv';
% Grande test file
% file_path = '/Users/pipe/work_dir/observed_synchrony/evolution_results/physics_sim_26sep2024/clust_size_sync_vs_async_26sep2024/test_grande_1000m_50n/edges_sim_files/edges_network_sim_2.csv';

%# Reading edges file and making the vectors of the node1 and node2 that form the edge added
edges_data = readtable(file_path, 'Delimiter', ',', 'ReadVariableNames', true);
node1_vals      = edges_data.node1;
node2_vals      = edges_data.node2;


%# Show overlap plot?
flag_overlap_plot=1; % 1 = Y, 0 = N
%# Save rotating video of the 3d plot of the overlaps
flag_save_overlap_video=0; % 1 = Y, 0 = N

%# This output file is for the size at fracture of each of the simulations and the edges
% that were not added
output_filename = 'test_output_file.csv';
write_to_output_file=0;

%# Calculate individual cells overlap? This is for the main loop
calc_ind_overlap=0; % 1 = Y, 0 = N
%# save individual overlaps to file?
save_ind_overlaps=0; % 1 = Y, 0 = N
overlaps_filename = 'test_overlaps.csv';

clust_size_output=0; %#initializing variable to avoid errors even if it is not used
if write_to_output_file==1
    clust_size_output = create_file(output_filename);

    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if clust_size_output == -1
        return;
    end

    % Add file header
    fprintf(clust_size_output, 'sim_number,size_fracture,failed_to_add,missing_mother_edge\n');
    %# Note: the total amount of edges not included is the sum of failed_to_add and missing_mother_edge,
    %I just counted separate in case it was necessary but I am thinking that the sum is the important value
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
    fprintf(overlaps_output, 'file_num,sim_number,node_id,degree,overlap_vol,num_small_overlaps,num_big_overlaps\n');
end

N               = 1; % number of clusters to generate
numGens         = 5; % number of generations of cell division in a group
diam            = 4.58; % smallest diameter of the cell
err_diam        = 0; % variation in cell size, taken from coefficient of variation data
AR              = 1.2;%:.1:2.8; % aspect ratio value
err_AR          = 0; % standard deviation in aspect ratio from Shane's measurements
pole_theta      = deg2rad(10); % buds nearest the pole will be chosen from between 0 and 10 degrees in polar angle
THETA           = deg2rad(45); % polar angle average from SEM data
thetaVariance   = deg2rad(0); % variation in polar angle from SEM data
distance_thresh = 1.1672; % minimum distance (um) separating bud scars
new_bud_prob    = .8; % probability that the first cell will bud near the pole
check_overlap   = 1; % do we check the overlap?
overlap_thresh  = 3.0e1; % threshold of total overlaps %#2.5e2


% Grow the cluster to the same amount of cells
check_number_of_cells=0; % if 0 then it is growing the cluster until fragmentation
grow_to_cells=200;


% Force input parameters:
include_forces  = 0; % do we include forces for rearrangements? 1 = Y, 0 = N
NEIGHBOR_THRESH = 3*diam; % radius at which you consider cells "neighbors"
STERIC_MAG      = 2; % relative strength of steric interactions(from cells colliding)
CHITIN_MAG      = 1; % relative strength of the bonds holding cells together
BOND_TORQ_MAG   = 1; % relative torsional strength of chitin bonds
mobility_pos    = .1;
mobility_rot    = .1;
dt              = .1; % timestep for mechanical relaxation
T               = 10; % number of timesteps to allow mechanical relaxation %# 25

% Do we visualize figures?
figure_viz = 0; % 1 = Y, 0 = N



save_edges_not_added=0; % 1 = Y, 0 = N
edges_filename = "edges_not_added.csv";
% clust_edges_filename='clust_size_edges_not_added.csv';

%# Creating file for saving edges not added
edges_output=0; %#initializing variable to avoid errors even if it is not used
edges_clust_size_output=0;
if save_edges_not_added==1
    edges_output = create_file(edges_filename);
    % Check if the file already existed, if it already existed the program
    % execution is going to be stopped
    if edges_output == -1
        return;
    end
    % Add file header
    fprintf(edges_output, 'sim_number,node1,node2,degree_cell,failed_to_add,missing_mother_edge\n');

    % edges_clust_size_output=create_file(clust_edges_filename);
    % if edges_clust_size_output == -1
    %     return;
    % end

    % fprintf(edges_clust_size_output, 'sim_number,fracture_size\n')

end

% Creating input variables dictionary
global input_variables;
input_variables=containers.Map();
input_variables('save_edges_not_added')=save_edges_not_added;
input_variables('edges_output')=edges_output;
input_variables('grow_to_cells')=grow_to_cells;
input_variables('check_number_of_cells')=check_number_of_cells;

% Initialization:
fprintf('Welcome to the cell simulator v2.0\n');
fprintf('Written by Thomas C. Day, 2020\n');

%% MAIN CODE: -------------------------------------------------------------
% -------------------------------------------------------------------------
% Parallel iteration over varied params:
% save('params.mat','N','numGens','diam','pole_theta','THETA','distance_thresh','new_bud_prob','overlap_thresh');
Ncells = zeros(length(AR),1);
Scells = Ncells;
for a = 1:length(AR)
    cells_sim = cell(N,1);
    aspRat    = AR(a);
    theta     = THETA;
    filename = ['sim','_AspRat=',num2str(aspRat)];
    filename = strrep(filename, '.','-');
    COM = zeros(3,N);
    RAD = zeros(1,N);
	ncells = zeros(1,N);
	
    for ii = 1:N
        
        % Generate one cluster:
        fprintf(['Cluster #: ',num2str(ii), ' / ', num2str(N), '\n']);
        [cell_list] = ELYES_SIM(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, numGens, NEIGHBOR_THRESH, check_overlap, overlap_thresh, node1_vals, node2_vals, clust_size_output, write_to_output_file, ii);
        ncells(ii) = length(cell_list);
		

        % Allow for force relaxations from steric interactions:
        if include_forces == 1
            warning('off');
            [~, cell_list, overlaps] = INCLUDE_FORCES(cell_list, NEIGHBOR_THRESH, STERIC_MAG, CHITIN_MAG, BOND_TORQ_MAG, mobility_pos, mobility_rot, dt, T, figure_viz, overlap_thresh, NEIGHBOR_THRESH);
            warning('on');
        end        
 
        % Calculating individual cell overlaps
        if calc_ind_overlap == 1
            % cell_overlaps = calculate_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);
            [cell_overlaps, big_overlaps_list, small_overlaps_list] = calculate_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);

            if save_ind_overlaps==1
                    save_ind_overlaps_to_file(overlaps_output, 1, ii, cell_list, cell_overlaps, big_overlaps_list, small_overlaps_list);
            end
        end

		
        
        % if save_edges_not_added==1
        %     fprintf(edges_clust_size_output, '%d,%d\n', ii, length(cell_list));
		% end

    end
    fprintf(['Mean number of cells per group: ',num2str(mean(ncells)),'\n']);
    fprintf(['Stdev number of cells per group: ',num2str(std(ncells)),'\n']);
    Ncells(a) = mean(ncells);
    Scells(a) = std(ncells);

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

    %# This is commented out because only one video can be made per code execution
    % save_video_3d_plot('rotating_normal_cluster.avi');

end

%# Function to make the plot of the individual overlaps
if flag_overlap_plot == 1
    cell_overlaps = plot_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);

    %# This part saves the rotating video of the 3d plot
    if flag_save_overlap_video==1
        save_video_3d_plot('rotating_overlaps.avi'); %# change this name if you want to save the movie to a different file
    end
end


%# First attempt plot to show the overlap volume
% Here the cells are shaded according to how high their overlap volume is, and
%it plots all cells with a black shaded area
% if flag_overlap_plot == 1
%     figure;
%     hold on; box on; set(gca,'linewidth',2);

%     % Normalize cell_overlaps values to range [0, 1]
%     norm_overlaps = cell_overlaps / max(cell_overlaps);

%     % Choose a constant color for all cells (e.g., blue)
%     constantColor = [0, 0, 1]; % This is an RGB triplet for blue

%     % Loop over each cell in the list
%     for k = 1:length(cell_list)
%         [x,y,z] = VISUALIZE_ELLIPSOID(cell_list(k), 30);

%         % Use normalized overlap value as transparency
%         alphaValue = norm_overlaps(k);

%         % Draw the ellipsoid with constant color and varying transparency
%         surf(x,y,z,'facecolor',constantColor,'edgecolor','none','facealpha',alphaValue);
%     end

%     constantColor = [0, 0, 0]; % This is an RGB triplet for blue

%     % Loop over each cell in the list
%     for k = 1:length(cell_list)
%         [x,y,z] = VISUALIZE_ELLIPSOID(cell_list(k), 30);

%         % Use normalized overlap value as transparency
%         alphaValue = 0.05;

%         % Draw the ellipsoid with constant color and varying transparency
%         surf(x,y,z,'facecolor',constantColor,'edgecolor','none','facealpha',alphaValue);
%     end

%     view(3); axis equal;
%     lighting gouraud;
%     lightangle(0,30);
%     material dull;
% end


%# Close output file
if write_to_output_file==1
    fclose(clust_size_output);
end

if save_ind_overlaps==1
    fclose(overlaps_output);
end



%% Added functions by Luis 

%# Function to make the movie of the 3D rotating plot
%This code could be made faster if not all degrees are saved for the video
% as right now it is rotating the image 360 degree one degree at a time
function save_video_3d_plot(output_filename)
    % Set up the video writer
    videoWriter = VideoWriter(output_filename);
    videoWriter.Quality = 95; % Set to desired quality level, from 0 to 100

    % Open the video writer
    open(videoWriter);

    %# this can be uncommented and modifed if you want to save the videos with a different resolution
    % figureHandle = gcf; % Get the handle to the current figure
    % set(figureHandle, 'Position', [100, 100, 850, 582]); % Set the figure size

    % Rotate and capture frames
    for angle = 1:360
        % Rotate the view
        view(angle, 30); % 30 is the elevation angle, you can adjust it

        % Capture the frame
        frame = getframe(gcf);

        % Write frame to video
        writeVideo(videoWriter, frame);
    end

    % Close the video writer
    close(videoWriter);
end

%# Function to create output file and avoid rewriting results
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

%# Function to calculate each cells total overlap volume
% The function is modification of get_overlaps and it only returns the cumulative value
% of that cell
% probably some of the stuff that is no longer used could be deleted to increase performance of
% the function
function [sum_overlap, num_big_overlaps, num_small_overlaps] = get_cell_overlap(cell_list, ix_o_i, nbors)
    
    % Obtain overlapping volume and center point of overlap:
    cell1       = cell_list(ix_o_i);
    Force_pts   = [];
    Overlap_pts = [];
    Overlaps    = [];
    Directions  = [];

    sum_overlap=0;
    num_small_overlaps=0;
    num_big_overlaps=0;

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

                sum_overlap=sum_overlap+overlap(2);
                num_small_overlaps=num_small_overlaps+1;
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

            sum_overlap=sum_overlap+overlap(2);
            num_big_overlaps=num_big_overlaps+1;
        end
        % Store info to cell array:
        Force_pts = [Force_pts, Fpt_1];
        Overlap_pts = [Overlap_pts, COM];
        Overlaps = [Overlaps, overlap];
        Directions = [Directions, direction];

    end

end



function save_ind_overlaps_to_file(overlaps_output, file_number, sim_number, cell_list, cell_overlaps, big_overlaps_list, small_overlaps_list)

    for i = 1:length(cell_list)
        %cell id (network)
        temp_cell_id=cell_list(i).networkID;

        %cell degree
        %The plus one is to take into account their connection to their mother
        temp_degree=length(cell_list(i).Daughters)+1;

        fprintf(overlaps_output, '%d,%d,%d,%d,%d,%d,%d\n', file_number, sim_number, temp_cell_id, temp_degree, cell_overlaps(i), big_overlaps_list(i), small_overlaps_list(i));

    end
end


%# function that coordinates the plot of the overlaps
% colors are coded using RGB vectors, and alpha values for transparency
function [final_overlaps] = plot_ind_cell_overlaps(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude)
    % Calculate the forces and torques acting on all cells in the group.
    % This function can get quite heavy.

    % Obtain distances between all particles:
    [D, ~, ~, ~] = get_particle_distances(cell_list);
    
    % Find all pairs of cells within the neighborhood threshold:
    N = D < neighbor_thresh; % finds all neighbor cells
    N = N - eye(size(N)); % we don't care about counting a cell as its own neighbor
    
    % Obtain forces/torques for each cell:
    % Forces = zeros(3, length(cell_list));
    % Torques = zeros(3, length(cell_list));

    figure; hold on; box on; grid off; set(gca,'linewidth',2);

    constantColor = [1, 0, 0];

    % Loop over each cell in the list
    for k = 1:length(cell_list)
        [x,y,z] = VISUALIZE_ELLIPSOID_scaled(cell_list(k), 30, 0.95);

        alphaValue = 0.05;

        % Draw the ellipsoid
        surf(x,y,z,'facecolor',constantColor,'edgecolor','none','facealpha',alphaValue);
    end

    final_overlaps = [];
    for n = 1:length(cell_list)
        neighbors = find(N(:,n) == 1); % only need to check the neighboring cells for interactions
        
        temp_cell_overlap = get_cell_overlap_plot(cell_list, n, neighbors);
        % final_overlaps = [final_overlaps, temp_cell_overlap]
    end
    %# not sure what all of this stuff do
    view(3); axis equal;
    lighting gouraud;
    lightangle(0,30);
    material dull;
end

%# Function to plot each cell overlap
% function modified from get_cell_overlap and uses INTERSECTION_plot to plot each of the intersections
% Note: each overlap is being plotted twice or more depending of how many cells share that overlap
% but for testing purposes, and because this is only done once, I don't think it matter that much
%# I just noticed that it is only plotting the overlaps when it is "a significant amount"
% whatever that means in Tom words lol
function [sum_overlap] = get_cell_overlap_plot(cell_list, ix_o_i, nbors)
    
    % Obtain overlapping volume and center point of overlap:
    cell1       = cell_list(ix_o_i);
    Force_pts   = [];
    Overlap_pts = [];
    Overlaps    = [];
    Directions  = [];

    sum_overlap=0;

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
                hasMultipleUniqueValues(i) = length(uniqueValues) > 1;
            end

            if all(hasMultipleUniqueValues)
                [~,vol] = convhull(inters(:,1),inters(:,2),inters(:,3));
                overlap = [nbors(j); vol];
                COM     = mean(inters)';
                Fpt_1   = get_force_point(COM, cell1);% Find force point on the surface of the cell of interest
                direc   = M1 * [0;0;0;1];
                direction = direc(1:3) - Fpt_1;
                direction = direction/norm(direction);

                sum_overlap=sum_overlap+overlap(2);
            else
                overlap = []; % add nothing to the overlaps
                COM     = [];
                Fpt_1   = [];
                direction = [];
            end

            sum_overlap=sum_overlap+overlap;
        else % they overlap by a significant amount
            v1       = [x1,y1,z1]; % vertices defining cell1
            v2       = [x2,y2,z2]; % vertices defining cell2
            r1_E2    = M2 \ r1;
            D_help   = vecnorm(r1_E2(1:3,:));
            ix_E2    = find(D_help < 1);
            inters   = [r1(1:3,ix_E2), r2(1:3,ix)]';
            COM      = mean(inters)'; % a point common to both cells
            [inters] = INTERSECTION_plot(v1, v2, COM');
            [~,vol]  = convhull(inters(:,1), inters(:,2), inters(:,3)); % returns the volume of overlap
            overlap  = [nbors(j); vol]; % add to overlaps array
            COM      = mean(inters)'; % center of the overlap volume
            Fpt_1    = get_force_point(COM, cell_list(ix_o_i));
            direc    = M1 * [0;0;0;1];
            direction = direc(1:3) - Fpt_1;
            direction = direction/norm(direction);

            sum_overlap=sum_overlap+overlap;
        end
        % Store info to cell array:
        Force_pts = [Force_pts, Fpt_1];
        Overlap_pts = [Overlap_pts, COM];
        Overlaps = [Overlaps, overlap];
        Directions = [Directions, direction];

        %# For testing purposes to make the code only plot some overlaps
        % if j>3
        %     break
        % end
    end

end


%# Function to plot each intersection
% This is modified from the INTERSECTION function
function [vertices_intersection] = INTERSECTION_plot(vertices1, vertices2, center)
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
    
    % For the plot, everything is commented out except for the line that plots
    % the initialization of the plot is made in plot_ind_cell_overlaps when performing the
    % plot of all overlaps in a cluster
    % figure; hold on; box on; grid on; set(gca,'linewidth',2);
    [k1,~] = convhull(vertices1);
    [k1d,~] = convhull(vdual1);
    [k2,~] = convhull(vertices2);
    [k2d,~] = convhull(vdual2);
    [k3,~] = convhull(vertices_intersection);
    [k3d,~] = convhull(vertices_dual);
    % plot3(center(1), center(2), center(3), '.','markersize',20,'color',[.8,.5,0]);
    % trisurf(k1, vertices1(:,1), vertices1(:,2), vertices1(:,3),'facecolor',[.7,.7,.7],'facealpha',.3,'edgecolor','k','linewidth',3);
    % trisurf(k2, vertices2(:,1), vertices2(:,2), vertices2(:,3),'facecolor',[.7,.1,.7],'facealpha',.2,'edgecolor',[.4,0,.4],'linewidth',3);
    trisurf(k3, vertices_intersection(:,1), vertices_intersection(:,2), vertices_intersection(:,3), 'facecolor',[0,0,1],'facealpha',.2, 'edgecolor', 'none');
    % Parameters I took out of the function: 'edgecolor',[.4,0,0],'linewidth',4

    % view(3);
    % axis equal;
    % title('Real space');
    
    % figure; hold on; box on; grid on; set(gca,'linewidth',2);
    % trisurf(k1d, vdual1(:,1), vdual1(:,2), vdual1(:,3), 'facecolor',[.7,.7,.7],'facealpha',.5,'edgecolor','b','linewidth',3);
    % trisurf(k2d, vdual2(:,1), vdual2(:,2), vdual2(:,3), 'facecolor',[.7,.1,.7],'facealpha',.5,'edgecolor','k','linewidth',2);
    % trisurf(k3d, vertices_dual(:,1), vertices_dual(:,2), vertices_dual(:,3),'facecolor', [.7,0,0],'facealpha',.5,'edgecolor',[.4,0,0],'linewidth',4);
    % view(3);
    % axis equal;
    % title('Dual space');

end

%# Function to visualize all cells of the cluster when doing the overlaps plot
% This is a modification of VISUALIZE_ELLIPSOID because there was a visual bug when using the original
% function because there were too many planes in the same plane, and it wasn't showing a smooth color, so
% to fix the problem I made this function that scales down the size of the ellipses to avoid the overlap between
% colored surfaces
function [xScaled,yScaled,zScaled] = VISUALIZE_ELLIPSOID_scaled(cell_of_interest, resolution, scaleFactor)
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

    centerX = mean(x(:));
    centerY = mean(y(:));
    centerZ = mean(z(:));

    % Scale the ellipse
    xScaled = centerX + scaleFactor * (x - centerX);
    yScaled = centerY + scaleFactor * (y - centerY);
    zScaled = centerZ + scaleFactor * (z - centerZ);

end


%% FUNCTIONS

%# note this function was also modified by Luis
function [cell_list] = ELYES_SIM(diam, err_diam, aspRat, err_AR, pole_theta, theta, thetaVariance, distance_thresh, new_bud_prob, numGens, neighbor_thresh, check_overlap, overlap_thresh, node1_vals, node2_vals, clust_size_output, write_to_output_file, cont_sim)

    global input_variables;
    save_edges_not_added=input_variables('save_edges_not_added');
    edges_output=input_variables('edges_output');
    grow_to_cells=input_variables('grow_to_cells');
    check_number_of_cells=input_variables('check_number_of_cells');

    cell_list       = [];  % this will eventually be the tabulated list
    
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

    % Iterate over generations:
    % for g = 1:numGens
    % Iterate over the edges of the loaded file
    num_rows_edges=length(node1_vals);
    for cont_rows = 1:num_rows_edges
        
        % First check the total amount of overlap in the cluster:
        %fprintf(['Gen = ',num2str(g),'\n']);
        if check_overlap == 1
            flag = CHECK_OVERLAPS2(cell_list, overlap_thresh);
        else
            flag = 0;
        end
        
        if check_number_of_cells==1
            if grow_to_cells<=length(cell_list)
                flag=1;
            end
        end

        % disp(flag);

        if flag ~= 1
            
            % If there is not too much overlap:
            % nPossible   = length(cell_list); % number of possible new cells
            %# What to do with this variables?
            % variedTheta = theta + thetaVariance * randn(nPossible, 1); % list of theta values corresponding to these cells

            % Assign positions for the newest cells:
            % for n = 1:nPossible
                % n is the mother cell

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

                while cont_attempts<10 & newBud == 0 %if one of this conditions is met then it needs to stop
                    %# This may need to be outside the loop, now we are only sampling 1 value because it is processing one cell at a time
                    variedTheta = theta + thetaVariance * randn(1); %theta value corresponding to the cell
                    %# To test if variedtheta gets different values or not when
                    %# thetaVariance is 0
                    % fprintf(['variedTheta: ', num2str(variedTheta), '\n']); %# it prints always the same value

                    [newBud, newBudRel, newAxis] = getDaughterPos(cell_list(n), pole_theta, variedTheta, distance_thresh, new_bud_prob); % finds the new bud xyz position and relative position to the old cell
                    cont_attempts=cont_attempts+1;
                end

                % If there wasn't space for the cell then it is just sampled randomly
                if newBud==0
                    variedTheta = theta + thetaVariance * randn(1);
                    [newBud, newBudRel, newAxis] = getDaughterPos_always(cell_list(n), pole_theta, variedTheta, distance_thresh, new_bud_prob); % finds the new bud xyz position and relative position to the old cell
                end

                if newBud ~= 0 % if the budding chance was successful

                    if save_edges_not_added==1
                        temp_degree=length(cell_list(n).Daughters)+1;

                        fprintf(edges_output, '%d,%d,%d,%d,F,F\n', cont_sim, mother_cell_id, daughter_cell_id, temp_degree);
                    end

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
        fprintf(clust_size_output, '%d,%d,%d,%d\n', cont_sim, length(cell_list), cont_failed_to_add_edge, cont_missing_mother_cell);
    end
    
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


function [newPos, relativePos, axis] = getDaughterPos(CELL, POLE_TH, TH, DISTANCE_THRESH, NEW_BUD_PROB_THRESH)

    % Surface definition for possible mother cell:
    [Sm, Rm, Tm] = GET_SURFACE_MATRICES(CELL, .1);

    % Determine where on cell body to bud the next scar:
    if isempty(CELL.Daughters) % reproduce near the pole with prob 0.7
        frac = rand;
        if frac < NEW_BUD_PROB_THRESH
            % bud at pole
            th = POLE_TH * rand;
            ph = 2 * pi * rand;
            ix_too_close = [];
            r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
            x = r * cos(th);
            y = r * sin(th) * cos(ph);
            z = r * sin(th) * sin(ph);
        else
            % bud on the polar angle
            ph = 2 * pi * rand;
            th = TH;
            ix_too_close = [];
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

    if isempty(ix_too_close) % new bud is successful
        relativePos = [x; y; z];
        newPos      = Tm * Rm * [relativePos; 1];
        normal      = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
        axis        = normal./norm(normal);
    else % this bud is too close, loses its chance
        relativePos = 0;
        newPos      = 0;
        axis        = 0;
    end

end


function [newPos, relativePos, axis] = getDaughterPos_always(CELL, POLE_TH, TH, DISTANCE_THRESH, NEW_BUD_PROB_THRESH)

    % Surface definition for possible mother cell:
    [Sm, Rm, Tm] = GET_SURFACE_MATRICES(CELL, .1);

    % Determine where on cell body to bud the next scar:
    if isempty(CELL.Daughters) % reproduce near the pole with prob 0.7
        frac = rand;
        if frac < NEW_BUD_PROB_THRESH
            % bud at pole
            th = POLE_TH * rand;
            ph = 2 * pi * rand;
            ix_too_close = [];
            r = sqrt(( (cos(th)/Sm(1,1))^2 + (sin(th)*cos(ph)/Sm(2,2))^2 + (sin(th)*sin(ph)/Sm(3,3))^2  ).^(-1));
            x = r * cos(th);
            y = r * sin(th) * cos(ph);
            z = r * sin(th) * sin(ph);
        else
            % bud on the polar angle
            ph = 2 * pi * rand;
            th = TH;
            ix_too_close = [];
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

    relativePos = [x; y; z];
    newPos      = Tm * Rm * [relativePos; 1];
    normal      = 2 * [x/Sm(1,1)^2; y/Sm(2,2)^2; z/Sm(3,3)^2];
    axis        = normal./norm(normal);

end

function [cell_list_in, cell_list_out, Overlaps] = INCLUDE_FORCES(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude, mobility_pos, mobility_rot, dt, T, figure_viz, overlap_thresh, NEIGHBOR_THRESH)
    % This function takes as input a grown snowflake and calculates the forces
    % and torques on each cell. From these forces, it the calculates any
    % rearrangements.

    cell_list_in = cell_list;

    %# To print the overlap as the include forces function is working
    flag=CHECK_OVERLAPS2_test(cell_list, overlap_thresh, 0);

    % cell_overlaps = plot_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);

    for t = 1:T
        % fprintf(['Time = ',num2str(t),'\n']);
        [Forces, Torques, flog1, flog2, Overlaps] = get_forces_torques(cell_list, neighbor_thresh, steric_magnitude, chitin_magnitude, bond_torque_magnitude);
        [cell_list] = UPDATE_POSITIONS(cell_list, Forces, Torques, mobility_pos, mobility_rot, dt);
        
        %# To print the overlap as the include forces function is working
        flag=CHECK_OVERLAPS2_test(cell_list, overlap_thresh, t);
        
        % if mod(t, 5)==0
        %     cell_overlaps = plot_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);
        % end
        % cell_overlaps = plot_ind_cell_overlaps(cell_list, NEIGHBOR_THRESH);

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

function [flag] = CHECK_OVERLAPS2_test(cell_list, overlap_thresh, t)

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
    fprintf([num2str(t),',',num2str(U), '\n']);
    
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