% Thomas C. Day
% CELL DIVISION SIMULATION
% Simulates populations of dividing cells and exports each simulation's
% raw division event times/generations to CSV for downstream analysis
% (binning, ACF-based synchrony index, statistical tests) -- see
% cell_division_freq_by_time_18sep2026.R.
%
% Start with N mother cells at t = 0.
%   Each cell gets a division time sampled from lognormal distribution.
%   When it divides, the individual is "replaced" by 2 daughter cells.
%       (Cell 1 = mother, Cell2 = true daughter), each given its own new
%       division time (drawn from a lognormal).
%   The new division time may experience a "delay" if its the first
%   daughter.
%
% Modifications:
% Made a change in the function "simulateCellDivisions" to have  the 
% standard deviation of the first cell division to be scaled by the 
% differences in means just like in the synthetic log-normal distributions, 
% this was made by adding the following line at the beginning of the 
% function and using that new variable in the sampling of the doubling 
% times:
%       sim_sigmaFirst=sigmaFirst*(muFirst/muSecond);


clearvars; clear; close all; clc;
rng(1); 

% MAX_EVENTS: with N0=30 starting cells and mu=105 min, the population
% roughly doubles every ~105 min (N(t) ~ N0*2^(t/mu)), and each division
% adds exactly 1 net cell. By t=1000 min that's on the order of
% 30*(2^(1000/105) - 1) =~ 2.2e4 division events for the fastest-dividing
% conditions (baseline sigma, all delay values -- their steady-state mu
% is unaffected). Set MAX_EVENTS comfortably above that so runs are only
% stopped by SIM_RUNTIME, not by the cap. NOTE: the event queue below is
% fully re-sorted every iteration (sortrows), so cost grows faster than
% linearly with MAX_EVENTS -- expect this to run substantially slower
% than the old MAX_EVENTS=1e3 version. Time a single run before launching
% the full NSIMS loop, and lower MAX_EVENTS (or SIM_RUNTIME) if it's too
% slow for your needs.
MAX_EVENTS  = 3e4;
SIM_RUNTIME = 1e3;          % total simulation time, minutes
NSIMS       = 30;

NCELLS_mo   = 15;      % starting number of mother cells
NCELLS_da   = 15;      % starting number of daughter cells
mu          = 105;      % mean doubling time, minutes
si          = 15;      % std of doubling time, minutes

%%
% Simulation 1:
% All cells (mothers and every subsequent daughter) draw their division
% time from the same lognormal distribution.

[mu1, sigma1] = lognParams(mu, si);
[divTimes1, divGens1] = simulateCellDivisions(NCELLS_mo, NCELLS_da, mu1, sigma1, mu1, sigma1, SIM_RUNTIME, MAX_EVENTS);

% --- Export raw division events for Simulation 1 ---
sim1_table = table(divTimes1, divGens1, 'VariableNames', {'divTime','divGen'});
writetable(sim1_table, 'sim1_divisions.csv');


%%
% Simulation 2
% All cells draw their division time from the same lognormal distribution,
% sweep the standard deviation by factors of 2

std_values = [15, 20, 25, 30];

% --- Accumulators for CSV export across all sigma factors & runs ---
sim2_sigmaFactor = [];
sim2_run         = [];
sim2_divTime     = [];
sim2_divGen      = [];

for ii = 1:numel(std_values)
    disp(std_values(ii));
    [mu_ii, sigma_ii] = lognParams(mu, std_values(ii));
    divTimes_ii = cell(NSIMS,1);
    divGens_ii  = cell(NSIMS,1);
    for mm = 1:NSIMS
        [divTimes_ii{mm}, divGens_ii{mm}] = simulateCellDivisions(NCELLS_mo, NCELLS_da, mu_ii, sigma_ii, mu_ii, sigma_ii, SIM_RUNTIME, MAX_EVENTS);

        % Append this run's events to the export accumulators
        n = numel(divTimes_ii{mm});
        sim2_sigmaFactor = [sim2_sigmaFactor; repmat(std_values(ii), n, 1)];
        sim2_run         = [sim2_run;         repmat(mm, n, 1)];
        sim2_divTime     = [sim2_divTime;     divTimes_ii{mm}];
        sim2_divGen      = [sim2_divGen;      divGens_ii{mm}];
    end
end

% --- Export raw division events for Simulation 2 ---
sim2_table = table(sim2_sigmaFactor, sim2_run, sim2_divTime, sim2_divGen, ...
    'VariableNames', {'std_value','run','divTime','divGen'});
writetable(sim2_table, 'sim2_divisions.csv');


%% 
% Simulation 3
% Each cell's first daughter's cell division is delayed by an average value
% of DT, with the same variance as the standard lognormal distribution used
% in simulation 1.
% All subsequent cell division times are drawn from the standard lognormal
% distribution used above.

DTvalues = 105*(0:12.5:50)/100; % Modified to have a range from 0 to 100% delay
[muSecond, sigmaSecond] = lognParams(mu, si); % standard dist, used for all divisions after the first

% --- Accumulators for CSV export across all delays & runs ---
sim3_delay    = [];
sim3_run      = [];
sim3_divTime  = [];
sim3_divGen   = [];

for tt = 1:numel(DTvalues)
    disp(DTvalues(tt));
    DT = DTvalues(tt);
    [muFirst, sigmaFirst] = lognParams(mu + DT, si*((mu + DT)/mu));
    divTimes_tt = cell(NSIMS,1);
    divGens_tt  = cell(NSIMS,1);
    for mm = 1:NSIMS
        [divTimes_tt{mm}, divGens_tt{mm}] = simulateCellDivisions(NCELLS_mo, NCELLS_da, muFirst, sigmaFirst, muSecond, sigmaSecond, SIM_RUNTIME, MAX_EVENTS);

        % Append this run's events to the export accumulators
        n = numel(divTimes_tt{mm});
        sim3_delay   = [sim3_delay;   repmat(DT, n, 1)];
        sim3_run     = [sim3_run;     repmat(mm, n, 1)];
        sim3_divTime = [sim3_divTime; divTimes_tt{mm}];
        sim3_divGen  = [sim3_divGen;  divGens_tt{mm}];
    end
end

% --- Export raw division events for Simulation 3 ---
sim3_table = table(sim3_delay, sim3_run, sim3_divTime, sim3_divGen, ...
    'VariableNames', {'delay','run','divTime','divGen'});
writetable(sim3_table, 'sim3_divisions.csv');


%%
% Simulation 4:
% Grande vs. petite, driven by BOOTSTRAP RESAMPLING of empirically
% measured doubling times (time-lapse microscopy, imaged every 5 min)
% rather than by fitted lognormal parameters. Each CSV has columns
% division_number (0 = a cell's first division since being born, 1 = every
% division after that), delta_t, minutes, and hours; only division_number
% and minutes are used here.

% NOTE: update these to your real file paths/names if they differ.
grandeCSV = '/Users/pipe/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/grande_doubling_time_dist_2024dec2.csv';
petiteCSV = '/Users/pipe/cell_division_timing_shapes_the_morphology_and_size_of_nascent_multicellular_organisms/Data/petite_doubling_time_dist_2023may30.csv';

[grandeFirstPool, grandeLaterPool] = loadDivisionTimePools(grandeCSV);
[petiteFirstPool, petiteLaterPool] = loadDivisionTimePools(petiteCSV);

% --- Accumulators for CSV export across Grande/Petite & runs ---
sim4_condition = {};
sim4_run       = [];
sim4_divTime   = [];
sim4_divGen    = [];

divTimes_g = cell(NSIMS,1);
divTimes_p = cell(NSIMS,1);

for mm = 1:NSIMS
    [divTimes_g{mm}, divGens_g{mm}] = simulateCellDivisionsBootstrap(NCELLS_mo, NCELLS_da, grandeFirstPool, grandeLaterPool, SIM_RUNTIME, MAX_EVENTS);
    [divTimes_p{mm}, divGens_p{mm}] = simulateCellDivisionsBootstrap(NCELLS_mo, NCELLS_da, petiteFirstPool, petiteLaterPool, SIM_RUNTIME, MAX_EVENTS);

    % Append this run's events to the export accumulators (Grande)
    n_g = numel(divTimes_g{mm});
    sim4_condition = [sim4_condition; repmat({'Grande'}, n_g, 1)];
    sim4_run       = [sim4_run;       repmat(mm, n_g, 1)];
    sim4_divTime   = [sim4_divTime;   divTimes_g{mm}];
    sim4_divGen    = [sim4_divGen;    divGens_g{mm}];

    % Append this run's events to the export accumulators (Petite)
    n_p = numel(divTimes_p{mm});
    sim4_condition = [sim4_condition; repmat({'Petite'}, n_p, 1)];
    sim4_run       = [sim4_run;       repmat(mm, n_p, 1)];
    sim4_divTime   = [sim4_divTime;   divTimes_p{mm}];
    sim4_divGen    = [sim4_divGen;    divGens_p{mm}];
end

% --- Export raw division events for Simulation 4 ---
sim4_table = table(sim4_condition, sim4_run, sim4_divTime, sim4_divGen, ...
    'VariableNames', {'condition','run','divTime','divGen'});
writetable(sim4_table, 'sim4_divisions.csv');


%% Functions
function [mu, sigma] = lognParams(m, s)
    % Convert a desired mean (m) and std (s) of a lognormal random
    % variable into the mu/sigma parameters of the underlying normal
    % distribution, so that draws can be generated as exp(mu + sigma*randn).
    sigma2 = log(1 + (s/m)^2);
    sigma  = sqrt(sigma2);
    mu     = log(m) - sigma2/2;
end

function [firstPool, laterPool] = loadDivisionTimePools(csvPath)
    % Reads an empirical doubling-time CSV with columns division_number,
    % delta_t, minutes, hours (only division_number and minutes are used).
    % division_number == 0 -> a cell's first division since being born
    % division_number == 1 -> every division after that
    % Returns the two pools of measured minutes-to-divide values to
    % bootstrap-sample from in simulateCellDivisionsBootstrap().
    T = readtable(csvPath);

    isFirst = T.division_number == 0;
    isLater = T.division_number == 1;

    firstPool = T.minutes(isFirst & ~isnan(T.minutes));
    laterPool = T.minutes(isLater & ~isnan(T.minutes));

    assert(~isempty(firstPool), ...
        'loadDivisionTimePools: no division_number==0 rows with a valid minutes value found in %s', csvPath);
    assert(~isempty(laterPool), ...
        'loadDivisionTimePools: no division_number==1 rows with a valid minutes value found in %s', csvPath);
end

function [divTimes, divGens] = simulateCellDivisionsBootstrap(nStart_mo, nStart_da, firstPool, laterPool, tMax, maxEvents)
    % Same event-driven queue logic as simulateCellDivisions() below, but
    % every new division time is drawn by bootstrap resampling (sampling
    % WITH REPLACEMENT) from an empirically measured vector of doubling
    % times, rather than from a fitted lognormal distribution.
    %   firstPool - measured doubling times for a cell's FIRST division
    %               since being born (division_number == 0)
    %   laterPool - measured doubling times for every division after that
    %               (division_number == 1)
    % Returns divTimes/divGens with the same meaning as
    % simulateCellDivisions(): divGens==1 marks a first division, 2+ marks
    % later divisions, matching the mother/daughter convention used
    % throughout this script.

    t0_mo = laterPool(randi(numel(laterPool), nStart_mo, 1));
    t0_da = firstPool(randi(numel(firstPool), nStart_da, 1));
    t0 = [t0_mo; t0_da];
    g0 = [2*ones(nStart_mo,1); ones(nStart_da,1)];
    q0 = [t0, g0];
    queue = sortrows(q0, 1);

    divTimes = zeros(0, 1);
    divGens  = zeros(0, 1);
    nEvents  = 0;
    while ~isempty(queue)

        t = queue(1, 1);
        g = queue(1, 2);

        if t > tMax
            break;
        end

        queue(1, :) = [];

        divTimes(end + 1, 1) = t;
        divGens(end + 1, 1)  = g;
        nEvents = nEvents + 1;
        if nEvents > maxEvents
            warning('simulateCellDivisionsBootstrap:maxEvents', ...
                'Max event cap (%d) reached; stopping early.', maxEvents);
            break;
        end

        % Mother cell's next division: bootstrap-sampled from laterPool
        newT = t + laterPool(randi(numel(laterPool)));
        newG = g + 1;
        queue(end+1,:) = [newT, newG];

        % New daughter cell's first division: bootstrap-sampled from firstPool
        newT = t + firstPool(randi(numel(firstPool)));
        newG = 1;
        queue(end+1,:) = [newT, newG];

        queue = sortrows(queue, 1);
    end
end

function [divTimes, divGens] = simulateCellDivisions(nStart_mo, nStart_da, muFirst, sigmaFirst, muSecond, sigmaSecond, tMax, maxEvents)
    % Event-driven simulation of a doubling cell population.
    %
    % nStart mother cells and daughter cells begin at t = 0. Each cell's FIRST division time
    % is drawn using (muFirst, sigmaFirst). Every division after that,
    % for every descendant, uses (muRest, sigmaRest): new division time =
    % parent's division time + a fresh lognormal draw.
    %
    % Returns:
    %   divTimes - column vector, the time of every division event that
    %              occurred at or before tMax (one entry per division).
    %   divGens  - column vector (same length as divTimes), the
    %              "generation" of each division: 1 = a mother cell's
    %              first division, 2 = a daughter's division, etc.
    %              This is how each cell's division history is tracked.

    % This modification is made for it to match all the other simulations
    % with synthetic log-normal distributions
    % sim_sigmaFirst=sigmaFirst*(muFirst/muSecond);

    % queue columns: [scheduled division time, generation of that division]
    t0_mo = exp(muSecond + sigmaSecond*randn(nStart_mo,1));
    t0_da = exp(muFirst + sigmaFirst * randn(nStart_da, 1));
    t0 = [t0_mo; t0_da];
    g0 = [2*ones(nStart_mo,1); ones(nStart_da,1)];
    q0 = [t0, g0];
    queue = sortrows(q0, 1);

    divTimes = zeros(0, 1);
    divGens  = zeros(0, 1);
    nEvents   = 0;
    while ~isempty(queue)

        t = queue(1, 1);
        g = queue(1, 2);

        if t > tMax
            break;
        end

        % Remove the mother cell:
        queue(1, :) = [];      % remove this event

        % Record the time and generation # of this cell division event
        divTimes(end + 1, 1) = t;
        divGens(end + 1, 1)  = g;
        nEvents = nEvents + 1;
        if nEvents > maxEvents
            warning('simulateCellDivisions:maxEvents', ...
                'Max event cap (%d) reached; stopping early.', maxEvents);
            break;
        end

        % Append the mother cell to the end of the list with new next
        % division time:
        newT = t + exp(muSecond + sigmaSecond*randn(1,1));
        newG = g+1;
        queue(end+1,:) = [newT, newG];

        % Append the daughter cell to the end of the list with the new next
        % division time:
        newT = t + exp(muFirst + sigmaFirst*randn(1,1)); % this cell is always a brand-spanking new cell.
        newG = 1;
        queue(end+1,:) = [newT, newG];
        
        % Sort the queue:
        queue = sortrows(queue, 1);
    end
end