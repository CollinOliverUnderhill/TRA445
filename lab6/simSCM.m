%% Battery Pack Simulation Using a 1RC Model
% This script simulates a battery pack where cells are connected in series
% to form modules and the modules are connected in parallel to form the pack.
% Each cell is modeled using a 1RC equivalent circuit. Cell-to-cell variability
% is introduced by randomizing the initial state-of-charge (SOC), capacity, 
% and internal resistance.
%
% Note: Ensure that the file 'soc_ocv.mat' (with two columns: [SOC, OCV]) 
% is in your working directory.
%
% Author: Xiaolei Bian (Improved Version)
% Date: 2025-02-13

clear; close all; clc;

%%% Simulation & Model Parameters

% Pack configuration: 
numCells = 6;      % Number of cells in series per module
numModules = 3;    % Number of modules in parallel

% Load the SOC-OCV lookup table (two columns: [SOC, OCV])
load soc_ocv

% Simulation time parameters
totalTime   = 3600;  % Total simulation time in seconds
restTime    = 3000;  % Time (in seconds) after which the pack rests (I = 0)

% Pre-allocate storage for simulation data
storeSOC     = zeros(totalTime, numCells, numModules);      % SOC history
storeCurrent = zeros(totalTime, numCells, numModules);      % Cell current history

%%% Initialize Cell States and Model Parameters

% Initial state-of-charge (will be randomized shortly)
initialSOC = 0.25;
SOC = initialSOC * ones(numCells, numModules);  

% Initial RC branch current (for the R1-C1 network)
i_R1 = zeros(numCells, numModules); 

% Nominal cell capacity (Ah) and randomize each cell's capacity between 4.5 and 5.5 Ah
Q_nominal = 5;
%cellCapacity = 4.5 + rand(numCells, numModules);
cellCapacity = 5;
% Time constant and update coefficient for the RC branch
tau = 2.41;  
rc_update = exp(-1/tau) * ones(numCells, numModules);

% Resistances (Ohms)
R1 = 0.0025;   % Resistance in the RC branch
% Initially, internal (ohmic) resistance is set but later randomized:
R0 = 0.0111 * ones(numCells, numModules);

% Tab resistance (each cell has 2 tabs)
R_tab = 0.000125;

% Introduce Cell Variability

% Randomize the initial SOC between 30% and 70%
SOC = 0.30 + 0.40 * rand(numCells, numModules);

% Randomize each cell's internal resistance between 0.005 and 0.025 Ohm,
% then add the tab resistance (2*R_tab) to each cell.
R0 = 0.005 + 0.020 * rand(numCells, numModules);
R0 = R0 + 2 * R_tab;

%%% (Optional) Inject Faults
% Uncomment one of the following lines to simulate faults:
%   Open-circuit fault in cell (1,1): set its R0 to Inf
% R0(1,1) = Inf;
%
%   Short-circuit fault in cell (1,1): mark its SOC as NaN
% SOC(1,1) = NaN;

%%% Determine Pack Capacity and Set Initial Current

% For a series-connected module, the cell with the lowest capacity limits the module.
% Since modules are in parallel, the overall pack capacity is determined by the module
% with the smallest effective capacity.
%
% Here, we use the sum of cell capacities per module (i.e., along the series direction)
% and select the minimum. (Note: In a real series string, capacity is limited by the
% lowest cell capacity. This approach is one way to capture module imbalance.)
q = 4.5 + rand(numCells, numModules);
moduleCapacities = sum(q, 1);  % 1 x numModules
packCapacity = min(moduleCapacities);     % Limiting module capacity (Ah)

% Set initial applied current (A) at an 8C rate based on the limiting capacity.
I_app = 8 * packCapacity;

%%% Main Simulation Loop

for t = 1:totalTime
    % 1. Compute Open-Circuit Voltage (OCV) for Each Cell
    % Flatten the SOC matrix for interpolation
    SOC_vector = SOC(:);
    deadCells  = isnan(SOC_vector);  % Identify cells that are short-circuited
    
    % Interpolate OCV using the SOC-OCV lookup table (linear interpolation)
    OCV_vector = interp1(soc_ocv(:,1), soc_ocv(:,2), SOC_vector, 'linear', 'extrap');
    OCV_vector(deadCells) = 0;  % For dead cells, set OCV to 0
    
    % Reshape to original dimensions
    OCV = reshape(OCV_vector, size(SOC));
    
    % 2. Compute the Stable Voltage of Each Cell
    % The stable voltage includes the OCV and the voltage across the R1 branch.
    v_stable = OCV + R1 .* i_R1;
    
    % For cells marked as short-circuited (SOC is NaN), force R0 to 2*R_tab
    R0(isnan(SOC)) = 2 * R_tab;
    
    % 3. Compute Module (Pack) Voltage and Cell Currents
    % For each module (series string), calculate:
    %   v_module  = sum of stable voltages of cells in series (1 x numModules)
    %   R_module  = sum of internal resistances in series (1 x numModules)
    v_module = sum(v_stable, 1);
    R_module = sum(R0, 1);
    
    % In a parallel configuration, the pack (node) voltage satisfies:
    %   I_app = sum_{m=1}^{numModules} (V_pack - v_module(m)) / R_module(m)
    % Solve for the pack voltage:
    V_pack = (I_app + sum(v_module ./ R_module)) / sum(1 ./ R_module);
    
    % Compute the current for each module:
    %   i_module = (V_pack - v_module) ./ R_module   (1 x numModules)
    i_module = (V_pack - v_module) ./ R_module;
    
    % Since cells in series carry the same current, replicate the module current
    % for each cell in the module.
    i_cell = repmat(i_module, numCells, 1);
    
    % 4. Update SOC and RC Branch Current for Each Cell
    % SOC update (convert A to Ah per second: 1/3600 factor)
    SOC = SOC + (1/3600) .* (i_cell ./ cellCapacity);
    
    % Mark cells with SOC falling below 0 as dead (short-circuited)
    SOC(SOC < 0) = NaN;
    
    % Update RC branch current using an exponential update rule
    i_R1 = rc_update .* i_R1 + (1 - rc_update) .* i_cell;
    
    % 5. Adjust Charging/Discharging Current Direction
    if min(SOC(:)) < 0.05
        I_app = abs(I_app);   % Force charging if any cell is very low
    end
    if max(SOC(:)) > 0.95
        I_app = -abs(I_app);  % Force discharging if any cell is near full
    end
    
    % After the specified rest time, set applied current to zero (simulate rest)
    if t > restTime
        I_app = 0;
    end
    
    % 6. Store Data for Post-Processing
    storeSOC(t, :, :)     = SOC;
    storeCurrent(t, :, :) = i_cell;
end

%%% Plotting Results

% Create a time vector in minutes for plotting
timeMin = (0:totalTime-1) / 60;

% Identify modules without open-circuit faults.
% (A module is considered open if its total internal resistance is Inf.)
nonOpenModules = ~isinf(R_module);

% Figure 1: Cell Currents Over Time in Each Module
figure(1); clf;
for m = 1:numModules
    % Extract current data for module m (each row corresponds to a cell)
    currentModule = squeeze(storeCurrent(:, :, m));
    
    subplot(numModules, 1, m);
    plot(timeMin, currentModule, 'LineWidth', 1.5);
    axis([0 ceil(totalTime/60) -101 101]);
    title(sprintf('Module %d: Cell Currents', m));
    xlabel('Time (min)');
    ylabel('Current (A)');
    grid on;
end

% Figure 2: Cell SOC Over Time in Each Module
figure(2); clf;
avgSOC = [];  % To store average SOC (in %) for each module

for m = 1:numModules
    % Extract and convert SOC data to percentage
    socModule = squeeze(100 * storeSOC(:, :, m));
    
    subplot(numModules, 1, m);
    plot(timeMin, socModule, 'LineWidth', 1.5);
    axis([0 ceil(totalTime/60) 0 100]);
    title(sprintf('Module %d: Cell SOC (%)', m));
    xlabel('Time (min)');
    ylabel('SOC (%)');
    grid on;
    avgSOC = [avgSOC; mean(socModule, 2, 'omitnan')'];
end

% Figure 3: Average SOC Over Time for Each Module
figure(3); clf;
plot(timeMin, avgSOC', 'LineWidth', 2);
grid on;
xlabel('Time (min)');
ylabel('Average SOC (%)');
title('Average SOC for Each Module');
% Create legend strings only for modules that are not open-circuited
legendStr = arrayfun(@(m) sprintf('Module %d', m), find(nonOpenModules), 'UniformOutput', false);
legend(legendStr, 'Location', 'Best');


%%% Capacity Waste Calculation

% Compute capacity waste percentage
capacityWaste = (Q_nominal*numCells - packCapacity) /Q_nominal/numCells * 100;

% Display result
fprintf('Capacity Waste: %.5f%%\n', capacityWaste);

