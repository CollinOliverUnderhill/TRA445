%% Battery Pack Simulation Using a 1RC Model
% This script simulates a battery pack comprised of several modules (PCMs) 
% connected in series, where each module consists of multiple cells in parallel.
% The simulation uses a 1RC equivalent circuit model to update each cell's 
% state-of-charge (SOC) and current.
%
% Features:
%   - Random initialization of each cell's SOC, capacity, and internal resistance
%     to simulate cell-to-cell variability.
%   - Optional fault injection (open-circuit and short-circuit) for testing pack behavior.
%   - Data logging and plotting of individual cell SOC and current over time.
%
% Note: Ensure that the file 'soc_ocv.mat' (with two columns: [SOC, OCV]) is in your working directory.

% Some part of this code is from Plett, Gregory L., "Battery Management Systems, Volume II, Equivalent-Circuit Methods," 
% This code is for educational purposes only.

% Author: Xiaolei Bian
% Date: 2025-02-13

%%% Clear Workspace and Initialize
clear; close all; clc;

%%% Simulation and Battery Model Parameters
Ns = 4;       % Number of series-connected modules (PCMs)
Np = 3;       % Number of parallel cells per module

% Load the SOC-OCV lookup data (matrix with columns: [SOC, OCV])
load soc_ocv  

% Simulation time parameters
maxtime = 3600;   % Total simulation time in seconds
restTime = 3000;   % Time (in seconds) after which the pack rests (current set to zero)

% Pre-allocate arrays to store simulation data
storeSOC = zeros(maxtime, Ns, Np);        % To record cell SOC over time
storeCurrent = zeros(maxtime, Ns, Np);    % To record cell current over time

%%% Initialize Cell States and Model Parameters
% Initial state-of-charge (SOC) for cells (default 25%, later randomized)
initialSOC = 0.25;
z = initialSOC * ones(Ns, Np);  % SOC matrix for each cell
ir1 = zeros(Ns, Np);                  % Initial current through the R1 (RC branch)

% Nominal cell capacity (Ah)
Q_nominal = 5;
q = Q_nominal * ones(Ns, Np);  % Cell capacity for each cell

% Parameters for the R1-C1 branch (RC network)
tau = 2.41;                                      % Time constant (seconds)
rc = exp(-1/tau) * ones(Ns, Np);     % Exponential decay coefficient for discrete update

% Resistances (in Ohms)
r1 = 0.0025;                            % R1 resistance
r0 = 0.0111 * ones(Ns, Np);    % Internal ohmic resistance

% Tab resistance for each cell (each cell has 2 tabs)
rt = 0.000125;

%----------------------------------------------------------------------------------------------------------------
%%% Introduce Cell Variability (Randomization)
%%% Randomly set each cell's initial SOC between 30% and 70% (one example)
 z = 0.1 + 0.8 * rand(Ns, Np);
% 
%%% Randomly assign each cell's capacity between 4.5 and 5.5 Ah
 q = 4.5 + rand(Ns, Np);
% 
% % Randomly assign each cell's internal ohmic resistance between 0.005 and 0.025 Ohm,
% % then add the tab resistance (2*rt) to each cell's r0.
r0 = 0.005 + 0.020 * rand(Ns, Np);
r0 = r0 + 2 * rt;  % Add the tab resistance
%----------------------------------------------------------------------------------------------------------------


% %----------------------------------------------------------------------------------------------------------------
% %%% Inject Faults
% %%% To simulate an open-circuit fault, set the corresponding cell's r0 to Inf:
% r0(1,1) = Inf;  % Example: Open-circuit fault in cell 1 of module 1

% %%% To simulate a short-circuit fault, set the corresponding cell's SOC to NaN:
% z(1,1) = NaN;   % Example: Short-circuit fault in cell 1 of module 1
%----------------------------------------------------------------------------------------------------------------


%%% Determine Pack Capacity and Set Initial Current
% The pack capacity is defined by the module with the minimum total capacity 
% (i.e., the sum of capacities of its parallel cells).
totalCap = min(sum(q, 2));    % Pack capacity in Ah
I = 8 * totalCap;                     % Set initial current at 8C rate (in A)

%%% Main Simulation Loop
for k = 1:maxtime
    % 1. Compute Open-Circuit Voltage (OCV) for Each Cell
    % Flatten the SOC matrix to a vector for interpolation
    socVector = z(:);
    deadCells = isnan(socVector);  % Identify cells that are short-circuited
    
    % Interpolate OCV values using the SOC-OCV lookup table (linear interpolation)
    ocvVector = interp1(soc_ocv(:,1), soc_ocv(:,2), socVector, 'linear', 'extrap');
    ocvVector(deadCells) = 0;  % Set OCV to 0 for dead cells
    
    % Reshape the OCV vector back into the original matrix form
    ocv = reshape(ocvVector, size(z));
    
    % 2. Calculate Stable Voltage and Cell Currents
    % Stable voltage is defined as the sum of OCV and the voltage drop across the R1 branch
    % because they are not change immediately when the applied current is changed
    vStable = ocv + r1 .* ir1;
    
    % For short-circuited cells, reset r0 to 2*rt
    r0(isnan(z)) = 2 * rt;
    

    % Compute the module voltage as a weighted average of cell voltages:
    moduleVoltage = (sum(vStable ./ r0, 2) + I) ./ sum(1 ./ r0, 2);
    
    % Calculate the current for each cell:
    iCell = (repmat(moduleVoltage, 1, Np) - vStable) ./ r0;
    
    % 3. Update SOC and R1 Branch Current
    % Update the SOC: SOC_new = SOC_old + (1/3600)*(iCell./q)
    % (1/3600 converts from Ah to the change per second since 1 hour = 3600 seconds)
    z = z + (1/3600) * (iCell ./ q);
    
    % Mark cells with SOC below 0 as dead (short-circuited)
    z(z < 0) = NaN;
    
    % Update the R1 branch current using an exponential update rule
    ir1 = rc .* ir1 + (1 - rc) .* iCell;
    
    % 4. Adjust Charging/Discharging Current Direction
    % If any cell's SOC falls below 5%, ensure the current is set for charging (positive)
    if min(z(:)) < 0.05
        I = abs(I);
    end
    % If any cell's SOC exceeds 95%, set the current for discharging (negative)
    if max(z(:)) > 0.95
        I = -abs(I);
    end
    
    % After the specified restTime, set the current to zero (simulate a rest period)
    if k > restTime
        I = 0;
    end
    
    % 5. Store Simulation Data
    storeSOC(k, :, :) = z;
    storeCurrent(k, :, :) = iCell;
end

%%% Plotting Results
% Create a time vector in minutes
timeMin = (0:maxtime-1) / 60;

% Identify modules without open-circuit faults (modules where the total r0 is not Inf)
nonOpenModules = ~isinf(sum(r0, 2));


% Figure 1: Plot Cell Currents Over Time for Each Module
figure(1); clf;
numRows = Ns;   % The module number, number of row subplot
numCols = 1;
for moduleIdx = 1:Ns
    currentModule = squeeze(storeCurrent(:, moduleIdx, :));
    
    subplot(numRows, numCols, moduleIdx);
    plot(timeMin, currentModule);
    axis([0 ceil(maxtime/60) -101 101]);
    title(sprintf('Cells in PCM %d', moduleIdx));
    xlabel('Time (min)');
    ylabel('Current (A)');
    grid on;
end


% Figure 2: Plot Cell SOC Over Time for Each Module
figure(2); clf;

avgSOC = [];  % To store the average SOC for each module

for moduleIdx = 1:Ns
    % Extract and convert the SOC data for the current module to percentages
    socModule = squeeze(100 * storeSOC(:, moduleIdx, :));
    
    subplot(numRows, numCols, moduleIdx);
    plot(timeMin, socModule);
    axis([0 ceil(maxtime/60) 0 100]);
    title(sprintf('Cells in PCM %d', moduleIdx));
    xlabel('Time (min)');
    ylabel('SOC (%)');
    grid on;
    
    % Replace NaN values (dead cells) with 0 for calculating the average
    socModule(isnan(socModule)) = 0;
    if nonOpenModules(moduleIdx)
        avgSOC = [avgSOC; mean(socModule, 2)'];  % Average SOC of cells for each module 
    end
end


% Figure 3: Plot Average SOC Over Time for Each Module
figure(3); clf;
plot(timeMin, avgSOC');
grid on;
xlabel('Time (min)');
ylabel('SOC (%)');
title('Average SOC for each PCM');
legendStrings = arrayfun(@(x) sprintf('PCM %d', x), find(nonOpenModules), 'UniformOutput', false);
legend(legendStrings, 'Location', 'Best');

% capacity wasted especially for S4P3
Waste_capacity = (15-totalCap)/15;
disp(['The waste of capacity is: ',num2str(Waste_capacity*100),'%'])
% %% Figure 4: Plot Difference Between Maximum and Minimum Average SOC
% figure(4); clf;
% plot(timeMin, max(avgSOC) - min(avgSOC));
% grid on;
% xlabel('Time (min)');
% ylabel('Difference in SOC (%)');
% title('Max-average SOC minus Min-average SOC');



