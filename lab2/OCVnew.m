%% SoC estimation using OCV method
clc; clear;

% Determine OCV-SOC Relationship
OCV_SOC_25C = OCV_SOC(); 
% OCV_SOC_25C represents the relationship between OCV and SOC at 25°C, simulating normal temperature conditions.

% Load Test Data and Preprocess
load('DST_80SOC_25C.mat');  % Load test data
t = num(:,2);         % Test time (s)
I = num(:,7);         % Current (A): positive for charge, negative for discharge
V = num(:,8);       % Terminal voltage (V)
Dc = num(:,10);   % Discharge capacity (Ah)

% Clear unnecessary variables
clear num txt;

% Plot the data 


% Apply the OCV method
% Select data from the start index, the data after fully charged
ST = 333;
I = I(ST:end);
V = V(ST:end);
t = t(ST:end) - t(ST); % Adjust time to start from 0
SOC_0 = 100;          % Initial SoC
Q = 2;                       % Cell capacity
R0 = 0.08191;          % Ohmic Resistance, obtained from pulse test 
R0_values = [0.01, 0.07, 0.08191, 0.1, 0.1]; 

% Initialize variables
n = length(I);
SOC = zeros(n, 1);
SOC(1) = SOC_0;

% Coulomb counting
for i = 1:n-1
    dt = t(i+1) - t(i); % Time step
    SOC(i+1) = SOC(i) + (dt * I(i) / (Q * 3600)) * 100;   % measured/reference SoC
end


% SOC-OCV interpolation setup
SOC_points = OCV_SOC_25C(:, 1); % SOC sample points
OCV_values = OCV_SOC_25C(:, 2); % Corresponding OCV values


for j = 1:length(R0_values)
    R0 = R0_values(j);  
    OCV = zeros(n, 1);
    
    for i = 1:n  
        OCV(i) = V(i) - R0 * I(i); 
        SoC_OCV(i, j) = interp1(OCV_values, SOC_points, OCV(i), 'linear', 'extrap'); 
    end
end
figure; hold on
colors = {'b', 'g', 'r', 'm', 'k'}; 

for j = 1:length(R0_values)
    plot(t/3600, SoC_OCV(:, j), 'Color', colors{j}, 'LineWidth', 2, ...
        'DisplayName', ['R0 = ', num2str(R0_values(j)), ' Ω']);
end


SOC = zeros(n, 1);
SOC(1) = SOC_0;
for i = 1:n-1
    dt = t(i+1) - t(i);
    SOC(i+1) = SOC(i) + (dt * I(i) / (Q * 3600)) * 100;
end
plot(t/3600, SOC, 'k--', 'LineWidth', 3, 'DisplayName', 'True SoC');

xlabel('Time (h)', 'FontSize', 18); 
ylabel('SoC (%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
grid on;
legend('show', 'FontSize', 15);
title('SoC Estimation with Different R₀ Values');
