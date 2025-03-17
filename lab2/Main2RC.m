% Extended Kalman Filter for State of Charge (SOC) Estimation using 2RC model
clear; clc;

% Load and Read Data
% Load dynamic driving data
[D_FUDS, D_HDS, D_BJDST] = Read_dynamic_data();
%   D_FUDS: Federal Urban Driving Schedule data [Time(s), Current(A), Voltage(V)]
%   D_HDS: Highway Driving Schedule data
%   D_BJDST: Beijing Dynamic Stress Test data

% Load OCV-SOC and build the dOCV-SOC relation
load('OCV_SOC_relation.mat');   % OCV-SOC look-up table (25°C)
dOCV_SOC = dOCV_SOC();           % dOCV-SOC look-up table

% Compute SOC using Coulomb counting (experimental results)
[SOC_FUDS, SOC_HDS, SOC_BJDST] = SOC_measured(D_FUDS, D_HDS, D_BJDST);   % Effective range: 10-80%

% Model parameters (2RC model)
% Parameter vector: [R0; R1; C1; R2; C2]
x_P = [0.06962; 0.00949; 987.42575; 0.00111; 655.93623];

% Extended Kalman Filter Implementation using 2RC model
% Choose dataset: SOC_FUDS, SOC_HDS, or SOC_BJDST for measured SOC and corresponding data
SOCdata = SOC_FUDS;       % Measured SOC from Coulomb counting
D = D_FUDS;               % Data for the EKF (choose D_FUDS, D_HDS, or D_BJDST)
[x_hat_plus] = SOC_model_EKF_2RC(D, OCV_SOC_25C, dOCV_SOC, x_P);

% Extract and process results
% In 2RC model, state vector is [V1; V2; SOC]
SOC_model = x_hat_plus(3, :)' * 100;    % Estimated SOC (as percentage)

t = SOCdata(:, 1);            % Time (s)
SOC_measure = SOCdata(:, 2);  % Measured SOC

% Calculate RMSE and MAE
n = length(SOC_model);
RMSE = sqrt(mean((SOC_measure(1:n) - SOC_model(1:n)).^2));
MAE = mean(abs(SOC_measure(1:n) - SOC_model(1:n)));

% Plot Results
figure;
plot(t, SOC_measure, 'LineWidth', 2, 'DisplayName', 'Measured SOC');
hold on;
plot(t(1:n), SOC_model(1:n), '--', 'LineWidth', 2, 'DisplayName', 'Estimated SOC');
hold off;
grid on;
xlabel('Time (s)', 'FontSize', 18);
ylabel('SOC (%)', 'FontSize', 18);
title(['SOC Estimation with EKF (MAE = ' num2str(MAE, '%.2f') '%)']);
legend('show', 'FontSize', 15);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 15);
