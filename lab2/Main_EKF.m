% Extended Kalman Filter for State of Charge (SOC) Estimation
clear; clc;

% Load and Read Data
% Load dynamic driving data
[D_FUDS, D_HDS, D_BJDST] = Read_dynamic_data();
%   D_FUDS: Federal Urban Driving Schedule data [Time(s), Current(A), Voltage(V)]
%   D_HDS: Highway Driving Schedule data
%   D_BJDST: Beijing Dynamic Stress Test data

% Load OCV-SOC and build the dOCV-SOC relation
load('OCV_SOC_relation.mat');   % OCV-SOC look-up table (25°C)

dOCV_SOC = dOCV_SOC();      % dOCV-SOC look-up table

% Compute SOC using Coulomb counting (experimental results)
[SOC_FUDS, SOC_HDS, SOC_BJDST] = SOC_measured(D_FUDS, D_HDS, D_BJDST);   % Effective range: 10-80%

% Model parameters (1RC model)
 x_P = [0.070248, 0.009953, 885.996888];    % PSO 1
%x_P = [0.07219, 0.01193, 94836.46270];    % PSO bounds from pulse test
% x_P = [0.07259, 0.03197, 1413.22938];    % from Least Squares method
% x_P = [0.08191, 0.02386, 47418.23135];  % from pulse test method

% x_P = [0.06962;0.00949;987.42575;0.00111;655.93623]; %2RC
% Extended Kalman Filter Implementation
% SOC estimation using EKF
SOCdata = SOC_FUDS;       % SOC_FUDS, SOC_HDS, or SOC_BJDST, Measured SOC
D = D_FUDS;                        % D_FUDS, D_HDS, or D_BJDST, the used data
[x_hat_plus] = SOC_model_EKF(D, OCV_SOC_25C, dOCV_SOC, x_P);

% Extract and process results
SOC_model = x_hat_plus(2, :)' * 100;    % SOC from model (as percentage)


t = SOCdata(:, 1);            % Time (s)
SOC_measure = SOCdata(:, 2);     % SOC_FUDS, SOC_HDS, or SOC_BJDST, Measured SOC
% Calculate Root Mean Square Error (RMSE)
n = length(SOC_model);    
RMSE = sqrt(mean((SOC_measure(1:n) - SOC_model(1:n)).^2));
% MAE = mean((SOC_measure(1:n) - SOC_model(1:n)));
% absolute_error =     abs(y_true - y_pred);
MAE = mean( abs(SOC_measure(1:n) - SOC_model(1:n)));


% Plot Results
figure;
plot(t, SOC_measure, 'LineWidth', 2, 'DisplayName', 'Measured SOC');
hold on;
plot(t(1:n), SOC_model(1:n), '--', 'LineWidth', 2, 'DisplayName', 'Estimated SOC');
hold off;
grid on; 
 % xlim([0 2000])
xlabel('Time (s)', 'FontSize', 18);
ylabel('SOC (%)', 'FontSize', 18);
title(['SOC Estimation with EKF (MAE = ' num2str(MAE, '%.2f') '%)']);
legend('show', 'FontSize', 15);
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 15);