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
figure;
subplot(2, 1, 1)
plot(t/3600, I, 'LineWidth', 2)
xlabel('Time (h)', 'FontSize', 18); 
ylabel('Current (A)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
xlim ([0 9])

subplot(2, 1, 2)
plot(t/3600, V, 'LineWidth', 2)
xlabel('Time (h)', 'FontSize', 18) 
ylabel('Voltage (V)', 'FontSize', 18)
xlim ([0 9])

set(gcf, 'Color', 'w'); 
grid on;                            % Enable grid
set(gca, 'FontSize', 15);  % Set axis font size
%%
% Apply the OCV method
% Select data from the start index, the data after fully charged
ST = 333;
I = I(ST:end);
V = V(ST:end);
t = t(ST:end) - t(ST); % Adjust time to start from 0
SOC_0 = 100;          % Initial SoC
Q = 2;                       % Cell capacity
R0 = 0.08191;          % Ohmic Resistance, obtained from pulse test 

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

% OCV Method Simulate
for i = 1:n  

    % OCV calculation based on the model
     OCV(i) = V(i) - R0 * I(i);

    % Interpolate SoC from OCV-SOC relationship
    SoC_est(i) = interp1(OCV_values, SOC_points,  OCV(i), 'linear','extrap');

end

%%
% Plot the result
figure; hold on

plot(t/3600, SoC_est, 'LineWidth', 2, 'DisplayName', 'SoC estimate');
plot(t/3600, SOC, 'LineWidth', 3, 'DisplayName', 'True');

xlabel('Time (h)', 'FontSize', 18); 
ylabel('SoC (%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);
%%

figure; hold on

plot(t/3600, V, 'LineWidth', 2, 'DisplayName', 'Terminal voltage');
plot(t/3600, OCV, 'LineWidth', 2, 'DisplayName', 'OCV');

xlabel('Time (h)', 'FontSize', 18); 
ylabel('Voltage (V)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);







