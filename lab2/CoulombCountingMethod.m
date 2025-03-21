%% SoC estimation using Coulomb Counting method
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

Q = 2;
dSOC = (diff(t(1053:1196)) / 3600) .* (I(1053:1195)+ I(1054:1196)) / 2 / Q * 100;  % SOC increment
SOC_0 = 100 + sum(dSOC);  % initial SOC

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
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size


% Apply the Coulomb Counting method
% From step 1917, the DST start
ST = 1917;
OCV = V (ST);  % The OCV is equal the terminal voltage after a long resting time

I = I(ST:end);            % set as ideal current signal
V = V(ST:end);
t = t(ST:end) - t(ST); % Adjust time to start from 0
Q = 2;                       % Cell capacity

%%% simulate measured current signal by adding noises
noise_std = 0.05;            % Standard deviation of the noise
noise_mean = 0.0001;        % Mean of the noise
noise = noise_std * randn(size(I)) + noise_mean;   % Generate noise with a non-zero mean
I_mea = I + noise;        % Add noise to the signal 

%%%


% Interpolate initial SoC from OCV-SOC relationship
SOC_points = OCV_SOC_25C(:, 1); % SOC sample points
OCV_values = OCV_SOC_25C(:, 2); % Corresponding OCV values
SOC0_est = interp1(OCV_values, SOC_points,  OCV, 'linear');

% Initialize variables
n = length(I);
SOC = zeros(n, 1);           % true SoC
SOC(1) = SOC_0;            % initial SoC

SOC_est = zeros(n, 1);         % estimated SoC
SOC_est(1) = SOC0_est;     % initial SoC

% Coulomb counting
for i = 1:n-1
    dt = t(i+1) - t(i); % Time step
    SOC(i+1) = SOC(i) + (dt * I(i) / (Q * 3600)) * 100;                % measured/reference SoC
    SOC_est(i+1) = SOC_est(i) + (dt * I_mea(i) / (Q * 3600)) * 100;   % estimated SoC
end

% Plot the result
figure; hold on

plot(t/3600, SOC_est, 'LineWidth', 2, 'DisplayName', 'SoC estimate');
plot(t/3600, SOC, 'LineWidth', 3, 'DisplayName', 'True');

xlabel('Time (h)', 'FontSize', 18); 
ylabel('SoC (%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);

figure; hold on

plot(t/3600, (SOC_est-SOC)./SOC, 'LineWidth', 2, 'DisplayName', 'SoC error');


xlabel('Time (h)', 'FontSize', 18); 
ylabel('Error(%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);
