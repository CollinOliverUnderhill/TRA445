clear all;
load ECM2;

%% 4f: Parameter identification without using SoC/OCV information (with bias)
% Model: v(k) = theta0 + theta1*v(k-1) + theta2*i(k) + theta3*i(k-1)

Ts = 0.1;          % Sampling time (s)
N = length(Time);  % Number of data points

% Use raw terminal voltage (do not subtract vOC)
Y = Voltage;

% Construct regression matrix Phi:
% Each row: [1, v(k-1), i(k), i(k-1)]
Phi = [ ones(N-1,1), Y(1:N-1), Current(2:N), Current(1:N-1) ];

% RLS parameter settings
lambda = 1;  % Forgetting factor (try 0.99 or 0.995 to observe dynamic effects)
theta(:,1) = [0.01; 0.9; 0.0015; 0];  % Initial parameter guess (4 parameters)
P = diag([1, 1, 1, 1]);  % Initial covariance matrix

% Recursive Least Squares (RLS) main loop
for k = 2:(N-1)
    % Prediction error: fit Voltage(k) (Y(k))
    eps = Y(k) - Phi(k,:) * theta(:,k-1);
    % Compute gain vector
    denom = lambda + Phi(k,:) * P * Phi(k,:)';
    K = (P * Phi(k,:)') / denom;
    % Update parameters
    theta(:,k) = theta(:,k-1) + K * eps;
    % Update covariance matrix
    P = (1/lambda) * (eye(4) - K * Phi(k,:)) * P;
end

% Mapping: According to the following relationships:
%   theta(1) is the bias (not used for physical mapping),
%   theta(2) corresponds to alpha (alpha_hat),
%   theta(3) corresponds to R0 (R0_hat),
%   theta(4) is used to compute R1 (R1_hat) as:
%       R1_hat = (theta(4) + R0_hat * alpha_hat) / (1 - alpha_hat)
theta0_hat = theta(1);          % Bias (not used further)
alpha_hat  = theta(2);          % Estimated alpha
R0_hat     = theta(3);          % Estimated R0
R1_hat     = (theta(4) + R0_hat * alpha_hat) / (1 - alpha_hat);  % Estimated R1

% Estimate capacitance using the mapping: 
%   C_hat = -Ts / (R1_hat * ln(alpha_hat))
if (alpha_hat <= 0) || (alpha_hat >= 1) || (abs(log(alpha_hat)) < 1e-6)
    C_hat = inf;
else
    C_hat = -Ts / (R1_hat * log(alpha_hat));
end

% Plot: Display evolution of the 4 estimated parameters over time
figure;
plot(Time(1:end-1), theta, 'LineWidth', 1.2);
legend('theta0 (bias)','theta1 (\alpha)','theta2 (R0 part)','theta3 (R1 part)','Location','best');
xlabel('Time (s)');
ylabel('Parameter Estimates');
title('RLS Parameter Estimates with Bias (4f, no SoC/OCV)');
grid on;
axis tight;

% Simple Check: Display final mapped physical parameters
fprintf('\n----- Simple Check -----\n');
fprintf('R0_hat = %.6f ohm\n', R0_hat);
fprintf('R1_hat = %.6f ohm\n', R1_hat);
fprintf('C_hat  = %.2f F\n', C_hat);
