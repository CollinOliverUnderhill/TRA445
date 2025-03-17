clear all;
load ECM2;

% --------------------------------------------------------------------
% 4f: Parameter identification without using SoC/OCV information,
%     but with an added bias term.
%     Regression model:
%       v(k) = theta0 + theta1*v(k-1) + theta2*i(k) + theta3*i(k-1)
% --------------------------------------------------------------------

Ts = 0.1;              % Sampling time (s)
N = length(Time);      % Data length

% Use the raw terminal Voltage (do not subtract OCV)
Y = Voltage;

% Construct regression matrix Phi:
% Each row corresponds to a data point with format:
% [1, v(k-1), i(k), i(k-1)]
Phi = [ ones(N-1,1), Y(1:N-1), Current(2:N), Current(1:N-1) ];

% RLS parameter settings
lambda = 1;    % Forgetting factor (try 0.99 or 0.995 to observe dynamic effects)

% Preallocate and set initial guess for theta (4 parameters)
% theta = [theta0; theta1; theta2; theta3]
theta = zeros(4, N-1);
theta(:,1) = [0; 0.9; 0.0015; -0.0015];

% Initialize covariance matrix P
P = eye(4);

% Recursive Least Squares (RLS) update loop
for k = 2:(N-1)
    % Prediction error: target is Voltage(k)
    eps_k = Y(k) - Phi(k,:) * theta(:, k-1);
    % Compute gain vector
    denom = lambda + Phi(k,:) * P * Phi(k,:)';
    K = (P * Phi(k,:)') / denom;
    % Update theta
    theta(:, k) = theta(:, k-1) + K * eps_k;
    % Update covariance matrix
    P = (1/lambda) * (eye(4) - K * Phi(k,:)) * P;
end

% Plot the evolution of the 4 estimated parameters over time
figure(5)
plot(Time(1:end-1), theta, 'LineWidth', 1.5);
legend('theta0 (bias)','theta1','theta2','theta3','Location','best');
xlabel('Time (s)');
ylabel('Parameter Estimates');
title('RLS Parameter Estimates with Bias (No SoC/OCV)');
grid on;
axis tight;

% --------------------------------------------------------------------
% Simple Check: Map final estimated theta to physical parameters
% --------------------------------------------------------------------
% For the final estimated theta, let:
%   theta(1): bias (ignored for mapping)
%   theta(2): estimated alpha, i.e., alpha_hat
%   theta(3): estimated R0 (R0_hat)
%   theta(4): used to compute R1_hat via:
%             R1_hat = (theta(4) + R0_hat * alpha_hat) / (1 - alpha_hat)
%
% And the mapping for capacitance is:
%   C_hat = -Ts / (R1_hat * log(alpha_hat))
%
% Note: The original mapping formulas are designed for the case when OCV is removed.
% Here we simply use them for comparison.
theta_final = theta(:, end);
alpha_hat_est = theta_final(2);
R0_hat_est    = theta_final(3);
R1_hat_est    = (theta_final(4) + R0_hat_est * alpha_hat_est) / (1 - alpha_hat_est);
if (alpha_hat_est <= 0) || (alpha_hat_est >= 1) || (abs(log(alpha_hat_est)) < 1e-6)
    C_hat_est = 1e9;
else
    C_hat_est = -Ts / (R1_hat_est * log(alpha_hat_est));
end

% For comparison, compute the theoretical alpha (from the original ECM2data parameters)
% Note: R1 and C here are from the ECM2 data file (as loaded), i.e., the nominal values.
alpha_theo = exp(-Ts/(R1*C));

fprintf('\n----- Simple Check -----\n');
fprintf('Estimated alpha_hat: %.6f\n', alpha_hat_est);
fprintf('Theoretical alpha  : %.6f\n', alpha_theo);
fprintf('Estimated R0_hat   : %.6f ohm\n', R0_hat_est);
fprintf('Estimated R1_hat   : %.6f ohm\n', R1_hat_est);
fprintf('Estimated C_hat    : %.2f F\n', C_hat_est);
