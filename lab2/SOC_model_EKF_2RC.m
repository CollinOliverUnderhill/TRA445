function [x_hat_plus] = SOC_model_EKF_2RC(Data, OCV_SOC_25C, dOCV_SOC, x_P)
% Extended Kalman Filter for SOC estimation using a 2RC model.
% This implementation estimates both transient voltages (V1, V2) and SOC.
%
% Inputs:
%   Data         - [Time, Current, Voltage] matrix.
%   OCV_SOC_25C  - OCV-SOC lookup table (for interpolation) at 25°C.
%   dOCV_SOC     - dOCV-SOC lookup table (for interpolation).
%   x_P          - Parameter vector [R0; R1; C1; R2; C2].
%
% Outputs:
%   x_hat_plus - State estimates [V1; V2; SOC] for each time step.

% Extract data
t = Data(:, 1);   % Time (s)
I = Data(:, 2);   % Current (A)
V = Data(:, 3);   % Measured Voltage (V)
n = length(I);
SOC_initial = 0.8; % Initial SOC guess

% 2RC model parameters
R0 = x_P(1); R1 = x_P(2); C1 = x_P(3); R2 = x_P(4); C2 = x_P(5);
capacity = 2;  % Battery capacity (Ah)

% Kalman Filter Initialization
V1_0 = 0;       % Initial transient voltage guess (V1)
V2_0 = 0;       % Initial transient voltage guess (V2)
SOC_0 = SOC_initial;
x_hat_0 = [V1_0; V2_0; SOC_0];  % Initial state

P0 = [5e-5 0 0; 0 5e-6 0; 0 0 2e-3];    % Initial estimation error covariance
Q = [1e-6 0 0; 0 1e-8 0; 0 0 1e-10];     % Process noise covariance
R = 1e-4;                                % Measurement noise variance

% Preallocate variables
x_hat_plus = zeros(3, n);   % State estimates: [V1; V2; SOC]
K = zeros(3, n);            % Kalman gain
err = zeros(1, n);          % Measurement error
V_est = zeros(n, 1);        % Estimated terminal voltage
P_plus = P0;                % Initialize error covariance

% EKF Loop
for i = 1:n
    % Time step calculation
    if i == 1
        dt = t(i);
    else
        dt = t(i) - t(i-1);
    end
    
    % State transition matrix F and input matrix G for 2RC model
    F = [exp(-dt/(R1*C1)), 0, 0;
         0, exp(-dt/(R2*C2)), 0;
         0, 0, 1];
    G = [R1 * (1 - exp(-dt/(R1*C1)));
         R2 * (1 - exp(-dt/(R2*C2)));
         dt/(capacity*3600)];
     
    % Time update (Prediction)
    if i == 1
        x_hat_minus = F * x_hat_0 + G * I(i);
        P_minus = F * P0 * F' + Q;
    else
        x_hat_minus = F * x_hat_plus(:, i-1) + G * I(i);
        P_minus = F * P_plus * F' + Q;
    end
    
    % Interpolate OCV and dOCV based on current SOC estimate
    SOC = x_hat_minus(3);
    soc_range = OCV_SOC_25C(:, 1) / 100;  % Convert percentage to fraction
    ocv_values = OCV_SOC_25C(:, 2);
    OCV_val = interp1(soc_range, ocv_values, SOC, 'spline');
    
    soc_range_d = dOCV_SOC(:, 1) / 100;
    docv_values = dOCV_SOC(:, 2);
    dOCV_val = interp1(soc_range_d, docv_values, SOC, 'spline');
    
    % Measurement model:
    % Terminal voltage: V_est = OCV + R0 * I + V1 + V2
    V_est(i) = OCV_val + R0 * I(i) + x_hat_minus(1) + x_hat_minus(2);
    
    % Observation matrix H (Jacobian of measurement function)
    H = [1, 1, dOCV_val];
    
    % Kalman gain calculation
    K(:, i) = P_minus * H' / (H * P_minus * H' + R);
    
    % Measurement update
    err(i) = V(i) - V_est(i);  % Voltage prediction error
    x_hat_plus(:, i) = x_hat_minus + K(:, i) * err(i);
    P_plus = (eye(3) - K(:, i) * H) * P_minus;
end

end
