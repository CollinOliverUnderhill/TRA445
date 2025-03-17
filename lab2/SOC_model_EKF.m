function [x_hat_plus] = SOC_model_EKF(Data, OCV_SOC_25C, dOCV_SOC, x_P)
% Extended Kalman Filter to predict SOC of a Lithium-ion cell.
% Outputs SOC and transient voltage (V1) from the state estimates.

% % Assign Data
t = Data(:, 1);   % Time (s)
I = Data(:, 2);  % Current (A)
V = Data(:, 3); % Measured Voltage (V)
n = length(I);
SOC_initial = 0.8; % Initial SOC guess


%%% simulate measured current signal by adding noises
noise_std = 0;             % Standard deviation of the noise
noise_mean = 0;            % Mean of the noise
noise = noise_std * randn(size(I)) + noise_mean;   % Generate noise with a non-zero mean
I = I + noise;             % Add noise to the signal 
%%%
% 
% 
% Model Parameters
R0 = x_P(1); R1 = x_P(2); C1 = x_P(3);
capacity = 2;              % Battery capacity (Ah)

% Kalman Initialization
V1_0 = 0;                  % Initial transient voltage guess (V1)
SOC_0 = SOC_initial;       % Initial SOC guess
x_hat_0 = [V1_0; SOC_0];   % Initial state

K_lfx = 1;
P0 = diag([5e-5, 5e-2]);   % Initial estimation error covariance
Q = diag([1e-6, 1e-5]);    % Process noise covariance
R = (6e-2);                  % Measurement noise variance

x_hat_plus = zeros(2, n);  % Preallocate state estimates
K = zeros(2, n);           % Preallocate Kalman gain
err = zeros(1, n);         % Preallocate error vector

% EKF Loop
for i = 1:n

% time interval    
   if i == 1; dt = t(i);
   else; dt = t(i)-t(i-1); end

    % State transition and input matrices
    A = [exp(-dt / (R1 * C1)), 0; 0, 1];
    B = [R1 * (1 - exp(-dt / (R1 * C1))); dt / (capacity * 3600)];

% Prediction step
   if i == 1 
      x_hat_minus = A*x_hat_0 + B*I(i);             % State estimate.  
      P_minus = A*P0*A' + Q;                        % Estimation-error covariance.
   else
      x_hat_minus = A*x_hat_plus(:,i-1) + B*I(i);   % State estimate.  
      P_minus = A*P_plus*A' + Q;                    % Estimation-error covariance.  
   end


% Update step
    % Interpolate OCV and dOCV
    SOC = x_hat_minus(2);
    OCV = interp1(OCV_SOC_25C(:, 1) / 100, OCV_SOC_25C(:, 2), SOC, 'spline');
    dOCV = interp1(dOCV_SOC(:, 1) / 100, dOCV_SOC(:, 2), SOC, 'spline');

     C = [1, dOCV];
%   C = [1, 0];

    % Kalman gain
    K(:, i) = P_minus * C' / (C * P_minus * C' + R);

   % Measurement model    
    y_est = OCV + R0 * I(i) + x_hat_minus(1);

    % Measurement update
    err(i) = V(i) - y_est;      % Voltage prediction error
    x_hat_plus(:, i) = x_hat_minus + K(:, i) * err(i);   % update the state estimate
    P_plus = (eye(2) - K(:, i) * C) * P_minus;           % update the error covariance
    Pv1(i) = P_plus(1,1);
    Psoc(i) = P_plus(2,2);
end

end
% 
% 



%  %% Extended Kalman filtering '2RC model'
% % 2RC model
%   R0 = x_P(1); R1 = x_P(2); C1 = x_P(3); R2 = x_P(4); C2 = x_P(5); 
% % Kalman initialization and internal parameters
%   V1_0 = 0;                               % Initial guess of transient voltage (V1).
%   V2_0 = 0;                               % Initial guess of transient voltage (V2).
%   SOC_0 = SOC_initial;                    % Initial guess of SOC.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   x_hat_0 = [V1_0; V2_0; SOC_0];          % Initial state.
% 
%   P0 = [5e-5 0 0; 0 5e-6 0; 0 0 2e-3];    % Initial estimation error covariance.
% 
%   Q = [1e-6 0 0; 0 1e-8 0; 0 0 1e-10];    % Initial process noise covariance.
%   R = 1e-4;                               % Measurement variance (covariance).
%   capacity = 2;                        % (Ah), set as unchanged.
%   [D_FUDS, D_HDS, D_BJDST] = Read_dynamic_data();
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  for i = 1:D_BJDST   %L_FUDS, L_HDS, L_BJDST% Filtering over all selected time steps.
% 
%    % Time update
%    if i == 1 
%       dt = t(i);
%    else
%       dt = t(i)-t(i-1);
%    end
% 
%    F = [exp(-dt/(R1*C1)) 0 0; 0 exp(-dt/(R2*C2)) 0; 0 0 1];  % Define the matrix of F, F will be a constant matrix if R1 and C1 are unchanged, in this code is unchanged.
%    G = [ R1 * (1-exp(-dt/(R1*C1))) ; R2 * (1-exp(-dt/(R2*C2))) ; dt/(capacity*3600) ];  % Define the matrix of G.
% 
%    % Time update of the state estimate and estimation-error covariance. 
%    % Time update of the estimation-error covariance.
%    if i == 1 
%       P_minus = F*P0*F' + Q; % Estimation-error covariance.
%    else
%       P_minus = F*P_plus*F' + Q; % Estimation-error covariance.  
%    end
% 
%    % Time update of the state estimate.
%    if i == 1  
%       x_hat_minus = F*x_hat_0 + G*I(i);           % State estimate.  
%    else
%       x_hat_minus = F*x_hat_plus(:,i-1) + G*I(i); % State estimate.   
%    end
% 
%  % Compute the following partial derivative matrices. 
%    SOC = x_hat_minus(3);
%    x = OCV_SOC_25C(:,1)/100;             % SOC (%)
%    y = OCV_SOC_25C(:,2);                 % OCV (V)
%    OCV = interp1 (x, y, SOC, 'spline');  % Get the corresponding OCV from the above points.   
% 
%    x1 = dOCV_SOC(:,1)/100;                  % SOC (%)
%    y1 = dOCV_SOC(:,2);                      % OCV (V)
%    dOCV = interp1 (x1, y1, SOC, 'spline');  % Get the corresponding OCV from the above points.   
% 
%    H = [1, 1, dOCV]; % Calculating H matrix. H is changed all the time, because SOC is changed.
% 
%    y(i) = OCV + R0*I(i) + x_hat_minus(1) + x_hat_minus(2);  % Output equation, 'y' is 'terminal voltage'.
% 
% % Measurement update.
%    K(:,i) = P_minus * H' * ( (H * P_minus * H' + R )^-1); % Kalman gain matrix.
% 
%    err(i) = V(i)-y(i); % Error of predicted 'terminal voltage'.
%    % 'step.voltage(i)' is measure result.
%    % 'y(i)' is estimate result.  
% 
%    % State estimate measurement update.
%    x_hat_plus(:,i) = x_hat_minus + K(:,i) * err(i); 
% 
%    % Estimation error covariance measurement update.
%    % P_plus = ( eye(2) - K(:,i)*H )*P_minus * ( eye(2) - K(:,i)*H )' + K(:,i)*R*(K(:,i))'; 
%    P_plus = ( eye(3) - K(:,i)*H )*P_minus;
%  end

 





