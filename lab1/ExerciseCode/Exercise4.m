clear all
load ECM2n.mat

% Recursive least squares for Method 1 in Exercise 2

% Fit polynomial model to OCV data
Pocv = polyfit(OCV(:,1),OCV(:,2),10);

% Open loop estimation of z
Ts = 0.1;           % Sampling time (s)
Q  = 60*3600;       % Capacity (As)
z(1) = 0.5;                                 % Initial SoC
vOC(1) = polyval(Pocv,z(1));                % Inital OCV
for k = 1:length(Time)-1
    z(k+1) = z(k)+(Ts/Q)*Current(k);        % Simulated SoC
    vOC(k+1) = polyval(Pocv,z(k+1));        % Corresponding OCV
end

% Standard LSQ
N   = length(Time);
Y   = Voltage - vOC';
Phi = [Y(1:N-1) Current(2:N) Current(1:N-1)];

theta = inv(Phi'*Phi)*Phi'*Y(2:N);

% Simple check
alpha_hat = theta(1)
alpha     = exp(-Ts/(R1*C))

R0_hat = theta(2)
R0

R1_hat = (theta(3)+theta(2)*theta(1))/(1-theta(1))
R1

C_hat  = -1/(R1_hat*log(alpha_hat))
% C_hat  = (-1/(R1_hat*log(alpha_hat)))/10


