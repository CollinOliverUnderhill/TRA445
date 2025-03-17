clear all
load ECM1

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

v0 = Voltage - vOC';

% RLS algorithm
Y   = v0;
Phi = Current;
lambda = 0.9999     % Try different values of lambda (forgetting factor) 

theta(1) = 0.001;   % Initial parameter guess
P(1)     = 1;    % Initialization of P
for k = 2:length(Time)
    eps      = Y(k,:) - Phi(k,:)*theta(k-1);
    K        = P(k-1)*Phi(k,:)/(lambda+Phi(k,:)*P(k-1)*Phi(k,:)');
    theta(k) = theta(k-1) + K*eps;
    P(k)     = (1/lambda)*(1-K*Phi(k,:))*P(k-1);
end

figure(4)
plot(Time,theta)