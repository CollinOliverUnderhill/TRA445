clear all
load ECM1

% Method 1

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

% Least square fit
Y   = v0;
Phi = Current;
R0hat = inv(Phi'*Phi)*Phi'*Y
figure(3)
plot(Current,v0,'.',Current,Current*R0hat,'r')
xlabel('Current'),ylabel('\it v_0')


% Method 2
% Assuming dvOC/dz = K the regression model can be written as
%
% y(k) = v(k) - v(k-1) = [i(k)  -i(k-1)]*[R0  R0-K*Ts/Q]'
%

N = length(Time);
Y = Voltage(2:N)-Voltage(1:N-1);
Phi = [Current(2:N) -Current(1:N-1)];

theta = inv(Phi'*Phi)*Phi'*Y