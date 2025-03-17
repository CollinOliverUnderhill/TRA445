clear all
% Creates data for the exercises with gradual parameter changes

% Read drive cycle data (Trips in BMW i3, see 
% https://ieee-dataport.org/open-access/battery-and-heating-data-real-driving-cycles)
Data_table = readtable('TripA27.csv'); 

Time    = Data_table{:,1};  % (s)
Current = Data_table{:,9};  % (A)
% Voltage = Data_table(:,8);  % (V)

% Sampling time
Ts = 0.1;                   % (s)

% Read OCV curve
load LG_ocv;

% Format OCV data [SOC, OCV] (SoC in (0,1)!)
OCV = [LG_soc/100 mean(LG_ocv_charge,2)];

% Fit polynomial model to OCV data
Pocv = polyfit(OCV(:,1),OCV(:,2),10);

% ---------------------------
% Define parameters that vary gradually over time
% 初始参数与结束参数
R0_start = 0.001;   R0_end = 0.002;    % Ω
R1_start = 0.0015;  R1_end = 0.002;      % Ω
C_start  = 10000;   C_end  = 8000;       % F

% ---------------------------
% Capacity
Q  = 60*3600; % (As)

% Voltage sensor noise (normally distributed)
w = 0.00001*randn(size(Time));

% Initialize state variables
z(1)   = 0.5;                % Initial SoC
% For initial polarization voltage, use initial R1 value:
R1_current = R1_start; 
v1(1)  = Current(1)*R1_current;   % Initial polarization voltage (assumed steady state)
x(:,1) = [v1(1); z(1)];      % Initial state
vOC    = polyval(Pocv,z(1)); % Initial OCV
v(1)   = vOC + R0_start*Current(1) + v1(1);  % Initial terminal voltage (using R0_start)

% Loop over time steps
n = length(Time);
for k = 2:n
    % Compute a normalized time factor from 0 to 1
    tau = (k-1) / (n-1);
    
    % Linearly interpolate parameter values
    R0 = R0_start + (R0_end - R0_start) * tau;
    R1 = R1_start + (R1_end - R1_start) * tau;
    C  = C_start  + (C_end  - C_start ) * tau;
    
    % Compute discrete-time state-space matrices using current parameters
    Ad = [exp(-Ts/(R1*C))  0; 0 1];
    Bd = [R1*(1-exp(-Ts/(R1*C))); Ts/Q];
    Cd = [1 0];
    Dd = R0;
    
    % State update: x(k) = Ad*x(k-1) + Bd*Current(k-1)
    x(:,k) = Ad*x(:,k-1) + Bd*Current(k-1);
    v1   = x(1,k);  % Polarization voltage
    z(k) = x(2,k);  % SoC (assumed measured or estimated)
    vOC  = polyval(Pocv,z(k));  % OCV at current SoC
    
    % Terminal voltage: v(k) = vOC + R0*Current(k) + v1 + measurement noise
    v(k) = vOC + R0*Current(k) + v1 + w(k);
end

Voltage = v';

save('ECM2n','Time','Current','Voltage','OCV','z','Ts','R0','R1','C')

% Optional: Plot results
figure;
subplot(3,1,1)
plot(Time, z)
xlabel('Time (s)')
ylabel('SoC')
title('SoC over time')

subplot(3,1,2)
plot(Time, Voltage)
xlabel('Time (s)')
ylabel('Voltage (V)')
title('Terminal Voltage')

subplot(3,1,3)
% 画出参数随时间的变化（示例只画 R0）
R0_values = R0_start + (R0_end - R0_start) * ((1:n)-1)'/(n-1);
plot(Time, R0_values)
xlabel('Time (s)')
ylabel('R0 (Ohm)')
title('R0 variation over time')

