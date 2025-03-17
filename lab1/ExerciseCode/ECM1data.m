clear all
% Creates data for Exercise 1

% Read drive cycle data (Trips in BMW i3, see 
% https://ieee-dataport.org/open-access/battery-and-heating-data-real-driving-cycles
Data_table = readtable('TripA27.csv'); 

%save as struct
%Data_struct = table2struct(Data_table); 

Time    = Data_table{:,1};  % (s)
Current = Data_table{:,9};  % (A)
% Voltage = Data_table(:,8);  % (V)

% Sampling time
Ts = 0.1;                   % (s)

% Read OCV curve
load LG_ocv;

% To format [SOC, OCV] (SoC in (0,1)!)
OCV = [LG_soc/100 mean(LG_ocv_charge,2)];

% Fit polynomial model to OCV data
Pocv = polyfit(OCV(:,1),OCV(:,2),10);

% Parameters
Q  = 60;      % Capacity (Ah)
R0 = 0.0025;  % Ohmic resistance (Ohm)
% R1 = 0.001; % Charge transfer resistance (Ohm)
% C1 = 10000; % Dual layer capacistance (Farad)

% Calculate SoC
z(1) = 0.5; % Initial SoC

for k = 1:length(Time)-1
    z(k+1) = z(k)+(Ts/Q/3600)*Current(k);
end

% Voltage sensor noise (normally distributed)
w = 0.01*randn(size(Time));

% Calculate voltage based only on R0 and OCV polynomial
for k = 1:length(Time)
    vOC  = polyval(Pocv,z(k));           % OCV at current SoC
    v(k) = vOC + R0*Current(k)+w(k);     % Terminal voltage
end

Voltage = v';
save('ECM1','Time','Current','Voltage','OCV')

% Commands for viewing
% plot(OCV(:,1),OCV(:,2),'o')
% SoC = 0:0.01:1;
% vOC = polyval(Pocv,SoC);
% hold on
% plot(SoC,vOC)
% hold off

% plot(Time,z)

% plot(Time,Vt)
% 
% for k = 1:length(Time)
%     vOC   = polyval(Pocv,0.5);
%     v0(k) = Vt(k) - vOC;
%     R0(k) = v0(k)/Current(k);
% end
% 
% plot(Time,R0,'.')  

    
%plot(OCV(:,1),OCV(:,2))



%plot(Time,Current)

