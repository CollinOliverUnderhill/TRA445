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
Q  = 60*3600; % Capacity (As)
R0o = 0.001;   % Ohmic resistance (Ohm)
R1o = 0.0015;  % Charge transfer resistance (Ohm)
Co  = 10000;   % Dual layer capacistance (Farad)

% Calculate SoC
% z(1) = 0.5; % Initial SoC
% for k = 1:length(Time)-1
%     z(k+1) = z(k)+(Ts/Q/3600)*Current(k);
% end

% Voltage sensor noise (normally distributed)
w = 0.01*randn(size(Time));

% Let the parameters change
R0 = R0o*(1+0.2*sin(Time/3600*10));
R1 = R1o*(1+0.0*sin(Time/3600*10));
C  = Co*(1+0*sin(Time/3600*10));
% Calculate voltage for known z
v1(1)  = Current(1)*R1(1);      % Initial polarization voltage (SS assumed)
z(1)   = 0.5;                % Initial SoC
x(:,1) = [v1(1);z(1)];       % Initial state
vOC    = polyval(Pocv,z(1)); % Initial OCV
v(1)   = vOC+R0(1)*Current(1)+v1(1);  % Initial terminal voltage

for k  = 2:length(Time)
    Ad = [exp(-Ts/(R1(k)*C(k))) 0;0 1];
    Bd = [R1(k)*(1-exp(-Ts/(R1(k)*C(k)))); Ts/Q];

    x(:,k) = Ad*x(:,k-1)+Bd*Current(k-1);
    v1   = x(1,k);
    z(k) = x(2,k);
    vOC  = polyval(Pocv,z(k));             % OCV at current SoC
    v(k) = vOC+R0(k)*Current(k)+v1+w(k);      % Terminal voltage
end
Voltage = v';

save('ECM2dyn','Time','Current','Voltage','OCV','z','Ts','R0','R1','C')

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

