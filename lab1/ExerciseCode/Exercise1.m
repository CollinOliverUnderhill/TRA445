clear all
load ECM1

% Fit polynomial model to OCV data and plot comparison
Pocv = polyfit(OCV(:,1),OCV(:,2),10);
SoC = 0:0.01:1;
vOC = polyval(Pocv,SoC);

figure(1)
plot(OCV(:,1),OCV(:,2),'o',SoC,vOC,'-')
xlabel('SoC'), ylabel('OCV')

% At SoC = 0.5 we can calculate the OCV
vOC   = polyval(Pocv,0.5);

% Subtract the OCV and use Ohm's law to get R0 at each time instant
for k = 1:length(Time)    
    v0(k) = Voltage(k) - vOC;
    R0(k) = v0(k)/Current(k);
end
R0mean = mean(R0)

figure(2)
plot(Time,R0,'.')
xlabel('Time (s)'), ylabel('\it R_0')



%% Least square fit
Phi = Current;
Y   = v0';
R0hat = inv(Phi'*Phi)*Phi'*Y
figure(3)
plot(Current,v0,'.',Current,Current*R0hat,'r')
xlabel('Current'),ylabel('\it v_0')

% In the next exercise you will see why the estimate was not accurate