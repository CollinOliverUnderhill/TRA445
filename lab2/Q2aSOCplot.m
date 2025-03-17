clc;clear all;
SOCmea01 = load('SOCestinitial50.mat')
SOCmea1 = load('SOCestinitial120.mat')
datatrue = load('turedata')
SOC_50 = SOCmea01.SOC_est;
SOC_st1000 = SOCmea1.SOC_est;
SOC = datatrue.SOC;
SOC_100 = datatrue.SOC_est
t=datatrue.t;

figure; hold on

plot(t/3600, SOC_100, 'LineWidth', 2, 'DisplayName', 'SoC Original');
plot(t/3600, SOC, 'LineWidth', 3, 'DisplayName', 'True');
plot(t/3600, SOC_50, 'LineWidth', 2, 'DisplayName', 'SoC initial lower');
plot(t/3600, SOC_st1000, 'LineWidth', 3, 'DisplayName', 'Soc initial higher');
xlabel('Time (h)', 'FontSize', 18); 
ylabel('SoC (%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);

%%

clc;clear all;
SOCmea01 = load('SOCnomea01.mat')
SOCmea1 = load('SOCnomea1.mat')
datatrue = load('turedata')
SOC_st100 = SOCmea01.SOC_est;
SOC_st1000 = SOCmea1.SOC_est;
SOC = datatrue.SOC;
SOC_100 = datatrue.SOC_est
t=datatrue.t;

figure; hold on

plot(t/3600, SOC_100, 'LineWidth', 2, 'DisplayName', 'noise-mean = 0.0001');
plot(t/3600, SOC, 'LineWidth', 3, 'DisplayName', 'True');
plot(t/3600, SOC_st100, 'LineWidth', 2, 'DisplayName', 'noise-mean = 0.1');
plot(t/3600, SOC_st1000, 'LineWidth', 3, 'DisplayName', 'noise-mean = 1');
xlabel('Time (h)', 'FontSize', 18); 
ylabel('SoC (%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);

%%

clc;clear all;
SOCst100 = load('SOCnost100.mat')
SOCst1000 = load('SOCnost1000.mat')
datatrue = load('turedata')
SOC_st100 = SOCst100.SOC_est;
SOC_st1000 = SOCst1000.SOC_est;
SOC = datatrue.SOC;
SOC_o = datatrue.SOC_est
t=datatrue.t;

figure; hold on

plot(t/3600, SOC_o, 'LineWidth', 2, 'DisplayName', 'noise-std = 0.05');
plot(t/3600, SOC, 'LineWidth', 3, 'DisplayName', 'True');
plot(t/3600, SOC_st100, 'LineWidth', 2, 'DisplayName', 'noise-std = 5');
plot(t/3600, SOC_st1000, 'LineWidth', 3, 'DisplayName', 'noise-std = 50');
xlabel('Time (h)', 'FontSize', 18); 
ylabel('SoC (%)', 'FontSize', 18); 
set(gca, 'FontSize', 15);
 % xlim ([4 8])
set(gcf, 'Color', 'w'); 
grid on;  % Enable grid
set(gca, 'FontSize', 15); % Set axis font size
legend('show', 'FontSize', 15);