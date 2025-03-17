clear; clc;
% Load data
SOC_1 = load('Q3aSOC1.mat');
SOC_4 = load('Q3aSOC4.mat');
SOC_6 = load('Q3aSOC6.mat');
initial_data = load('Q3aSOCinitial.mat');

% 提取变量
t = initial_data.t;
SOC_measure = initial_data.SOC_measure;
SOC_model = initial_data.SOC_model;
SOC_std_1 = SOC_1.SOC_model;
SOC_initial_4 = SOC_4.SOC_model;
SOC_initial_6 = SOC_6.SOC_model;

% 绘制图形
figure;
hold on;
h(1) = plot(t, SOC_measure, 'LineWidth', 2);
h(2) = plot(t, SOC_std_1, 'LineWidth', 2);
h(3) = plot(t, SOC_initial_4, 'LineWidth', 2);
h(4) = plot(t, SOC_initial_6, 'LineWidth', 2);
h(5) = plot(t, SOC_model, '--', 'LineWidth', 2);
hold off;
grid on;

% 设置图形属性
xlabel('Time (s)', 'FontSize', 18);
ylabel('SOC (%)', 'FontSize', 18);
title('SOC initial difference', 'FontSize', 18);
legend(h, {'initial SOC = 80%', 'initial SOC = 100%', 'initial SOC = 40%', 'initial SOC = 60%', 'Estimated SOC'}, 'FontSize', 15, 'Location', 'best');

set(gcf, 'Color', 'w');
set(gca, 'FontSize', 15);
%%
clear; clc;
% Load data
SOC_1 = load('Q3bmea1.mat');
SOC_10 = load('Q3bmea10.mat');
initial_data = load('Q3aSOCinitial.mat');

t = initial_data.t;
SOC_measure = initial_data.SOC_measure;
SOC_model = initial_data.SOC_model;
SOC_std_1 = SOC_1.SOC_model;
SOC_std_10 = SOC_10.SOC_model;


figure;
hold on;
h(1) = plot(t, SOC_measure, 'LineWidth', 2);
h(2) = plot(t, SOC_std_1, 'LineWidth', 2);
h(3) = plot(t, SOC_std_10, 'LineWidth', 2);
h(4) = plot(t, SOC_model, '--', 'LineWidth', 2);
hold off;
grid on;

xlabel('Time (s)', 'FontSize', 18);
ylabel('SOC (%)', 'FontSize', 18);
title('effect of SOC noise mean', 'FontSize', 18);
legend(h, {'SOC noise mean = 0', 'SOC noise mean = 1', 'SOC noise mean = 10', 'Estimated SOC'}, 'FontSize', 15, 'Location', 'best');

set(gcf, 'Color', 'w');
set(gca, 'FontSize', 15);

%%
clear; clc;
% Load data
SOC_1 = load('Q3bstdm10.mat');
SOC_10 = load('Q3bstd10.mat');
initial_data = load('Q3aSOCinitial.mat');

t = initial_data.t;
SOC_measure = initial_data.SOC_measure;
SOC_model = initial_data.SOC_model;
SOC_std_1 = SOC_1.SOC_model;
SOC_std_10 = SOC_10.SOC_model;


figure;
hold on;
h(1) = plot(t, SOC_measure, 'LineWidth', 2);
h(2) = plot(t, SOC_std_1, 'LineWidth', 2);
h(3) = plot(t, SOC_std_10, 'LineWidth', 2);
h(4) = plot(t, SOC_model, '--', 'LineWidth', 2);
hold off;
grid on;

xlabel('Time (s)', 'FontSize', 18);
ylabel('SOC (%)', 'FontSize', 18);
title('Effect of SOC noise std', 'FontSize', 18);
legend(h, {'SOC noise std = 0', 'SOC noise std = -10', 'SOC noise std = 10', 'Estimated SOC'}, 'FontSize', 15, 'Location', 'best');

set(gcf, 'Color', 'w');
set(gca, 'FontSize', 15);