clc; clear; close all;

% 设置 Ns 和 Np 范围
maxNs = 10;
maxNp = 10;
numTrials = 100; % 每个 Ns, Np 组合运行 10 次

% 预分配存储浪费比例的矩阵
waste = zeros(maxNs, maxNp);

% 执行不同 Ns 和 Np 组合的仿真
for i = 1:maxNs
    for j = 1:maxNp
        Ns = i;  % 设置串联模块数
        Np = j;  % 设置并联电池数
        tempResults = []; % 存储多次仿真结果

        for trial = 1:numTrials
            run('simPCMtest.m'); % 运行仿真脚本
            
            % 只存入非负数结果
            if Waste_capacity >= 0
                tempResults = [tempResults, Waste_capacity]; % 累积合理的仿真结果
            end
        end
        
        % 计算均值并存入 waste 矩阵
        if ~isempty(tempResults)
            waste(i,j) = mean(tempResults);
        else
            waste(i,j) = NaN; % 若所有仿真值都无效，设为 NaN
        end
    end
end

% 绘制三维柱状图
figure;
bar3(waste); 
xlabel('Parallel Cells per Module (Np)');
ylabel('Series Modules (Ns)');
zlabel('Average Wasted Capacity (%)');
title('Battery Pack Capacity Waste for Different Configurations');
colorbar;
grid on;
set(gca, 'XTick', 1:maxNp, 'YTick', 1:maxNs);
view(135, 30); % 角度调整，方便观察

% 显示 NaN 数据（如果有）
disp('Configurations with NaN values (invalid cases due to negative waste):');
[row, col] = find(isnan(waste));
for idx = 1:length(row)
    fprintf('Ns = %d, Np = %d\n', row(idx), col(idx));
end
