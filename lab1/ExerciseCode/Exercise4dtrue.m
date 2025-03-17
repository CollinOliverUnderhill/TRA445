%% 清空工作区并载入数据
clear; clc;
load ECM2

%% 1) 确保 Voltage、Current、Time 为列向量
if size(Voltage,1) < size(Voltage,2)
    Voltage = Voltage';
end
if size(Current,1) < size(Current,2)
    Current = Current';
end
if size(Time,1) < size(Time,2)
    Time = Time';
end

Time_len = length(Time);

%% 2) 物理参数的目标值
R0_ref = 1.0e-3;    % 目标 R0 = 0.001 ohms
R1_ref = 1.5e-3;    % 目标 R1 = 0.0015 ohms
C_ref  = 1.0e4;     % 目标 C  = 10000 F

%% 3) 多项式阶数的离散搜索范围 (例如从 6 到 15)
polyDegrees = 6:15;

%% 4) 定义连续变量的搜索范围
% x(1): 遗忘因子 lambda, 范围 [0.85, 1]
% x(2): 初始 theta(1), 范围 [0.75, 1.15]
% x(3): 初始 theta(2), 范围 [-0.005, 0.005]
% x(4): 初始 theta(3), 范围 [-0.005, 0.003]
lb = [0.85,  0.75, -0.005, -0.005];
ub = [1,     1.15,  0.005,  0.003];
nvars = 4;

%% 5) 设置粒子群算法选项
options = optimoptions('particleswarm',...
    'Display','iter',...
    'SwarmSize',100,...   % 增大群体规模
    'MaxIterations',200); % 增加迭代次数

%% 6) 固定参数
Ts = 0.1;         % 采样时间 (s)
Q  = 60*3600;     % 电池容量 (As)
z_init = 0.5;     % 初始 SoC

%% 7) 开始离散搜索：对于每个多项式阶数，利用 PSO 优化连续变量
best_global_err = Inf;
best_global_result = struct();

for deg = polyDegrees
    % 定义目标函数，固定当前多项式阶数为 deg
    objFun = @(x) objective_func_phys_only( ...
        deg, x, R0_ref, R1_ref, C_ref, Ts, Q, z_init);
    
    % 使用粒子群算法搜索连续变量
    [x_opt, fval] = particleswarm(objFun, nvars, lb, ub, options);
    
    fprintf('Degree = %d | best x_opt = [%f, %f, %f, %f], err_phys = %f\n',...
        deg, x_opt, fval);
    
    if fval < best_global_err
        best_global_err = fval;
        best_global_result.polyDegree = deg;
        best_global_result.x_opt = x_opt;
        best_global_result.err = fval;
    end
end

%% 8) 输出最终的全局最佳结果
disp('===== Best overall across all degrees =====');
disp(best_global_result);

%% ========== 定义目标函数 (仅物理参数惩罚项) ==========
function err_total = objective_func_phys_only(...
    polyDegree, x, R0_ref, R1_ref, C_ref, Ts, Q, z_init)
    % 输入 x:
    %   x(1): lambda
    %   x(2): 初始 theta(1)
    %   x(3): 初始 theta(2)
    %   x(4): 初始 theta(3)
    lambda = x(1);
    initTheta = [x(2); x(3); x(4)];
    
    %% 获取全局数据
    global Voltage Current Time OCV;
    if isempty(Voltage) || isempty(Current) || isempty(Time) || isempty(OCV)
        Voltage = evalin('base','Voltage');
        Current = evalin('base','Current');
        Time    = evalin('base','Time');
        OCV     = evalin('base','OCV');
    end
    Time_len = length(Time);
    
    %% 1) 对 OCV 数据进行预处理：标准化 X 轴 & 去重
    X_norm = (OCV(:,1) - mean(OCV(:,1))) / std(OCV(:,1));
    [X_unique, ia] = unique(X_norm);
    Y_unique = OCV(ia,2);
    
    % 限制多项式阶数不超过唯一数据点数 - 1
    maxDegree = min(polyDegree, length(X_unique) - 1);
    Pocv = polyfit(X_unique, Y_unique, maxDegree);
    
    %% 2) 计算 SoC 和对应的 vOC（开路电压）
    z = zeros(Time_len,1);
    vOC = zeros(Time_len,1);
    z(1) = z_init;
    vOC(1) = polyval(Pocv, (z(1) - mean(OCV(:,1))) / std(OCV(:,1)));
    for k = 1:Time_len-1
        z(k+1) = z(k) + (Ts/Q)*Current(k);
        vOC(k+1) = polyval(Pocv, (z(k+1) - mean(OCV(:,1))) / std(OCV(:,1)));
    end
    
    %% 3) 构造回归数据：Y = Voltage - vOC
    Y = Voltage - vOC;
    Phi = [Y(1:Time_len-1), Current(2:Time_len), Current(1:Time_len-1)];
    
    %% 4) RLS 递归更新
    theta_hat = zeros(3, Time_len-1);
    theta_hat(:,1) = initTheta;
    P_mat = 10 * eye(3);  % 初始 P 矩阵
    
    for kk = 2:Time_len-1
        eps_k = Y(kk) - Phi(kk,:) * theta_hat(:,kk-1);
        denom = lambda + Phi(kk,:) * P_mat * Phi(kk,:)';
        K = (P_mat * Phi(kk,:)') / denom;
        theta_hat(:,kk) = theta_hat(:,kk-1) + K * eps_k;
        P_mat = (1/lambda) * (eye(3) - K * Phi(kk,:)) * P_mat;
    end
    theta_final = theta_hat(:, end);
    
    %% 5) 映射到物理参数
    % 定义: alpha_hat = theta_final(1)
    %      R0_hat    = theta_final(2)
    %      R1_hat    = (theta_final(3) + R0_hat * alpha_hat) / (1 - alpha_hat)
    %      C_hat     = -Ts / (R1_hat * log(alpha_hat))
    alpha_hat = theta_final(1);
    R0_hat = theta_final(2);
    R1_hat = (theta_final(3) + R0_hat * alpha_hat) / (1 - alpha_hat);
    
    if (alpha_hat <= 0) || (alpha_hat >= 1) || (abs(log(alpha_hat)) < 1e-6)
        C_hat = 1e9;  % 给极大惩罚
    else
        C_hat = -Ts / (R1_hat * log(alpha_hat));
    end
    
    %% 6) 计算物理参数误差 (目标函数)
    err_phys = abs(R0_hat - R0_ref) + abs(R1_hat - R1_ref) + abs(C_hat - C_ref);
    
    % 目标函数即物理参数误差（也可以加权调整）
    err_total = err_phys;
end
