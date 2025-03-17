clear all;
load ECM2;

% --------------------------------------------------------------------
% 4f (改进示例): 不使用 SoC/OCV 信息，增加偏置项
% 回归模型: v(k) = theta0 + theta1*v(k-1) + theta2*i(k) + theta3*i(k-1)
% 并在 "Simple Check" 输出 alpha, R0, R1, C
% --------------------------------------------------------------------

Ts = 0.1;           % 采样时间 (s)
N  = length(Time);  % 数据长度

% 直接使用原始端电压，不减去 OCV
Y = Voltage;

% 构造回归矩阵 Phi
% 每行: [1, v(k-1), i(k), i(k-1)]
Phi = [ ones(N-1,1), ...
        Y(1:N-1), ...
        Current(2:N), ...
        Current(1:N-1) ];

% ============= RLS 参数设置 =============
lambda = 0.995;     % 遗忘因子
theta_init = [2; 0.9; 0.001; 0.0005];  % 初始参数 (theta0, alpha, R0, partial for R1)
P_init = diag([10, 1, 10, 10]);       % 初始协方差矩阵

% 预分配 theta
theta = zeros(4, N-1);
theta(:,1) = theta_init;

% 初始化 P
P = P_init;

% ============= 递归最小二乘 (RLS) =============
for k = 2:(N-1)
    eps_k = Y(k) - Phi(k,:) * theta(:, k-1);      % 预测误差
    denom = lambda + Phi(k,:) * P * Phi(k,:)';
    K = (P * Phi(k,:)') / denom;                 % 增益
    theta(:,k) = theta(:,k-1) + K * eps_k;       % 参数更新
    P = (1/lambda)*(eye(4) - K*Phi(k,:))*P;      % 协方差更新
end

% ============= 绘图：观察 4 个参数随时间变化 =============
figure;
plot(Time(1:end-1), theta, 'LineWidth', 1.2);
legend('theta0 (bias)','theta1 (\alpha)','theta2 (R0 part)','theta3 (R1 part)','Location','best');
xlabel('Time (s)');
ylabel('Parameter Estimates');
title('RLS with Bias (No SoC/OCV) + Simple Check');
grid on; axis tight;

% ============= Simple Check: 映射为 alpha, R0, R1, C =============
theta_final = theta(:,end);

% 对应前述映射关系
theta0_hat = theta_final(1);  % 偏置
alpha_hat  = theta_final(2);  % alpha
R0_hat     = theta_final(3);  % R0
R1_hat     = (theta_final(4) + R0_hat * alpha_hat) / (1 - alpha_hat);  % R1

if (alpha_hat <= 0) || (alpha_hat >= 1) || (abs(log(alpha_hat)) < 1e-6)
    C_hat = inf;
else
    C_hat = -Ts / (R1_hat * log(alpha_hat));
end

% 打印结果
fprintf('\n----- Simple Check (Bias method) -----\n');
fprintf('theta0 (bias)  = %.6f\n', theta0_hat);
fprintf('alpha_hat      = %.6f\n', alpha_hat);
fprintf('R0_hat         = %.6f ohm\n', R0_hat);
fprintf('R1_hat         = %.6f ohm\n', R1_hat);
fprintf('C_hat          = %.2f F\n', C_hat);



