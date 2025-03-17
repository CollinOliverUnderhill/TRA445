%% task2_relaxation.m
% 独立运行的 Task 2 脚本：浓度松弛模拟
% 本脚本假设参数结构体 para 已在工作区中存在，且相关函数（spm_ode_system、diffusion_sphere）已定义。
%
% 工作流程：
% 1. 利用已有的充电结束状态（c_n_end, c_p_end）作为松弛初始条件
% 2. 模拟无外加电流情况下的松弛过程（例如 1 小时）
% 3. 绘制松弛过程中表面浓度的变化和整个径向浓度分布的动画
% 4. （可选）计算松弛结束时的平均浓度

%% 请确保参数结构体 para 已经存在于工作区中，
% 如果没有，可以调用你的 mlx 文件或加载参数文件：
%   load('para.mat');  % 例如

%% 
function dcdt = diffusion_sphere(c, R, D_ref, N)
    % c: vector of length N+1
    % R: particle radius
    % D_ref: diffusion coefficient
    % N: # of sub-intervals in radial discretization
    
    dr = R / N;
    r  = linspace(0, R, N+1)';
    
    dcdt = zeros(size(c));
    
    % Interior nodes
    for i = 2:N
        dcdt(i) = D_ref * ( ...
            (c(i+1) - 2*c(i) + c(i-1)) / dr^2 + ...
            (2/r(i)) * (c(i+1) - c(i-1)) / (2*dr) );
    end
    dcdt(N+1)=D_ref*((2*c(N)-2*c(N+1))/dr^2 );
    % Boundary at r=0 (due to symmetry)
    dcdt(1) = 2 * D_ref * ...
    (c(2) - c(1)) / dr^2;  % or a specialized second-order approach
end

%% 
function dcdt = spm_ode_system(~, x, para)
    % Unpack states
    Nn = para.Nr_anode;
    Np = para.Nr_cathode;
    
    c_n = x(1:Nn+1);                       % anode concentration profile
    c_p = x(Nn+2 : Nn+1 + (Np+1));         % cathode concentration profile
    
    % Initialise dxdt
    dcdt = zeros(size(x));
    
    % Compute diffusion in anode
    dc_n = diffusion_sphere(c_n, para.Rs1, para.Ds1_ref, Nn);
    
    % Compute diffusion in cathode
    dc_p = diffusion_sphere(c_p, para.Rs3, para.Ds3_ref, Np);
    
    % Flux boundary condition from the applied current
    
    % Negative sign for anode if current is extraction of Li
    flux_n = para.I_app / para.F / para.thick1;  % [mol/(m^2·s)]
    flux_p = -para.I_app / para.F / para.thick3;  % same magnitude, but direction is different
    
    % Apply boundary condition at the outer radius:
    dc_n(end) = dc_n(end) - flux_n;  % anode
    dc_p(end) = dc_p(end) - flux_p;  % cathode
    
    % Combine
    dcdt(1:Nn+1) = dc_n;
    dcdt(Nn+2 : Nn+1+(Np+1)) = dc_p;
end
%% Step 0: 检查充电结束时的浓度状态是否存在
if ~exist('c_n_end','var') || ~exist('c_p_end','var')
    error('请先运行充电过程，获得变量 c_n_end 和 c_p_end 作为松弛初始状态。');
end

%% Step 1: 设置松弛初始条件与时间区间
x0_relax = [c_n_end, c_p_end];  % 充电结束状态作为松弛初始状态

% 松弛过程中无外加电流
para.I_app = 0;
para.time_span = [0, 3600];  % 模拟松弛1小时

%% Step 2: 求解松弛过程的ODE系统
opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
[t_relax, x_relax] = ode15s(@(t,x) spm_ode_system(t,x,para), para.time_span, x0_relax, opts);

% 分离负极和正极浓度数据
c_n_relax = x_relax(:, 1:para.Nr_anode+1);
c_p_relax = x_relax(:, para.Nr_anode+2:end);

%% Step 3: 绘制松弛过程中表面浓度随时间的变化
figure('Name','Surface Concentration During Relaxation');
subplot(2,1,1);
plot(t_relax, c_n_relax(:, end), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Anode Surface Concentration (mol/m^3)');
title('Anode Surface Concentration During Relaxation');
grid on;

subplot(2,1,2);
plot(t_relax, c_p_relax(:, end), 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Cathode Surface Concentration (mol/m^3)');
title('Cathode Surface Concentration During Relaxation');
grid on;

%% Step 4: 动画展示松弛过程中径向浓度分布的演化
% 生成径向坐标（负极和正极）
r1 = linspace(0, para.Rs1, para.Nr_anode+1);
r3 = linspace(0, para.Rs3, para.Nr_cathode+1);

figure('Name','Relaxation Profile Evolution');
for k = 1:2:length(t_relax)
    subplot(2,1,1);
    plot(r1, c_n_relax(k,:), 'LineWidth', 2);
    xlabel('r (m)');
    ylabel('Anode Concentration (mol/m^3)');
    title(sprintf('Anode Profile at t = %.1f s', t_relax(k)));
    xlim([0 para.Rs1]); ylim([0 para.cs1_max]);
    grid on;
    
    subplot(2,1,2);
    plot(r3, c_p_relax(k,:), 'LineWidth', 2);
    xlabel('r (m)');
    ylabel('Cathode Concentration (mol/m^3)');
    title(sprintf('Cathode Profile at t = %.1f s', t_relax(k)));
    xlim([0 para.Rs3]); ylim([0 para.cs3_max]);
    grid on;
    
    drawnow;
    pause(0.1);  % 根据需要调整暂停时间以控制动画速度
end

%% Step 5: （可选）计算松弛结束时的平均浓度
% 平均浓度计算公式： avg = (3/R^3)*∫0^R c(r)*r^2 dr
c_n_final = x_relax(end, 1:para.Nr_anode+1);
c_p_final = x_relax(end, para.Nr_anode+2:end);

avg_anode = (3 / para.Rs1^3) * trapz(r1, c_n_all(end,:) .* (r1.^2));
avg_cathode = (3 / para.Rs3^3) * trapz(r3, c_p_all(end,:) .* (r3.^2));


fprintf('Relaxation Final Average Concentrations:\n');
fprintf('Anode: %.2f mol/m^3\n', avg_anode);
fprintf('Cathode: %.2f mol/m^3\n', avg_cathode);


