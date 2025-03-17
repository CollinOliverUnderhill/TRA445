
% 
% function run_spm_thermal_part3()
%     clear; close all; clc;
% 
%    para = define_parameters_var();
% 
%     % 给定热模型额外参数
%     para.Ea_Ds1 = 35e3;  % Activation energy for anode diffusion [J/mol]
%     para.Ea_Ds3 = 29e3;  % Activation energy for cathode diffusion [J/mol]
%     para.Ea_k1  = 20e3;  % Activation energy for anode reaction rate [J/mol]
%     para.Ea_k3  = 58e3;  % Activation energy for cathode reaction rate [J/mol]
% 
% 
%     % 18650 热模型参数
%     para.Cp     = 750;         % [J/kg/K]
%     para.rho    = 1626;        % [kg/m^3]
%     para.height = 65e-3;       % [m]
%     para.diam   = 18e-3;       % [m]
%     para.Vc     = pi*(para.diam/2)^2*para.height;  % 圆柱体体积
%     para.h      = 30;          % [W/m^2/K]
%     para.T_amb  = 25 + 273.15; % [K]
% 
%     % =============== 1) 定义充放电倍率、初始化结果存储 ==================
%     C_rates = [1, 4];  % 1C, 4C
%     colors  = {'b','r'};
%     legendTxt = {};
% 
%     figure('Name','Temperature vs Capacity');
% 
%     for i = 1:length(C_rates)
%         % 1. 设置倍率
%         Crate = C_rates(i);
%         % 假设放电 I_app>0, 充电 I_app<0(也可反之)
%         % 下面假设: Crate>0 => 放电
%         I_app = Crate * para.C_nom / para.As;  
%         para.I_app = I_app;  % [A/m^2 * m^2 => A]
% 
%         % 2. 设置仿真时间(理论上一小时 * 1/C_rate)
%         %    但4C时，满电量仅需15min，所以 time_span=15min=900s
%         sim_time = 3600 / Crate;  
%         para.time_span = [0 sim_time];
% 
%         % 3. 初始化浓度 + 温度
%         para.Nr_anode   = 30;
%         para.Nr_cathode = 30;
% 
%         % 如果是放电(假设 SOC=100% 开始)，则 anode 满, cathode 空
%         c_n0 = para.soc1_a * para.cs1_max;  % anode
%         c_p0 = para.soc1_c * para.cs3_max;  % cathode
% 
%         % 初始温度 (与环境相同)
%         T0 = para.T_amb;
% 
%         % 拼接初始状态
%         x0 = [repmat(c_n0, para.Nr_anode+1,1);
%               repmat(c_p0, para.Nr_cathode+1,1);
%               T0];
% 
%         % 4. 数值求解
%         opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
%         [t_sol, X_sol] = ode15s(@(t,x) spm_ode_system_thermal_var(t, x, para), ...
%                                 para.time_span, x0, opts);
% 
%         % 5. 计算容量 & 温度
%         %    容量(Ah) = (I[A]) * t[s]/3600
%         cap  = (I_app * t_sol) / 3600;  % [Ah]
%         T_all = X_sol(:, end);         % 最后一列是温度 [K]
% 
%         % 转换为摄氏度
%         T_degC = T_all - 273.15;
% 
%         % 6. 绘图: 温度 vs. 容量
%         plot(cap, T_degC, 'LineWidth',2, 'Color', colors{i}); hold on;
%         legendTxt{end+1} = [num2str(Crate) 'C'];
%     end
% 
%     xlabel('Delivered Capacity (Ah)');
%     ylabel('Cell Temperature (^{\circ}C)');
%     title('Temperature vs. Capacity for 1C and 4C Discharge');
%     legend(legendTxt, 'Location','Best');
%     grid on;
% end










function run_spm_thermal_var()
    % run_spm_thermal_var:
    %   Demonstrates an SPM + thermal model where
    %   (1) Reaction rates depend on temperature (Arrhenius)
    %   (2) Overpotentials depend on local surface concentration
    %   (3) Temperature is solved lumped, but heat generation depends on c_surf, T, I.
    %
    %   We compare e.g. 1C and 4C (discharge or charge).

    clc; close all;

    %% 1) Load parameters
    para = define_parameters_var();

    % Two C-rates for demonstration
    Crates = [1, 4];

    % Structure to store results
    results = struct('crate', {}, 'time', {}, 'tempK', {}, 'capAh', {}, 'voltage', {});

    for idx = 1:length(Crates)
        crate = Crates(idx);

        %% (A) Decide if we do discharge or charge
        % Example: discharge => negative current
        I_app_Ap = -crate * para.C_nom / para.As;
        para.I_app = I_app_Ap;

        %% (B) Set simulation time:
        if crate == 1
            para.time_span = [0, 3600];  % up to 1 hour
        else
            para.time_span = [0, 900];   % up to 15 min
        end

        %% (C) Initial conditions: assume "discharge" from full stoich
        c_n0 = para.soc1_a * para.cs1_max;  % anode is near "full Li"
        c_p0 = para.soc1_c * para.cs3_max;  % cathode is near "low Li"

        x_c_n0 = repmat(c_n0, para.Nr_anode+1, 1);
        x_c_p0 = repmat(c_p0, para.Nr_cathode+1, 1);

        % Initial temperature
        T0 = para.T_amb;

        % Combine states
        x0 = [x_c_n0; x_c_p0; T0];

        %% (D) Solve PDE + thermal ODE system
        opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
        [t_sol, x_sol] = ode15s(@(t,x) spm_ode_system_thermal_var(t, x, para), ...
                                para.time_span, x0, opts);

        %% (E) Post-processing
        T_all = x_sol(:, end);      % last state is temperature
        I_abs = abs(I_app_Ap);
        I_total = I_abs * para.As;  % A
        cap_used = (I_total .* t_sol)/3600;  % [Ah]

        % If we want the terminal voltage at each step:
        V_cell = zeros(size(t_sol));
        for k = 1:length(t_sol)
            V_cell(k) = compute_terminal_voltage_thermal_var(x_sol(k,1:end-1)', x_sol(k,end), para);
        end

        % Store
        results(idx).crate   = crate;
        results(idx).time    = t_sol;
        results(idx).tempK   = T_all;
        results(idx).capAh   = cap_used;
        results(idx).voltage = V_cell;
    end

    %% 2) Plot temperature vs. capacity
    figure('Color','w'); hold on; grid on;
    for i = 1:length(results)
        plot(results(i).capAh, results(i).tempK - 273.15, 'LineWidth', 2, ...
             'DisplayName', [num2str(results(i).crate),'C']);
    end
    xlabel('Delivered Capacity (Ah)','FontSize',14);
    ylabel('Cell Temperature (°C)','FontSize',14);
    title('Temperature vs. Capacity (SPM+Thermal with c,T-dependent params)','FontSize',16);
    legend('Location','best');
    set(gca,'FontSize',13);

    %% 3) (Optional) Plot voltage vs. capacity
    figure('Color','w'); hold on; grid on;
    for i = 1:length(results)
        plot(results(i).capAh, results(i).voltage, 'LineWidth', 2, ...
             'DisplayName', [num2str(results(i).crate),'C']);
    end
    xlabel('Delivered Capacity (Ah)','FontSize',14);
    ylabel('Voltage (V)','FontSize',14);
    title('Terminal Voltage vs. Capacity','FontSize',16);
    legend('Location','best');
    set(gca,'FontSize',13);

end