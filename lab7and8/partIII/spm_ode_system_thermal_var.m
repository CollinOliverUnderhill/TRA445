function dXdt = spm_ode_system_thermal_var(~, X, para)
    % spm_ode_system_thermal_var:
    %   1) radial diffusion in anode/cathode
    %   2) lumped thermal model with T depending on Q_gen - Q_loss
    %   3) Q_gen depends on ohmic + reaction + reversible
    %
    % States X = [c_n(1..Nn+1), c_p(1..Np+1), T]^T

    Nn = para.Nr_anode;
    Np = para.Nr_cathode;

    % 解包状态
    c_n = X(1:Nn+1);
    c_p = X(Nn+2 : Nn+1 + (Np+1));
    T   = X(end);

    dXdt = zeros(size(X));

    %% (A) 计算温度依赖的扩散系数 (Arrhenius)
    Ds1 = para.Ds1_ref * exp( -para.Ea_Ds1 / para.R * (1/T - 1/para.T_ref) );
    Ds3 = para.Ds3_ref * exp( -para.Ea_Ds3 / para.R * (1/T - 1/para.T_ref) );

    % (A1) 在负极/正极计算径向扩散 (使用Ghost Node处理 r=0)
    dc_n = diffusion_sphere_withGhost_var(c_n, para.Rs1, Ds1, Nn);
    dc_p = diffusion_sphere_withGhost_var(c_p, para.Rs3, Ds3, Np);

    % (A2) 外表面通量边界 (由 I_app 决定)
    flux_n =  para.I_app / (para.F * para.thick1);  
    flux_p = -para.I_app / (para.F * para.thick3);

    dc_n(end) = dc_n(end) - flux_n;
    dc_p(end) = dc_p(end) - flux_p;

    dXdt(1:Nn+1) = dc_n;
    dXdt(Nn+2 : Nn+1 + (Np+1)) = dc_p;

    %% (B) 热力学方程 dT/dt
    % 计算总热产生 Q_gen
    [~, Q_rxn, Q_rev] = compute_terminal_voltage_thermal_var([c_n; c_p], T, para);

    % (B1) 欧姆热 (可以放在 compute_terminal_voltage_thermal_var 内或外部)
    I_total = para.I_app * para.As;  % [A]
    IR_drop = I_total * para.Rc;     % 欧姆导线/集流体的电阻
    Q_ohmic = abs(I_total * IR_drop);

    % (B2) 总热产生
    Q_gen  = Q_ohmic + Q_rxn + Q_rev;

    % (B3) 对流散热
    %   假设整个圆柱表面都以 h 进行对流散热
    %   对于18650: 表面积 = pi*D*height (略去顶部/底部或合并到 h 中)
    A_surf  = pi * para.diam * para.height; 
    Q_loss  = para.h * A_surf * (T - para.T_amb);

    % (B4) 温度微分方程
    %   m_cell = ρ * V_cell
    m_cell = para.rho * para.Vc;  
    dTdt   = (1/(m_cell * para.Cp)) * (Q_gen - Q_loss);

    % 将 dTdt 放入状态导数
    dXdt(end) = dTdt;
end
