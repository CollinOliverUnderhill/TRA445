function [V_cell, Q_rxn, Q_rev] = compute_terminal_voltage_thermal_var(c_full, T, para)
% compute_terminal_voltage_thermal_var:
%   c_full = [c_n(1..end), c_p(1..end)]
%   T = cell temperature
%   returns:
%       V_cell = terminal voltage
%       Q_rxn  = irreversible reaction heat
%       Q_rev  = reversible heat

    Nn = para.Nr_anode;
    Np = para.Nr_cathode;

    % 拆分负极 & 正极浓度
    c_n = c_full(1:Nn+1);
    c_p = c_full(Nn+2 : Nn+1 + (Np+1));

    % 表面浓度 (末端节点)
    c_n_surf = c_n(end);
    c_p_surf = c_p(end);

    % 化学计量比
    x1 = c_n_surf / para.cs1_max;  % anode
    y3 = c_p_surf / para.cs3_max;  % cathode

    % OCP (与温度无关的部分)
    U_n = para.U_ocp_anode(x1);
    U_p = para.U_ocp_cathode(y3);

    % Arrhenius 形式的反应速率
    k1 = para.k1_ref * exp( -para.Ea_k1 / para.R * (1/T - 1/para.T_ref) );
    k3 = para.k3_ref * exp( -para.Ea_k3 / para.R * (1/T - 1/para.T_ref) );

    % 交换电流 i_0
    c_e = para.ce_avg;  % 假设电解液浓度不变
    i_0n = k1 * para.F * (c_e^0.5) * (c_n_surf^0.5) * ((para.cs1_max - c_n_surf)^0.5);
    i_0p = k3 * para.F * (c_e^0.5) * (c_p_surf^0.5) * ((para.cs3_max - c_p_surf)^0.5);

    % 总电流
    I_tot = para.I_app * para.As;  % [A]

    % 过电位 (对称BV, alpha=0.5)
    RT_F = (para.R * T) / para.F;
    eta_n = 2 * RT_F * asinh( I_tot / (2 * i_0n + eps) );
    % 注意正极的电流方向相反
    eta_p = 2 * RT_F * asinh( -I_tot / (2 * i_0p + eps) );

    % 端电压 (不含Rc, 你也可以放在这里减去 I_tot*Rc)
    V_cell = (U_p + eta_p) - (U_n + eta_n);

    % 不可逆反应热 Q_rxn = I * (|η_n| + |η_p|)
    %   或更严格地：Q_rxn = I*(η_n + η_p) 的符号与放电/充电方向有关
    Q_rxn = abs(I_tot) * ( abs(eta_n) + abs(eta_p) );

    % 可逆热 Q_rev = I * T * dU/dT
    %   需要对负极、正极分别计算 dU/dT，再相减
    dU_n_dT = para.dU_an(x1);  % anode dU/dT
    dU_p_dT = para.dU_ca(y3);  % cathode dU/dT
    Q_rev = I_tot * T * (dU_p_dT - dU_n_dT);
end
