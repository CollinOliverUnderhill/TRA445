function para = define_parameters_var()
    % define_parameters_var:
    %   sets up the SPM+Thermal parameters, including
    %   - Arrhenius for diffusion & reaction rate
    %   - stoichiometry bounds
    %   - geometry & thermal constants
    %   - etc.

    % Basic constants
    para.R       = 8.314;     % J/(mol.K)
    para.F       = 96487;     % C/mol
    para.T_ref   = 298.15;    % K (25°C)
    para.C_nom   = 4.2;       % [Ah] nominal capacity

    % geometry
    para.As      = 0.0982;    % electrode area [m^2]
    para.thick1  = 70e-6;     % anode thickness
    para.thick2  = 25e-6;     % separator thickness
    para.thick3  = 70e-6;     % cathode thickness

    % solid-phase (single particle)
    para.Rs1     = 12.5e-6;   % anode particle radius
    para.Rs3     = 8.5e-6;    % cathode particle radius
    para.cs1_max = 30556;
    para.cs3_max = 51555;

    % stoichiometry limits
    para.soc0_a  = 0.0068;
    para.soc1_a  = 0.7560;
    para.soc0_c  = 0.8933;
    para.soc1_c  = 0.4650;

    % Diffusion (reference) + Activation energy
    para.Ds1_ref = 3.9e-14;
    para.Ds3_ref = 1.0e-14;
    para.Ea_Ds1  = 35e3;      % [J/mol]
    para.Ea_Ds3  = 29e3;

    % Reaction rate (reference) + Activation energy
    para.k1_ref  = 1.764e-11;  
    para.k3_ref  = 6.667e-11;
    para.Ea_k1   = 20e3;      % [J/mol]
    para.Ea_k3   = 58e3;

    % OCP polynomials (unchanged)
    para.U_ocp_anode = @(theta) ( ...
        0.7222 + 0.1387*theta + 0.0290*sqrt(theta) - 0.0172./theta + ...
        0.0019./(theta.^(1.5)) + 0.2808.*exp(0.90 - 15.*theta) - ...
        0.7984.*exp(0.4465.*theta - 0.4108) );

    para.U_ocp_cathode = @(theta) ( ...
       (-4.656 + 88.669.*theta.^2 - 401.119.*theta.^4 + 342.909.*theta.^6 ...
        - 462.471.*theta.^8 + 433.434.*theta.^10) ./ ...
       (-1 + 18.933.*theta.^2 - 79.532.*theta.^4 + 37.311.*theta.^6 ...
        - 73.083.*theta.^8 + 95.96.*theta.^10) );

    % (Optional) if you want dU/dT polynomials for reversible heat:
    para.dU_an = @(theta) 0;  % for demonstration
    para.dU_ca = @(theta) 0;  % can be replaced with actual fitted polynomials

    % Thermal parameters
    para.rho   = 1626;        % [kg/m^3]
    para.Cp    = 750;         % [J/(kg.K)]
    para.T_amb = 298.15;      % 25°C in Kelvin
    para.h     = 10;          % [W/(m^2.K)]
    para.Rc    = 0.005;       % ohmic conduction (m^2 ohm), can be tuned

    % geometry for cooling
    para.diam   = 18e-3;
    para.height = 65e-3;
    para.Vc     = pi*(para.diam/2)^2 * para.height;

    % discretization
    para.Nr_anode   = 30;
    para.Nr_cathode = 30;

    % average electrolyte concentration (assume constant for SPM)
    para.ce_avg = 1000;  
end
