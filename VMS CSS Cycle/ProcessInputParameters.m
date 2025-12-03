

function InputParameters = ProcessInputParameters(N,RH,L)

% 1) Define operating bed parameters
% 2) Define flue gas parameters and constants
% 3) Define adsorbent parameters 
% 4) Define isotherm parameters by calling Isothermdata.mat file you loaded
% in PSACycle
% 5) Distribute the values defined above to new variables e.g. Param(1) =
% R (gas constant)
% 6) Define the time for each cycle. (NA for one cycle)
% 7) Create new variable with all Isotherm Parameters. 
% 8) Combine all lists into one variable of cells that can be used in other
% scripts. Call this variable InputParams.

%``````````````````````````Column Properties`````````````````````````````
    r_in    = 0.040     ;            % Inner column radius (units:m)
    r_out   = 0.041     ;            % Outer column radius (units:m)
    eb      = 0.35      ;            % bed voidage voidage (units: dimensionless)commandhistory
    ep      = 0.238     ;            % particle voidage (units: dimensionless)
    rp      = 4e-4      ;            % IRA 900 particle radius (units: m)
    tau     = 3      ;            % adsobent tortuosity (units: dimensionless);

%`````````````````````````` Properties and Constants````````````````````````
    rho_s    = 1060      ;        % IRA 900 adsorbent density {units: kg/m3}
    rho_w    = 7800      ;        % column wall density {units: kg/m3}
    C_pg     = 29.17     ;        % specific heat capacity of gas phase {units: J/mol*K}
    C_pa_CO2 = 37.0221   ;        % CO2 specific heat capacity of the sorbed phase {units: J/mol*K}
    C_pa_H2O = 75        ;        % H2O specific heat capacity of the sorbed phase {units: J/mol*K}
    C_ps     = 1580      ;        % specific heat capacity of adsorbent {units: J/kg*K}
    C_pw     = 502       ;        % specific heat capacity of column wall {units: J/kg*K}. Corelates to stainless steel 304
    mu       = 1.72e-5   ;        % fluid viscosity {units: kg/ms}
    Dm       = 2.05e-5   ;        % effective macropore diffusivity {units: m2/s}
    gamma    = 1.4       ;        % adiabatic constant {units: dimensionless}
    Kz       = 0.09      ;        % effective gas thermal conductivity {units: J/m*K*s}
    Kw       = 16        ;        % thermal conductivity of column wall {J/m*K*s}
    h_in     = 0         ;        % inside heat transfer coefficient {J/m^2*K*s}
    h_out    = 0         ;        % outside heat transfer coefficient {J/m^2*K*s}
    R        = 8.314     ;        % universal gas constant {units: m3*Pa/mol*K}.
    eta      = 0.72      ;        % compression/evacuation efficiency (units: dimensionless)
    lambda   = 0.5       ;        % parameter used to dertemine rate of pressure change {units: 1/s}
    dp       = rp*2      ;        % diameter of the adsorbent particle {units: m}
    tho      = 10e-10    ;        % constant used in van leer expression 
    MW_air   = 0.02896   ;
    MW_water = 0.018015  ; 

%``````````````````````````Operating Conditons for case study```````````````
    Tsat     = 300     ;          % feed temperature of water during vapor stripping {units: K}
    Tatm     = 298.15  ;          % Ambient temperature {units: K}
    PH       = 1e05    ;          % High Pressure {units: Pa}
    Pi       = 500     ;          % Intermediate Pressure. {units: Pa}
    PL       = 1000    ;          % Low Pressure {units: Pa}. 
    vair     = 1       ;          % feed air velocity {units: m/s}
    vvs      = 1       ;          % vapor stripping vacuum velocity {units: m/s}

    PH2O_sat = PressureVaporSat(Tatm) ; % H2O saturation pressure at given T
    yA0      = 4e-04                  ; % feed CO2 concentration 400ppm 
    yB0      = (RH*PH2O_sat)/(100*PH) ; % feed H2O concentration 
    yC0      = 1-(yA0+yB0);            % feed mole fraction N2 

    %%% 0.5x Baseline Kinetics
       % kCO2_ads = 2.5e-4  ; 
       % kCO2_des = 2.0e-4  ; 
       % kH2O_ads = 2.25e-4 ;
       % kH2O_des = 7.0e-4  ;

    %%% Baseline Kinetics
       kCO2_ads = 5.11e-4  ; 
       kCO2_des = 4.04e-4  ; 
       kH2O_ads = 4.52e-4  ;
       kH2O_des = 1.38e-3  ;

    %%% 2x Baseline Kinetics
       % kCO2_ads = 1.0e-3  ; 
       % kCO2_des = 8.0e-4  ; 
       % kH2O_ads = 9.0e-4  ;
       % kH2O_des = 2.8e-3  ;


%```````````````````````Isotherm Parameters for CO2````````````````````` 
    IEC     = 1.91          ;       % Ion Exchange Capacity of CO2 [mol/kg]
    kms     = 0.00061236    ;       % equilibrium constant [kg/molPa]
    n       = 4.0991        ;       % exponential factor on H2O [-]
    dUb_CO2 = -70000        ;       % heat of adsorption of CO2 [J/mol]

%```````````````````````Isotherm Parameters for H20````````````````````` 
    qm_H2O  = 11.4360   ;         % H20 lodading capacity monolayer for IRA900 [mol/kg]
    E_H2O   = 38.6      ;         % Experimentally fitted constant used in 1st layer heat of asdsorption [J/mol]
    F_H2O   = 59338     ;         % Experimentally fitted constant used in 1st layer heat of asdsorption [K^-1] 
    G_H2O   = 50.4      ;         % Experimentally fitted constant used in 1st layer heat of asdsorption [J/mol]
    H_H2O   = 58322     ;         % Experimentally fitted constant used in 1st layer heat of asdsorption [J/mol/K]
    dUb_H2O = -46000    ;         % heat of adsorption of H20 [J/mol]
 
%``````````````````````` Times (seconds) `````````````````````
    t_pz  = 10  ;          % Pressurization duration
    t_cs  = 12000 ;          % CO2 Sorption duration
    t_ev  = 30  ;          % Air Evacuation duration
    t_vs  = 24000 ;          % VaporStripping duration
    t_des = 200 ;          % Desorption duration

%% Distribute values to necessary variables
    ColumnParam      = zeros(8,1);
    ColumnParam(1)   = N;
    ColumnParam(2)   = L;
    ColumnParam(3)   = r_in;
    ColumnParam(4)   = r_out;
    ColumnParam(5)   = eb;
    ColumnParam(6)   = ep;
    ColumnParam(7)   = rp;
    ColumnParam(8)   = tau;
    
    CondParam        = zeros(28,1);
    CondParam(1)     = rho_s;
    CondParam(2)     = rho_w;
    CondParam(3)     = C_pg;
    CondParam(4)     = C_pa_CO2;
    CondParam(5)     = C_ps;
    CondParam(6)     = C_pw;
    CondParam(7)     = mu;
    CondParam(8)     = Dm;
    CondParam(9)     = gamma;
    CondParam(10)    = Kz;
    CondParam(11)    = Kw;
    CondParam(12)    = h_in;
    CondParam(13)    = h_out;
    CondParam(14)    = R;
    CondParam(15)    = eta;
    CondParam(16)    = lambda;
    CondParam(17)    = dp;
    CondParam(18)    = tho;
    CondParam(19)    = C_pa_H2O;
    CondParam(20)    = PH;
    CondParam(21)    = PL;
    CondParam(22)    = Pi;
    CondParam(23)    = vair;
    CondParam(24)    = Tsat;
    CondParam(25)    = Tatm;
    CondParam(26)    = yA0;
    CondParam(27)    = yB0;
    CondParam(28)    = yC0;
    CondParam(29)    = vvs;
    CondParam(30)    = kCO2_ads;
    CondParam(31)    = kCO2_des;
    CondParam(32)    = kH2O_ads;
    CondParam(33)    = kH2O_des;
    CondParam(34)    = MW_air;
    CondParam(35)    = MW_water; 
    
    isotherm_parameters(1)  = IEC;
    isotherm_parameters(2)  = kms;
    isotherm_parameters(3)  = n;
    isotherm_parameters(4)  = dUb_CO2;
    isotherm_parameters(5)  = qm_H2O;
    isotherm_parameters(6)  = E_H2O;
    isotherm_parameters(7)  = F_H2O;
    isotherm_parameters(8)  = G_H2O;
    isotherm_parameters(9)  = H_H2O;
    isotherm_parameters(10) = dUb_H2O;
  
    Times = [t_pz; t_cs; t_ev; t_vs; t_des];

    InputParameters{1} = ColumnParam;
    InputParameters{2} = CondParam;
    InputParameters{3} = isotherm_parameters;
    InputParameters{4} = Times;

end 
