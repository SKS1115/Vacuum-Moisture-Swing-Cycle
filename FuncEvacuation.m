
%% function to solve Temperature, Pressure and Mole fraction rate during Evacuation Stage
% NB: This function does  include half cells.

function derivatives = FuncEvacuation(t_evac,EVS,dz,N,L,ColumnParam,CondParam,isotherm_parameters)

% Extract Input Parameters
    r_in    = ColumnParam(3);
    r_out   = ColumnParam(4);
    eb      = ColumnParam(5);
    ep      = ColumnParam(6);
    rp      = ColumnParam(7);
    tau     = ColumnParam(8);
    rho_s   = CondParam(1);
    rho_w   = CondParam(2);
    Cpg     = CondParam(3);
    Cpa_CO2 = CondParam(4);
    Cps     = CondParam(5);
    Cpw     = CondParam(6);
    mu      = CondParam(7);
    Dm      = CondParam(8);
    Kz      = CondParam(10);
    Kw      = CondParam(11);
    h_in    = CondParam(12);
    h_out   = CondParam(13);
    R       = CondParam(14);
    lambda  = CondParam(16);
    dp      = CondParam(17);
    Cpa_H2O = CondParam(19);
    Pi      = CondParam(22);
    Ta      = CondParam(25);
    kCO2_des = CondParam(31);
    kH2O_des = CondParam(33); 
    MW_air   = CondParam(34);
    MW_water = CondParam(35); 

% Extract Relevant Isotherm Parameters
    dUb_CO2     = isotherm_parameters(4);
    dUb_H20     = isotherm_parameters(10);

% Calculate Constant Inputs
    DL         = 0.7*Dm;
    viscous_term  = (150*(1-eb)^2*mu)/(eb^2*dp^2) ; 
    e_ratio    = (1-eb)/eb;
    Dp         = Dm/tau;            
    mdiffusion = (15*ep*Dp)/rp^2;
    Tw1        = 1/(rho_w*Cpw);
    Tw2        = (2*r_in*h_in)/(r_out^2-r_in^2);
    Tw3        = (2*r_out*h_out)/(r_out^2-r_in^2);
    Tw4        = (2*h_in)/(eb*r_in);
    dz2        = dz/2;

% Initialize all main parameters
    % Temporal derivatives
    derivatives   = zeros(length(EVS),1);
    dyAdt   = zeros(N+2,1);
    dyAdt1  = zeros(N+2,1);
    dyAdt2  = zeros(N+2,1);
    dyAdt3  = zeros(N+2,1);
    dyAdt4  = zeros(N+2,1);
    dyBdt   = zeros(N+2,1);
    dyBdt1  = zeros(N+2,1);
    dyBdt2  = zeros(N+2,1);
    dyBdt3  = zeros(N+2,1);
    dyBdt4  = zeros(N+2,1);
    dpdt    = zeros(N+2,1);
    dpdt1   = zeros(N+2,1);
    dpdt2   = zeros(N+2,1);
    dpdt3   = zeros(N+2,1);
    dTdt    = zeros(N+2,1);
    dTdt1   = zeros(N+2,1);
    dTdt2   = zeros(N+2,1);
    dTdt3   = zeros(N+2,1);
    dTdt4   = zeros(N+2,1);
    dTdt5   = zeros(N+2,1);
    dTdt6   = zeros(N+2,1);
    dTWdt   = zeros(N+2,1);
    dqdt1   = zeros(N+2,1);
    dqdt2   = zeros(N+2,1);
    
% Spatial derivatives
    u          = zeros(N+1,1); 
    rho_gh     = zeros(N+1,1);
    kinetic_term = zeros(N+1,1); 
    LDF1       = zeros(N+2,1);
    LDF2       = zeros(N+2,1);   
    dpdz       = zeros(N+2,1);
    dpdzh      = zeros(N+1,1);
    dTdz       = zeros(N+2,1);
    dyAdz      = zeros(N+2,1);
    dyAdzh     = zeros(N+2,1);
    d2yAdz2    = zeros(N+2,1);
    dyBdz      = zeros(N+2,1);
    dyBdzh     = zeros(N+2,1);
    d2yBdz2    = zeros(N+2,1);
    d2Tdz2     = zeros(N+2,1);
    dTW2dz2    = zeros(N+2,1);
    
% Define values to output matrix
    YA  = EVS(1:N+2);
    YB  = EVS(N+3:2*N+4);
    P   = EVS(2*N+5:3*N+6);
    q1  = EVS(3*N+7:4*N+8);
    q2  = EVS(4*N+9:5*N+10);
    T   = EVS(5*N+11:6*N+12);
    TW  = EVS(6*N+13:7*N+14);

% Boundary Conditions
    %%% Inlet 
    YA(1)     = YA(2);    
    YB(1)     = YB(2);    
    P(1)      = P(2);     
    T(1)      = T(2);       
    TW(1)     = Ta;


    %%% Outlet
    YA(end)   = YA(end-1);
    YB(end)   = YB(end-1);
    T(end)    = T(end-1);    
    TW(end)   = Ta;   
    if P(end) > P(end-1)
       P(end) = P(end-1);
    end

    LDF1(:) = kCO2_des;
    LDF2(:) = kH2O_des;

%% Internal Nodes Spatial Derivatives
% 1) Pressure at the center of the volumes and walls
    Ph           = WENO(P)                ;  % solves the half cell approximation
    dpdz(2:N+1)   = (Ph(2:N+1)-Ph(1:N))./dz ;  % pressure gradient at the walls   
    dpdzh(1)      =  2*(P(2)-P(1))/dz      ;  
    dpdzh(2:N)    = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(N+1)    =  2*(P(N+2)-P(N+1))/dz  ;

% 2) Mole fraction at the center of the volumes
    yAh          = WENO(YA)                 ;
    dyAdz(2:N+1) = (yAh(2:N+1)-yAh(1:N))/dz ;
    yBh          = WENO(YB)                 ;
    dyBdz(2:N+1) = (yBh(2:N+1)-yBh(1:N))/dz ;

% 3) Temperature at the center of the volumes
    Th          = WENO(T)      ;
    dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;

% 4) 2nd derivatives
   % CO2 Mole Fraction
    d2yAdz2(2)   = (YA(3)-YA(2))/dz/dz                  ;
    d2yAdz2(3:N) = (YA(4:N+1)+YA(2:N-1)-2*YA(3:N))/dz/dz;
    d2yAdz2(N+1) = (YA(N)-YA(N+1))/dz/dz                ;

   % H2O Mole Fraction
    d2yBdz2(2)   = (YB(3)-YB(2))/dz/dz                  ;
    d2yBdz2(3:N) = (YB(4:N+1)+YB(2:N-1)-2*YB(3:N))/dz/dz;
    d2yBdz2(N+1) = (YB(N)-YB(N+1))/dz/dz                ;
    
   % Temperature
    d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;

 % 5) Velocities at the boundary interface
   rho_gh(1:N+1)   = Ph(1:N+1)./(R.*Th(1:N+1));
   kinetic_term(1:N+1) = (1.75.*(1-eb).*rho_gh(1:N+1).*(MW_air+(MW_water-MW_air).*yBh(1:N+1)))./(eb*dp) ;    
   discriminant = abs(viscous_term^2 + (4*kinetic_term(1:N+1).*abs(dpdzh(1:N+1)))) ;
   u(1:N+1)       = -sign(dpdzh(1:N+1)).*(-viscous_term + sqrt(discriminant))./(2.*kinetic_term(1:N+1)) ; 

%% Internal Nodes Temporal Derivatives

% 1) Adsorbed Mass Balance (molar loading for CO2 and N2)*
%       Linear Driving Force (LDF) is used to calculate the uptake rate of
%       the gas. The LDF mass transfer coefficients are assumed to be
%       constant over the pressure and temperature range of importance but can be varied.

%    1.1 Calculate equilibrium molar loading        
        qstar = IsothermMS(YA, YB, P, T, CondParam, isotherm_parameters);
        qstar1 = qstar(:, 1);
        qstar2 = qstar(:, 2);        
        
%    1.2 Calculate the uptake rate  
        dqdt1(2:N+1) = LDF1(2:N+1).*(qstar1(2:N+1)-q1(2:N+1));
        dqdt2(2:N+1) = LDF2(2:N+1).*(qstar2(2:N+1)-q2(2:N+1));

    
% 2) Column Energy Balance
%       For the energy balance in the column, it is assumed that:
%       * There is thermal equilibrium between the solid and gas phase
%       * Conduction occurs through both the solid and gas phase 
%       * Advection only occurs in the gas phase 

%     2.1 Calculate the Wall Temperature 
        dTW2dz2(2:N+1) = (TW(3:N+2) - 2.*TW(2:N+1) + TW(1:N))/dz/dz;  % 2nd order derivative central difference
        dTWdt(2:N+1) = (Kw.*Tw1.*dTW2dz2(2:N+1)) + (Tw2.*Tw1.*(T(2:N+1)-TW(2:N+1))) - (Tw3.*Tw1.*(TW(2:N+1)-Ta));
    
%     2.2 Calculate the Total Temperature   
        % sink        = (rho_s*C_ps) + C_pa.*(q1(2:N+1)+q2(2:N+1));        
         sink        = (rho_s*Cps) + (Cpa_CO2.*q1(2:N+1)) + (Cpa_H2O.*q2(2:N+1)); 

        Pv  = Ph(1:N+1).*u(1:N+1);
        PvT = Ph(1:N+1).*u(1:N+1)./Th(1:N+1);
       
        dHCO2    = dUb_CO2; % - (R.*(T(2:N+1)-T(1:N))); 
        dHH20    = dUb_H20;% - (R.*(T(2:N+1)-T(1:N)));
            
        dTdt1(2:N+1) = (Kz/(1-eb)./sink).*d2Tdz2(2:N+1)                     ; % Temperature change due to conduction
        dTdt2(2:N+1) = (Cpg/R/dz/e_ratio./sink).*(Pv(2:N+1) - Pv(1:N))     ; % Temperature change due to advection
        % dTdt3(2:N+1) = (C_pa./sink).*T(2:N+1).*(dqdt1(2:N+1)+dqdt2(2:N+1))  ; % Temperature change due to reaction
        dTdt3(2:N+1) = (1./sink).*T(2:N+1).*(Cpa_CO2.*dqdt1(2:N+1) + Cpa_H2O.*dqdt2(2:N+1))  ; % Temperature change due to reaction
        dTdt4(2:N+1) = (Cpg/R/e_ratio./sink).*dpdt(2:N+1)                  ; % Temperature change due to pressure
        dTdt5(2:N+1) = ((-dHCO2./sink).*dqdt1(2:N+1)) + ((-dHH20./sink).*dqdt2(2:N+1)) ;% Temperature change due to adsorption/desorption enthalpy
        dTdt6(2:N+1) = (Tw4./(e_ratio.*sink)).*(T(2:N+1)-TW(2:N+1))         ; % Temperature change due to the wall

        dTdt(2:N+1) = dTdt1(2:N+1) - dTdt2(2:N+1) - dTdt3(2:N+1) - dTdt4(2:N+1) + dTdt5(2:N+1) - dTdt6(2:N+1);
       
% 3) Total mass balance
        
        dpdt1(2:N+1)  = -T(2:N+1).*(PvT(2:N+1) - PvT(1:N))./dz            ;  % change in pressure due to advection
        dpdt2(2:N+1)  = -R.*T(2:N+1).*e_ratio.*(dqdt1(2:N+1)+dqdt2(2:N+1));  % change in pressure due to adsorption/desorption
        dpdt3(2:N+1)  = (P(2:N+1).*dTdt(2:N+1))./T(2:N+1)                ;  % change in pressure due to temperature

        dpdt(2:N+1)   = dpdt1(2:N+1) + dpdt2(2:N+1) + dpdt3(2:N+1);  
    
    
% 4)  Calculate CO2 Mole fraction
        
        PT      = Ph(1:N+1)./Th(1:N+1);                                   % Pressure divided by Temperature at all wall nodes
        TP      = T(1:N+1)./P(1:N+1);                                     % Temperature divided by Pressure at all cehter nodes
        YA_PVT  = (yAh(1:N+1).*Ph(1:N+1).*u(1:N+1))./Th(1:N+1);              % simplification of (Y*P*V)/T, evaluated at at all nodes
                
        dyAdt1(2:N+1)  = DL.*(d2yAdz2(2:N+1)+(dyAdz(2:N+1).*dpdz(2:N+1)./P(2:N+1))-(dyAdz(2:N+1).*dTdz(2:N+1)./T(2:N+1)))                        ;
        dyAdt2(2:N+1)  = TP(2:N+1).*(1/dz).*(YA_PVT(2:N+1) - YA_PVT(1:N));
        dyAdt3(2:N+1)  = R.*TP(2:N+1).*e_ratio.*dqdt1(2:N+1);
        dyAdt4(2:N+1)  = -(YA(2:N+1).*dpdt(2:N+1)./P(2:N+1)) + (YA(2:N+1).*dTdt(2:N+1)./T(2:N+1));
    
        dyAdt(2:N+1)   = dyAdt1(2:N+1) - dyAdt2(2:N+1) - dyAdt3(2:N+1) + dyAdt4(2:N+1); 
        
% 5) Calculate H2O Mole fraction 
        YB_PVT  = (yBh(1:N+1).*Ph(1:N+1).*u(1:N+1))./Th(1:N+1);              % simplification of (Y*P*V)/T, evaluated at at all nodes
                       
        dyBdt1(2:N+1)  = DL*(d2yBdz2(2:N+1)+(dyBdz(2:N+1).*dpdz(2:N+1)./P(2:N+1))-(dyBdz(2:N+1).*dTdz(2:N+1)./T(2:N+1)))                        ;
        dyBdt2(2:N+1)  = TP(2:N+1).*(1/dz).*(YB_PVT(2:N+1) - YB_PVT(1:N));
        dyBdt3(2:N+1)  = R.*TP(2:N+1).*e_ratio.*dqdt2(2:N+1);
        dyBdt4(2:N+1)  = -(YB(2:N+1).*dpdt(2:N+1)./P(2:N+1)) + (YB(2:N+1).*dTdt(2:N+1)./T(2:N+1));
    
        dyBdt(2:N+1)   = dyBdt1(2:N+1) - dyBdt2(2:N+1) - dyBdt3(2:N+1) + dyBdt4(2:N+1); 
        
        
   % Boundary Derivatives 
    dyAdt(1)   = dyAdt(2);
    dyAdt(end) = dyAdt(end-1);
    dyBdt(1)   = dyBdt(2);
    dyBdt(end) = dyBdt(end-1);
    dpdt(1)    = dpdt(2);
    % dpdt(end) = -lambda*(PH-Pi)*exp(-lambda*t_evac);
    dpdt(end) = (-lambda)*(P(end)-Pi);
    dqdt1(1)   = 0;
    dqdt1(end) = 0;
    dqdt2(1)   = 0;
    dqdt2(end) = 0;
    dTdt(1)    = dTdt(2);
    dTdt(end)  = dTdt(end-1) ;
    dTWdt(1)   = 0;
    dTWdt(end) = 0;

 % Export Derivatives to output  
    derivatives(1:N+2)         = dyAdt(1:N+2)  ; 
    derivatives(N+3:2*N+4)     = dyBdt(1:N+2)  ;
    derivatives(2*N+5:3*N+6)   = dpdt(1:N+2)   ;
    derivatives(3*N+7:4*N+8)   = dqdt1(1:N+2)  ;
    derivatives(4*N+9:5*N+10)  = dqdt2(1:N+2)  ;
    derivatives(5*N+11:6*N+12) = dTdt(1:N+2)   ;
    derivatives(6*N+13:7*N+14) = dTWdt(1:N+2)  ;

end