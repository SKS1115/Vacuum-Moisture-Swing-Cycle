
%% function to solve Temperature, Pressure and Mole fraction rate during Evacuation Stage
% NB: This function does include half cells.

function derivatives = FuncDesorption(t_des,DES,x40,dz,N,L,ColumnParam,CondParam,isotherm_parameters)

% Extract Input Parameters
    r_in        = ColumnParam(3); 
    r_out       = ColumnParam(4); 
    eb          = ColumnParam(5); 
    rho_s       = CondParam(1); 
    rho_w       = CondParam(2); 
    C_pg        = CondParam(3); 
    C_pa_CO2    = CondParam(4); 
    C_ps        = CondParam(5); 
    C_pw        = CondParam(6); 
    mu          = CondParam(7); 
    Dm          = CondParam(8); 
    Kz          = CondParam(10); 
    Kw          = CondParam(11); 
    h_in        = CondParam(12); 
    h_out       = CondParam(13); 
    R           = CondParam(14); 
    lambda      = CondParam(16); 
    dp          = CondParam(17); 
    C_pa_H2O    = CondParam(19); 
    PL          = CondParam(21); 
    Tatm        = CondParam(25); 
    kCO2_des    = CondParam(31); 
    kH2O_des    = CondParam(33); 
    MW_air      = CondParam(34); 
    MW_water    = CondParam(35);

% Extract Relevant Isotherm Parameters
    dUb_CO2     = isotherm_parameters(4);
    dUb_H2O     = isotherm_parameters(10);

% Calculate Constant Inputs
    DL      = 0.7*Dm;    
    e_ratio = (1-eb)/eb;    
    Tw1     = 1/(rho_w*C_pw);
    Tw2     = (2*r_in*h_in)/(r_out^2-r_in^2);
    Tw3     = (2*r_out*h_out)/(r_out^2-r_in^2);
    Tw4     = (2*h_in)/(eb*r_in);
    beta    = (150*(1-eb)^2*mu)/(eb^2*dp^2) ; 

%% Initialize all main parameters
    %%% Temporal derivatives
    derivatives   = zeros(length(DES),1);
    dyAdt    = zeros(N+2,1);
    dyAdt1   = zeros(N+2,1);
    dyAdt2   = zeros(N+2,1);
    dyAdt3   = zeros(N+2,1);
    dyAdt4   = zeros(N+2,1);
    dyBdt    = zeros(N+2,1);
    dyBdt1   = zeros(N+2,1);
    dyBdt2   = zeros(N+2,1);
    dyBdt3   = zeros(N+2,1);
    dyBdt4   = zeros(N+2,1);
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
    u       = zeros(N+1,1); 
    rho_g   = zeros(N+1,1);
    dpdz    = zeros(N+2,1);
    dpdzh   = zeros(N+1,1);
    dTdz    = zeros(N+2,1);
    dyAdz   = zeros(N+2,1);
    d2yAdz2 = zeros(N+2,1);
    dyBdz   = zeros(N+2,1);
    d2yBdz2 = zeros(N+2,1);
    d2Tdz2  = zeros(N+2,1);
    dTW2dz2 = zeros(N+2,1);
    alpha   = zeros(N+1,1);

    
% Define values to output matrix
    YA  = DES(1:N+2)        ; % CO2 mole fraction
    YB  = DES(N+3:2*N+4)    ; % H2O mole fraction
    P   = DES(2*N+5:3*N+6)  ; % Pressure
    q1  = DES(3*N+7:4*N+8)  ; % CO2 loading
    q2  = DES(4*N+9:5*N+10) ; % H2O loading
    T   = DES(5*N+11:6*N+12); % Temperature
    TW  = DES(6*N+13:7*N+14); % Wall Temperature

%% Boundary Conditions
    P_last   = x40(3*N+6) ; % Last pressure value from step4
    YA(1)    = YA(2);    
    YB(1)    = YB(2);
    P(1)     = P(2);    
    T(1)     = T(2);   
    TW(1)    = Tatm;    

    rho_g(1) = P(1)/(R*T(1));  
    MW_in    = MW_air + (MW_water-MW_air)*YB(1);
    alpha(1) = (1.75*(1-eb).*rho_g(1)*MW_in)./(eb*dp) ;
    dsmt1    = beta^2 - (4*alpha(1).*(P(2)-P(1))*(2/dz)) ;
    u(1)     = (-beta + sqrt(dsmt1))./(2*alpha(1)) ; 

    %%% Outlet
    YA(end)   = YA(end-1);
    YB(end)   = YB(end-1);
    T(end)    = T(end-1);
    TW(end)   = Tatm;
    if P(end) > P(end-1)
        P(end) = P(end-1);
    end
    
    MW_out = MW_air + (MW_water-MW_air)*YB(end);
    rho_g  = P(end)/(R*T(end));
    alpha(end) = (1.75*(1-eb).*rho_g*MW_out)./(eb*dp) ;
    dsmt2      = abs(beta^2 + (4*alpha(end).*abs(P(end)-P(end-1))*(2/dz))) ;
    u(end)     = -sign(P(end)-P(end-1))*(-beta + sqrt(dsmt2))./(2*alpha(end)) ; 
  
%% Internal Nodes Spatial Derivatives
% 1) Pressure at the center of the volumes and walls
    Ph           = WENO(P)                 ;  % solves the half cell approximation
    dpdz(2:N+1)  = (Ph(2:N+1)-Ph(1:N))./dz ;  % pressure gradient at the walls   
    dpdzh(1)     =  2*(P(2)-P(1))/dz      ;  
    dpdzh(2:N)   = (P(3:N+1)-P(2:N))/dz   ;
    dpdzh(N+1)   =  2*(P(N+2)-P(N+1))/dz  ;

% 2) CO2 and H2O Mole fraction at the center of the volumes
    yAh          = WENO(YA)                 ;
    dyAdz(2:N+1) = (yAh(2:N+1)-yAh(1:N))/dz ;
    yBh          = WENO(YB)                 ;
    dyBdz(2:N+1) = (yBh(2:N+1)-yBh(1:N))/dz ;

% 3) Temperature at the center of the volumes
    Th          = WENO(T)      ;
    dTdz(2:N+1) = (Th(2:N+1)-Th(1:N))/dz ;

% 4) 2nd derivatives
    %%% CO2 Mole Fraction
    d2yAdz2(2)   = (YA(3)-YA(2))/dz/dz                  ;
    d2yAdz2(3:N) = (YA(4:N+1)+YA(2:N-1)-2*YA(3:N))/dz/dz;
    d2yAdz2(N+1) = (YA(N)-YA(N+1))/dz/dz                ;

    %%% H2O Mole Fraction
    d2yBdz2(2)   = (YB(3)-YB(2))/dz/dz                  ;
    d2yBdz2(3:N) = (YB(4:N+1)+YB(2:N-1)-2*YB(3:N))/dz/dz;
    d2yBdz2(N+1) = (YB(N)-YB(N+1))/dz/dz                ;
    
    %%% Temperature
    d2Tdz2(3:N) = (T(4:N+1)+T(2:N-1)-2*T(3:N))/dz/dz ;
    d2Tdz2(2)   =  4*(Th(2)+T(1)-2*T(2))/dz/dz       ;
    d2Tdz2(N+1) =  4*(Th(N)+T(N+2)-2*T(N+1))/dz/dz   ;

 % 5) Velocities at the boundary interface
    rho_g(2:N) = Ph(2:N)./(R.*Th(2:N));
    alpha(2:N) = (1.75*(1-eb).*rho_g(2:N))./(eb*dp) ;    
    disc       = beta^2 - (4*alpha(2:N).*dpdzh(2:N)) ;
    u(2:N)     = (-beta + sqrt(disc))./(2.*alpha(2:N)) ; 

%% Internal Nodes Temporal Derivatives

% 1) Sorbed Mass Balance (molar loading for CO2 and H2O)
%    1.1 Calculate equilibrium molar loading        
        qstar = IsothermMS(YA, YB, P, T, CondParam, isotherm_parameters);
        qstar1 = qstar(:, 1);
        qstar2 = qstar(:, 2);
        
%    1.2 Calculate the uptake rate  
        dqdt1(2:N+1) = kCO2_des.*(qstar1(2:N+1)-q1(2:N+1));
        dqdt2(2:N+1) = kH2O_des.*(qstar2(2:N+1)-q2(2:N+1));

% % 2) Column Energy Balance
% %       For the energy balance in the column, it is assumed that:
% %       * There is thermal equilibrium between the solid and gas phase
% %       * Conduction occurs through both the solid and gas phase 
% %       * Advection only occurs in the gas phase 
% 
% %     2.1 Calculate the Wall Temperature 
%         dTW2dz2(2:N+1) = (TW(3:N+2) - 2.*TW(2:N+1) + TW(1:N))/dz/dz;  % 2nd order derivative central difference
%         dTWdt(2:N+1)   = (Kw.*Tw1.*dTW2dz2(2:N+1)) + (Tw2.*Tw1.*(T(2:N+1)-TW(2:N+1))) - (Tw3.*Tw1.*(TW(2:N+1)-Tatm));
% 
% %     2.2 Calculate the Total Temperature   
%         sink = (rho_s*C_ps) + (C_pa_CO2.*q1(2:N+1)) + (C_pa_H2O.*q2(2:N+1)); 
%         Pv   = Ph(1:N+1).*u(1:N+1);
%         PvT  = Ph(1:N+1).*u(1:N+1)./Th(1:N+1);
% 
%         dHCO2    = dUb_CO2 - (R.*(T(2:N+1)-T(1:N))); 
%         dHH2O    = dUb_H2O - (R.*(T(2:N+1)-T(1:N)));
% 
%         dTdt1(2:N+1) = (Kz/(1-eb)./sink).*d2Tdz2(2:N+1)                     ; % Temperature change due to conduction
%         dTdt2(2:N+1) = (C_pg/R/dz/e_ratio./sink).*(Pv(2:N+1) - Pv(1:N))     ; % Temperature change due to advection
%         dTdt3(2:N+1) = (1./sink).*T(2:N+1).*(C_pa_CO2.*dqdt1(2:N+1) + C_pa_H2O.*dqdt2(2:N+1))  ; % Temperature change due to reaction
%         dTdt4(2:N+1) = (C_pg/R/e_ratio./sink).*dpdt(2:N+1)                  ; % Temperature change due to pressure
%         dTdt5(2:N+1) = ((-dHCO2./sink).*dqdt1(2:N+1)) + ((-dHH20./sink).*dqdt2(2:N+1)) ;% Temperature change due to adsorption/desorption enthalpy
%         dTdt6(2:N+1) = (Tw4./(e_ratio.*sink)).*(T(2:N+1)-TW(2:N+1))         ; % Temperature change due to the wall
% 
%         dTdt(2:N+1) = dTdt1(2:N+1) - dTdt2(2:N+1) - dTdt3(2:N+1) - dTdt4(2:N+1) + dTdt5(2:N+1) - dTdt6(2:N+1);
% 
% % 3) Total mass balance        
%         dpdt1(2:N+1)  = -T(2:N+1).*(PvT(2:N+1) - PvT(1:N))./dz            ;  % change in pressure due to advection
%         dpdt2(2:N+1)  = -R.*T(2:N+1).*e_ratio.*(dqdt1(2:N+1)+dqdt2(2:N+1));  % change in pressure due to adsorption/desorption
%         dpdt3(2:N+1)  = (P(2:N+1).*dTdt(2:N+1))./T(2:N+1)                ;  % change in pressure due to temperature
% 
%         dpdt(2:N+1)   = dpdt1(2:N+1) + dpdt2(2:N+1) + dpdt3(2:N+1);  
%% Trying this expression for mass balance
% =========================================================================
% 2 & 3) Decoupled Column Energy & Total Mass Balances (Fixed Algebraic Loop)
% =========================================================================

% % 2.1 Calculate Wall Temperature
dTW2dz2(2:N+1) = (TW(3:N+2) - 2.*TW(2:N+1) + TW(1:N))/dz/dz; 
dTWdt(2:N+1) = (Kw.*Tw1.*dTW2dz2(2:N+1)) + (Tw2.*Tw1.*(T(2:N+1)-TW(2:N+1))) - (Tw3.*Tw1.*(TW(2:N+1)-Tatm)); 

% 2.2 Pre-calculate independent parts of Temperature and Pressure
sink = (rho_s*C_ps) + (C_pa_CO2.*q1(2:N+1)) + (C_pa_H2O.*q2(2:N+1)); 
Pv = Ph(1:N+1).*u(1:N+1); 
PvT = Ph(1:N+1).*u(1:N+1)./Th(1:N+1); 

dHCO2    = dUb_CO2 - (R.*(T(2:N+1)-T(1:N))); 
dHH2O    = dUb_H2O - (R.*(T(2:N+1)-T(1:N)));

% Independent components of dT/dt
dTdt1(2:N+1) = (Kz/(1-eb)./sink).*d2Tdz2(2:N+1); 
dTdt2(2:N+1) = (C_pg/R/dz/e_ratio./sink).*(Pv(2:N+1) - Pv(1:N)); 
dTdt3(2:N+1) = (1./sink).*T(2:N+1).*(C_pa_CO2.*dqdt1(2:N+1) + C_pa_H2O.*dqdt2(2:N+1)); 
dTdt5(2:N+1) = ((-dHCO2./sink).*dqdt1(2:N+1)) + ((-dHH2O./sink).*dqdt2(2:N+1)); 
dTdt6(2:N+1) = (Tw4./(e_ratio.*sink)).*(T(2:N+1)-TW(2:N+1)); 

% Group independent temperature terms: dTdt_independent = dTdt1 - dTdt2 - dTdt3 + dTdt5 - dTdt6
dT_ind = dTdt1(2:N+1) - dTdt2(2:N+1) - dTdt3(2:N+1) + dTdt5(2:N+1) - dTdt6(2:N+1);

% Independent components of dP/dt
dpdt1(2:N+1) = -T(2:N+1).*(PvT(2:N+1) - PvT(1:N))./dz; 
dpdt2(2:N+1) = -R.*T(2:N+1).*e_ratio.*(dqdt1(2:N+1)+dqdt2(2:N+1)); 

% Group independent pressure terms
dp_ind = dpdt1(2:N+1) + dpdt2(2:N+1);

% 2.3 Compute Multipliers for the Implicit Link
B_coeff = (C_pg / R / e_ratio ./ sink); % Factor of dp/dt inside dT/dt
D_coeff = P(2:N+1) ./ T(2:N+1);        % Factor of dT/dt inside dp/dt

% 2.4 Calculate Final Decoupled Derivatives Explicitly
dTdt(2:N+1) = (dT_ind - B_coeff .* dp_ind) ./ (1 + B_coeff .* D_coeff);
dpdt(2:N+1) = dp_ind + D_coeff .* dTdt(2:N+1);

% 2.5 Back-calculate individual terms for consistency diagnostics if needed
dTdt4(2:N+1) = B_coeff .* dpdt(2:N+1);
dpdt3(2:N+1) = D_coeff .* dTdt(2:N+1);


%%    
% 4)  Calculate CO2 Mole fraction
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
        
        
 %% Boundary Time Derivatives 
    dyAdt(1)   = dyAdt(2);
    dyAdt(end) = dyAdt(end-1);
    dyBdt(1)   = dyBdt(2);
    dyBdt(end) = dyBdt(end-1);            
    dpdt(1)    = dpdt(2);
    dpdt(end)  = lambda*(PL-P_last)*exp(-lambda*t_des);
    dqdt1(1)   = 0;
    dqdt1(end) = 0;
    dqdt2(1)   = 0;
    dqdt2(end) = 0;
    dTdt(1)    = dTdt(2);
    dTdt(end)  = dTdt(end-1);
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

