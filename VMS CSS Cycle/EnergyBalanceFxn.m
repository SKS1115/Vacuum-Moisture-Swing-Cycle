function EnergyBalance = EnergyBalanceFxn(a,b,c,d,e,t1,t2,t3,t4,t5,N,L,ColumnParam,CondParam,isotherm_parameters) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Energy Balance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%% FUNCTIONS USED
% j*trapz(i)    : trapz is a numerical integration function that integrates
                % the data set i based on j spacing. 

% Input parameters
    % InputParameters = ProcessInputParameters(N);        % Call input parameter function into this script of code
    % ColumnParam     = InputParameters{1};               % Retrieve respective set of saved parameters from input parameter function
    % CondParam       = InputParameters{2};               % Retrieve respective set of saved parameters from input parameter function
    % isotherm_parameters = InputParameters{3};           % Retrieve respective set of saved parameters from input parameter function

% Retrieve Necessary parameters and calculate constants
    r_in    = ColumnParam(3) ;
    eb      = ColumnParam(5) ;          
    rho_s   = CondParam(1)   ;
    C_pg    = CondParam(3)   ;
    C_pa_CO2= CondParam(4)   ;
    C_ps    = CondParam(5)   ;
    mu      = CondParam(7)   ;
    Dm      = CondParam(8)   ;
    R       = CondParam(14)  ;           
    dp      = CondParam(17)  ;
    C_pa_H2O= CondParam(19);
    dUb_CO2 = isotherm_parameters(4) ;
    dUb_H2O =  isotherm_parameters(10) ;
        

    % Calculate relevant values    
    Area  = (pi.*r_in^2)          ;               
    dz    = L/N                   ;
    dz2   = dz/2                  ;    
    beta  = (150*(1-eb)^2*mu)/(eb^2*dp^2) ; 
    z     = [0 dz2:dz:L-dz2 L]    ;      % axis from inlet-outlet boundary

 % Extract Data Function
    function [yCO2, yH2O, Press, Temp, qCO2, qH2O, yN2] = extract_params(matrix)
        yCO2 = matrix(:, 1:N+2);
        yH2O = matrix(:, N+3:2*N+4);
        Press= matrix(:, 2*N+5:3*N+6);
        Temp = matrix(:, 5*N+11:6*N+12);
        qCO2 = matrix(:, 3*N+7:4*N+8);
        qH2O = matrix(:, 4*N+9:5*N+10);
        yN2  = max(1 - (yCO2 + yH2O), 0); % Ensure no negative values
    end

    % Extract Data
    [yCO2_A, yH2O_A, Press_A, Temp_A, qCO2_A, qH2O_A, ~] = extract_params(a);
    [yCO2_B, yH2O_B, Press_B, Temp_B, qCO2_B, qH2O_B, ~] = extract_params(b);
    [yCO2_C, yH2O_C, Press_C, Temp_C, qCO2_C, qH2O_C, ~] = extract_params(c);
    [yCO2_D, yH2O_D, Press_D, Temp_D, qCO2_D, qH2O_D, ~] = extract_params(d);
    [yCO2_E, yH2O_E, Press_E, Temp_E, qCO2_E, qH2O_E, ~] = extract_params(e); 

    function [Heat_in, Heat_out, Heat_gen, Heat_SP, Heat_GP, Heat_SB, JEB] =StepEnergyBalance(Press, Temp, ~, yH2O, qCO2, qH2O, time)

        MW_in    = 0.02896 + (-0.0109)*yH2O(:,1);
        rho_g1   = Press(:,1)./(R*Temp(:,1));  
        alpha_in = (1.75*(1-eb).*rho_g1.*MW_in)./(eb*dp) ;
        disc1   = beta^2 - (4*alpha_in.*(Press(:,2) - Press(:,1)).*(2/dz)) ;
        Uin     = max(((-beta + sqrt(disc1)) ./ (2*alpha_in)), 0) ; 

        MW_out    = 0.02896 + (-0.0109)*yH2O(:,end);
        rho_g2    = Press(:,end)./(R*Temp(:,end));  
        alpha_out = (1.75*(1-eb).*rho_g2.*MW_out)./(eb*dp) ;
        disc2   = beta^2 - (4*alpha_out.*(Press(:,end) - Press(:,end-1)).*(2/dz)) ;
        Uout    = max(((-beta + sqrt(disc2)) ./ (2*alpha_out)), 0) ;    

        %%% Heat entering the column 
        VTrho_in     = (Uin.*Press(:,1).*R^-1)                        ; % [v*T*rho], where rho = P/RT
        Heat_in      = eb*Area*C_pg*(trapz(time, VTrho_in)) ; % [Joules]

        %%% Heat leaving the column 
        VTrho_out     = (Uout.*Press(:,end).*R^-1)                        ; % [v*T*rho], where rho = P/RT
        Heat_out      = eb*Area*C_pg*(trapz(time, VTrho_out)) ; % [Joules]
        
        %%% Calculate the Heat generated        
        dH_qCO2     = -dUb_CO2.*(qCO2(end,:) - qCO2(1,:))     ; % dH*qCO2 @t=final - dH*qCO2 @t=initial 
        dH_qH2O     = -dUb_H2O.*(qH2O(end,:) - qH2O(1,:))    ; % dH*qH2O @t=final - dH*qH2O @t=initial     
        Heat_gen    = (1-eb)*Area*(trapz(z, dH_qCO2) + trapz(z,dH_qH2O));
    
        %%% Heat in the solid phase 
        Temp_diff  = Temp(end,:) - Temp(1,:)           ; % temp @t=final - temp @t=initial       
        Heat_SP    = (1-eb)*Area*C_ps*rho_s*(trapz(z, Temp_diff))   ;

        %%% Heat in the gas phase
        T_rhog_diff = (Press(end,:) - Press(1,:))./R       ; % temp*rhog @t=final - initial, where rhog =P/RT
        Heat_GP     = eb*Area*C_pg*(trapz(z, T_rhog_diff))    ; 
    
        %%% Heat required for sorption 
        Tq_CO2_diff   = C_pa_CO2.*((Temp(end,:).*qCO2(end,:)) - (Temp(1,:).*qCO2(1,:)))       ;  % temp*qCO2 @t=final - initial 
        Tq_H2O_diff   = C_pa_H2O.*((Temp(end,:).*qH2O(end,:)) - (Temp(1,:).*qH2O(1,:)))       ;  % temp*qH2O @t=final - initial 
        Heat_SB       = (1-eb)*Area*(trapz(z, Tq_CO2_diff + Tq_H2O_diff));
    
        %%% Calculate the Energy error balance    
        Heat_gain   = Heat_in + Heat_gen                                       ; % total heat energy in and generated
        Heat_lost   = Heat_out + Heat_SP + Heat_GP + Heat_SB ; % total heat energy used and lost
        JEB         = abs((Heat_gain - Heat_lost)/Heat_lost)*100      ;    % energy balance error
    end

    % Compute Energy Balance for Each Step
    [Heat_in_pz, Heat_out_pz, Heat_gen_pz, Heat_SP_pz, Heat_GP_pz, Heat_SB_pz, JEB_pz]    =StepEnergyBalance(Press_A, Temp_A, yCO2_A, yH2O_A, qCO2_A, qH2O_A, t1);
    [Heat_in_ad, Heat_out_ad, Heat_gen_ad, Heat_SP_ad, Heat_GP_ad, Heat_SB_ad, JEB_ad]    =StepEnergyBalance(Press_B, Temp_B, yCO2_B, yH2O_B, qCO2_B, qH2O_B, t2);
    [Heat_in_ev, Heat_out_ev, Heat_gen_ev, Heat_SP_ev, Heat_GP_ev, Heat_SB_ev, JEB_ev]    =StepEnergyBalance(Press_C, Temp_C, yCO2_C, yH2O_C, qCO2_C, qH2O_C, t3);
    [Heat_in_vs, Heat_out_vs, Heat_gen_vs, Heat_SP_vs, Heat_GP_vs, Heat_SB_vs, JEB_vs]    =StepEnergyBalance(Press_D, Temp_D, yCO2_D, yH2O_D, qCO2_D, qH2O_D, t4);
    [Heat_in_ds, Heat_out_ds, Heat_gen_ds, Heat_SP_ds, Heat_GP_ds, Heat_SB_ds, JEB_ds]    =StepEnergyBalance(Press_E, Temp_E, yCO2_E, yH2O_E, qCO2_E, qH2O_E, t5);

    Step = ["Pressurization"; "CO2 Sorption"; "Air Evacuation"; "Vapor Stripping"; "Final Desorption"];
    Heat_in   = round([Heat_in_pz ; Heat_in_ad ; Heat_in_ev ; Heat_in_vs ; Heat_in_ds ], 3, 'significant') ; 
    Heat_out  = round([Heat_out_pz; Heat_out_ad; Heat_out_ev; Heat_out_vs; Heat_out_ds], 3, 'significant') ; 
    Heat_gen  = round([Heat_gen_pz; Heat_gen_ad; Heat_gen_ev; Heat_gen_vs; Heat_gen_ds], 3, 'significant');
    Heat_SP   = round([Heat_SP_pz ; Heat_SP_ad ; Heat_SP_ev ; Heat_SP_vs ; Heat_SP_ds ], 3, 'significant');
    Heat_GP   = round([Heat_GP_pz ; Heat_GP_ad ; Heat_GP_ev ; Heat_GP_vs ; Heat_GP_ds ], 3, 'significant');
    Heat_SB   = round([Heat_SB_pz ; Heat_SB_ad ; Heat_SB_ev ; Heat_SB_vs ; Heat_SB_ds ], 3, 'significant');
    JEB       = round([JEB_pz; JEB_ad; JEB_ev; JEB_vs; JEB_ds], 3, 'significant');

    EnergyBalanceTable = table(Step, Heat_in, Heat_out, Heat_gen, Heat_SP, Heat_GP, Heat_SB, JEB) ;
    EnergyBalanceTable.Properties.VariableNames{'JEB'} = 'JEB(%)'; % rename this column to have a % sign
    % Save the table to a file
    writetable(EnergyBalanceTable, 'EnergyBalanceResults.csv');

    % Final OUtput
    EnergyBalance = EnergyBalanceTable ;
end