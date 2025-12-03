

function MassBalance = MassBalanceFxn(a,b,c,d,e,t1,t2,t3,t4,t5, N, L, ColumnParam, CondParam)
    % Extract Necessary Parameters
    r_in = ColumnParam(3);
    eb   = ColumnParam(5);
    mu   = CondParam(7);
    Dm   = CondParam(8);
    R    = CondParam(14);
    dp   = CondParam(17);
    vair = CondParam(23);

    % Calculate relevant values
    Area = pi * r_in^2;
    beta = (150*(1-eb)^2*mu)/(eb^2*dp^2) ; 
    dz   = L / N;
    dz2  = dz/2 ;
    z    = [0 dz2:dz:L-dz2 L] ;
   

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
    [yCO2_A, yH2O_A, Press_A, Temp_A, qCO2_A, qH2O_A, yN2_A] = extract_params(a);
    [yCO2_B, yH2O_B, Press_B, Temp_B, qCO2_B, qH2O_B, yN2_B] = extract_params(b);
    [yCO2_C, yH2O_C, Press_C, Temp_C, qCO2_C, qH2O_C, yN2_C] = extract_params(c);
    [yCO2_D, yH2O_D, Press_D, Temp_D, qCO2_D, qH2O_D, yN2_D] = extract_params(d);
    [yCO2_E, yH2O_E, Press_E, Temp_E, qCO2_E, qH2O_E, yN2_E] = extract_params(e);

    % CO2 Mass Balance Function
    function [Uin, Uout, JMB_CO2, molCO2in, molCO2out, netCO2, netCO2_solid, netCO2_gas, molCO2_sorbed, nCO2in, nCO2out] = CO2_balance(Press, Temp, yCO2, yH2O, qCO2, ~, time)
        
        MW_in    = 0.02896 + (-0.0109)*yH2O(:,1);
        rho_g1   = Press(:,1)./(R*Temp(:,1));  
        alpha_in = (1.75*(1-eb).*rho_g1.*MW_in)./(eb*dp) ;
        disc1    = beta^2 - (4*alpha_in.*(Press(:,2) - Press(:,1)).*(2/dz)) ;
        Uin      = max(((-beta + sqrt(disc1)) ./ (2*alpha_in)), 0) ; 

        MW_out   = 0.02896 + (-0.0109)*yH2O(:,end);
        rho_g2   = Press(:,end)./(R*Temp(:,end));  
        alpha_out= (1.75*(1-eb).*rho_g2.*MW_out)./(eb*dp) ;
        disc2    = beta^2 - (4*alpha_out.*(Press(:,end) - Press(:,end-1)).*(2/dz)) ;
        Uout     = max(((-beta + sqrt(disc2)) ./ (2*alpha_out)), 0) ;         

        % Compute Moles IN & OUT
        nCO2in   = (eb * Area / R) .* (Uin .* Press(:,1) .* yCO2(:,1) ./Temp(:,1));
        molCO2in = trapz(time, nCO2in);
        nCO2out  = (eb * Area / R) .* (Uout .* Press(:,end) .* yCO2(:,end) ./Temp(:,end));
        molCO2out= trapz(time, nCO2out);

        % CO2 moles in Solid phase and Gas Phase
        netCO2_solid = trapz(z, (1 - eb) * Area .* (qCO2(end, :) - qCO2(1, :)));       
        netCO2_gas   = trapz(z, (eb*Area*R^-1).*(((Press(end,:).*yCO2(end,:))./Temp(end,:)) - ((Press(1,:).*yCO2(1,:))./Temp(1,:)))) ; 

        % Net sorbed CO2 and H2O
        netCO2 = molCO2in - molCO2out;
        molCO2_sorbed = netCO2_gas + netCO2_solid;

        JMB_CO2 = abs(((netCO2) - molCO2_sorbed) / molCO2_sorbed) * 100;
    end

    % Compute Mass Balance for Each Step
    [Uin_pz , Uout_pz , JMB_pz_CO2 , molCO2in_pz , molCO2out_pz , netCO2_pz , netCO2_SP_pz , netCO2_GP_pz , molCO2_SB_pz , nCO2in_pz , ~  ]     = CO2_balance(Press_A, Temp_A, yCO2_A, yH2O_A, qCO2_A, qH2O_A, t1);
    [Uin_ads, Uout_ads, JMB_ads_CO2, molCO2in_ads, molCO2out_ads, netCO2_ads, netCO2_SP_ads, netCO2_GP_ads, molCO2_SB_ads, nCO2in_ads, nCO2out_ads  ]     = CO2_balance(Press_B, Temp_B, yCO2_B, yH2O_B, qCO2_B, qH2O_B, t2);
    [Uin_ev , Uout_ev , JMB_ev_CO2 , molCO2in_ev , molCO2out_ev , netCO2_ev , netCO2_SP_ev , netCO2_GP_ev , molCO2_SB_ev , ~ , nCO2out_ev ]     = CO2_balance(Press_C, Temp_C, yCO2_C, yH2O_C, qCO2_C, qH2O_C, t3);
    [Uin_vp , Uout_vp , JMB_vp_CO2 , molCO2in_vp , molCO2out_vp , netCO2_vp , netCO2_SP_vp , netCO2_GP_vp , molCO2_SB_vp , ~ , nCO2out_vp ]     = CO2_balance(Press_D, Temp_D, yCO2_D, yH2O_D, qCO2_D, qH2O_D, t4);
    [Uin_des, Uout_des, JMB_des_CO2, molCO2in_des, molCO2out_des, netCO2_des, netCO2_SP_des, netCO2_GP_des, molCO2_SB_des, ~, nCO2out_des ]     = CO2_balance(Press_E, Temp_E, yCO2_E, yH2O_E, qCO2_E, qH2O_E, t5);

    save("CO2_molarflow.mat", 'nCO2in_pz', 'nCO2in_ads', 'nCO2out_ads', 'nCO2out_ev', 'nCO2out_vp', 'nCO2out_des');
    save("CO2_moles.mat", 'molCO2in_pz', 'molCO2in_ads', 'molCO2out_ads', 'molCO2out_ev', 'molCO2out_vp', 'molCO2out_des');
    
    Step = ["Pressurization"; "CO2 Adsorption"; "Air Evacuation"; "Vapor Stripping"; "Final Desorption"];
    molCO2in    = round([molCO2in_pz; molCO2in_ads; molCO2in_ev; molCO2in_vp; molCO2in_des], 3, 'significant') ; 
    molCO2out   = round([molCO2out_pz; molCO2out_ads; molCO2out_ev; molCO2out_vp; molCO2out_des],3, 'significant') ; 
    netCO2      = round([netCO2_pz; netCO2_ads; netCO2_ev; netCO2_vp; netCO2_des], 3, 'significant');
    netCO2_GP   = round([netCO2_GP_pz; netCO2_GP_ads; netCO2_GP_ev; netCO2_GP_vp; netCO2_GP_des],3, 'significant');
    netCO2_SP   = round([netCO2_SP_pz; netCO2_SP_ads; netCO2_SP_ev; netCO2_SP_vp; netCO2_SP_des], 3, 'significant');
    molCO2_acc  = round([molCO2_SB_pz; molCO2_SB_ads; molCO2_SB_ev; molCO2_SB_vp; molCO2_SB_des], 3, 'significant');
    JMB_CO2     = round([JMB_pz_CO2; JMB_ads_CO2; JMB_ev_CO2; JMB_vp_CO2; JMB_des_CO2], 3, 'significant');

    CO2_MassBalanceTable = table(Step, molCO2in, molCO2out, netCO2,netCO2_GP, netCO2_SP, molCO2_acc, JMB_CO2) ;
    CO2_MassBalanceTable.Properties.VariableNames{'JMB_CO2'} = 'JMB(%)_CO2'; % rename this column to have a % sign

    % Save the table to a file
    writetable(CO2_MassBalanceTable, 'CO2MassBalanceResults.csv');


% H2O Mass Balance Function
    function [JMB_H2O, molH2Oin, molH2Oout, netH2O, netH2O_solid, netH2O_gas, molH2O_sorbed, nH2Oin, nH2Oout] = H2O_balance(Press, Temp, ~, yH2O, ~, qH2O, time)
        
        MW_in    = 0.02896 + (-0.0109)*yH2O(:,1);
        rho_g1   = Press(:,1)./(R*Temp(:,1));  
        alpha_in = (1.75*(1-eb).*rho_g1.*MW_in)./(eb*dp) ;
        disc1    = beta^2 - (4*alpha_in.*(Press(:,2) - Press(:,1)).*(2/dz)) ;
        Uin      = max(((-beta + sqrt(disc1)) ./ (2*alpha_in)), 0) ; 

        MW_out   = 0.02896 + (-0.0109)*yH2O(:,end);
        rho_g2   = Press(:,end)./(R*Temp(:,end));  
        alpha_out= (1.75*(1-eb).*rho_g2.*MW_out)./(eb*dp) ;
        disc2    = beta^2 - (4*alpha_out.*(Press(:,end) - Press(:,end-1)).*(2/dz)) ;
        Uout     = max(((-beta + sqrt(disc2)) ./ (2*alpha_out)), 0) ;         


        % Compute Moles IN & OUT
        nH2Oin    = (eb * Area / R) .* (Uin .* Press(:,1) .* yH2O(:,1) ./Temp(:,1));
        molH2Oin  = trapz(time, nH2Oin);
        nH2Oout   = (eb * Area / R) .* (Uout .* Press(:,end) .* yH2O(:,end) ./Temp(:,end));
        molH2Oout = trapz(time, nH2Oout);

        % H2O moles in Solid phase and Gas Phase
        netH2O_solid = trapz(z, (1 - eb) * Area .* (qH2O(end, :) - qH2O(1, :)));
        netH2O_gas   = trapz(z, (eb*Area*R^-1).*(((Press(end,:).*yH2O(end,:))./Temp(end,:)) - ((Press(1,:).*yH2O(1,:))./Temp(1,:)))) ; 
                
        % Net sorbed CO2 and H2O
        netH2O = molH2Oin - molH2Oout;
        molH2O_sorbed = netH2O_gas + netH2O_solid;

        JMB_H2O = abs(((netH2O) - molH2O_sorbed) / molH2O_sorbed) * 100;
    end

% Compute Mass Balance for Each Step
    [JMB_pz_H2O , molH2Oin_pz , molH2Oout_pz , netH2O_pz , netH2O_SP_pz , netH2O_GP_pz , molH2O_SB_pz , nH2Oin_pz , ~        ]     = H2O_balance(Press_A, Temp_A, yCO2_A, yH2O_A, qCO2_A, qH2O_A, t1);
    [JMB_ads_H2O, molH2Oin_ads, molH2Oout_ads, netH2O_ads, netH2O_SP_ads, netH2O_GP_ads, molH2O_SB_ads, nH2Oin_ads , nH2Oout_ads]     = H2O_balance(Press_B, Temp_B, yCO2_B, yH2O_B, qCO2_B, qH2O_B, t2);
    [JMB_ev_H2O , molH2Oin_ev , molH2Oout_ev , netH2O_ev , netH2O_SP_ev , netH2O_GP_ev , molH2O_SB_ev , ~, nH2Oout_ev        ]     = H2O_balance(Press_C, Temp_C, yCO2_C, yH2O_C, qCO2_C, qH2O_C, t3);
    [JMB_vp_H2O , molH2Oin_vp , molH2Oout_vp , netH2O_vp , netH2O_SP_vp , netH2O_GP_vp , molH2O_SB_vp , nH2Oin_vp, nH2Oout_vp]     = H2O_balance(Press_D, Temp_D, yCO2_D, yH2O_D, qCO2_D, qH2O_D, t4);
    [JMB_des_H2O, molH2Oin_des, molH2Oout_des, netH2O_des, netH2O_SP_des, netH2O_GP_des, molH2O_SB_des, ~, nH2Oout_des       ]     = H2O_balance(Press_E, Temp_E, yCO2_E, yH2O_E, qCO2_E, qH2O_E, t5);

    % H2O Outputs
    save("H2O_molarflow.mat", 'nH2Oin_pz', 'nH2Oin_ads', 'nH2Oout_ads', 'nH2Oout_ev', 'nH2Oout_vp', 'nH2Oout_des');
    save("H2O_moles.mat",'molH2Oin_pz', 'molH2Oin_ads', 'molH2Oin_vp', 'molH2Oout_ads','molH2Oout_ev', 'molH2Oout_vp', 'molH2Oout_des');  

    Step = ["Pressurization"; "CO2 Adsorption"; "Air Evacuation"; "Vapor Stripping"; "Final Desorption"];
    molH2Oin    = round([molH2Oin_pz; molH2Oin_ads; molH2Oin_ev; molH2Oin_vp; molH2Oin_des], 3, 'significant') ; 
    molH2Oout   = round([molH2Oout_pz; molH2Oout_ads; molH2Oout_ev; molH2Oout_vp; molH2Oout_des],3, 'significant') ; 
    netH2O      = round([netH2O_pz; netH2O_ads; netH2O_ev; netH2O_vp; netH2O_des], 3, 'significant');
    netH2O_GP   = round([netH2O_GP_pz; netH2O_GP_ads; netH2O_GP_ev; netH2O_GP_vp; netH2O_GP_des],3, 'significant');
    netH2O_SP   = round([netH2O_SP_pz; netH2O_SP_ads; netH2O_SP_ev; netH2O_SP_vp; netH2O_SP_des], 3, 'significant');
    molH2O_acc  = round([molH2O_SB_pz; molH2O_SB_ads; molH2O_SB_ev; molH2O_SB_vp; molH2O_SB_des], 3, 'significant');
    JMB_H2O     = round([JMB_pz_H2O; JMB_ads_H2O; JMB_ev_H2O; JMB_vp_H2O; JMB_des_H2O], 3, 'significant');

    H2O_MassBalanceTable = table(Step, molH2Oin, molH2Oout, netH2O,netH2O_GP, netH2O_SP, molH2O_acc, JMB_H2O) ;
    H2O_MassBalanceTable.Properties.VariableNames{'JMB_H2O'} = 'JMB(%)_H2O'; % rename this column to have a % sign

    % Save the table to a file
    writetable(H2O_MassBalanceTable, 'H2OMassBalanceResults.csv');


% N2 Mass Balance Function
    function [molN2in, molN2out, nN2in, nN2out] = N2_balance(Press, Temp, ~, yH2O, ~, ~, yN2, time)
        
        MW_in    = 0.02896 + (-0.0109)*yH2O(:,1);
        rho_g1   = Press(:,1)./(R*Temp(:,1));  
        alpha_in = (1.75*(1-eb).*rho_g1.*MW_in)./(eb*dp) ;
        disc1    = beta^2 - (4*alpha_in.*(Press(:,2) - Press(:,1)).*(2/dz)) ;
        Uin      = max(((-beta + sqrt(disc1)) ./ (2*alpha_in)), 0) ; 

        MW_out    = 0.02896 + (-0.0109)*yH2O(:,end);
        rho_g2    = Press(:,end)./(R*Temp(:,end));  
        alpha_out = (1.75*(1-eb).*rho_g2.*MW_out)./(eb*dp) ;
        disc2     = beta^2 - (4*alpha_out.*(Press(:,end) - Press(:,end-1)).*(2/dz)) ;
        Uout      = max(((-beta + sqrt(disc2)) ./ (2*alpha_out)), 0) ;         

        % Compute Moles IN & OUT
        nN2in   = (eb * Area / R) .* (Uin .* Press(:,1) .* yN2(:,1) ./Temp(:,1));
        molN2in = trapz(time, nN2in);
        nN2out  = (eb * Area / R) .* (Uout .* Press(:,end) .* yN2(:,end) ./Temp(:,end));
        molN2out = trapz(time, nN2out);

    end

% Compute Mass Balance for Each Step
    [molN2in_pz , ~ , nN2in_pz , ~ ]     = N2_balance(Press_A, Temp_A, yCO2_A, yH2O_A, qCO2_A, qH2O_A, yN2_A, t1);
    [molN2in_ads, ~, nN2in_ads , ~ ]     = N2_balance(Press_B, Temp_B, yCO2_B, yH2O_B, qCO2_B, qH2O_B, yN2_B, t2);
    [~ , molN2out_ev , ~, nN2out_ev]     = N2_balance(Press_C, Temp_C, yCO2_C, yH2O_C, qCO2_C, qH2O_C, yN2_C, t3);
    [~ , molN2out_vp , ~, nN2out_vp]     = N2_balance(Press_D, Temp_D, yCO2_D, yH2O_D, qCO2_D, qH2O_D, yN2_D, t4);
    [~, molN2out_des, ~, nN2out_des]     = N2_balance(Press_E, Temp_E, yCO2_E, yH2O_E, qCO2_E, qH2O_E, yN2_E, t5);

    % N2 Outputs
    save("N2_molarflow.mat" , 'nN2in_pz', 'nN2in_ads', 'nN2out_ev', 'nN2out_vp', 'nN2out_des');
    save("N2_moles.mat",'molN2in_pz', 'molN2in_ads', 'molN2out_ev', 'molN2out_vp', 'molN2out_des');    
    
 % Final Output
    CO2_cycle_balance = (molCO2out_ads + molCO2out_ev + molCO2out_vp + molCO2out_des) / (molCO2in_pz + molCO2in_ads);
    H2O_cycle_balance = (molH2Oout_ads + molH2Oout_ev + molH2Oout_vp + molH2Oout_des) / (molH2Oin_pz + molH2Oin_ads);
    % Velocities          = [Uin_pz, Uin_ads, Uin_vp, Uout_ads, Uout_ev, Uout_vp, Uout_des] ; 

    MassBalance.CO2_cycle_balance      = CO2_cycle_balance     ;
    MassBalance.CO2_MassBalanceTable   = CO2_MassBalanceTable  ;
    MassBalance.H2O_cycle_balance      = H2O_cycle_balance     ;
    MassBalance.H2O_MassBalanceTable   = H2O_MassBalanceTable  ;
end
