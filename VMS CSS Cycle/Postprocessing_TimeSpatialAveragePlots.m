

% Input parameters
    InputParameters = ProcessInputParameters(N,RH,ID,OD,L,vair,vvs,t_cs,t_vs,t_des)  ;      % Call input parameter function into this script of code
    ColumnParam     = InputParameters{1};               % Retrieve respective set of saved parameters from input parameter function
    CondParam       = InputParameters{2};               % Retrieve respective set of saved parameters from input parameter function
    isotherm_parameters = InputParameters{3};           % Retrieve4 respective set of saved parameters from input parameter function

% Retrieve Necessary parameters and calculate constants
    rho_s   = CondParam(1)  ;  
    r_in    = ColumnParam(3) ;
    eb      = ColumnParam(5) ;          
    rp      = ColumnParam(7) ;
    mu      = CondParam(7)   ;
    R       = CondParam(14)  ;          
    Viscous = (4/150)*(eb/(1-eb))^2*(rp^2/mu);
    Area    = (pi.*r_in^2)          ;               
    dz      = L/N                   ;
    dz2     = dz/2                  ;
    z       = [0 dz2:dz:L-dz2 L]    ;      % axis from inlet-outlet boundary
    z2      = 0:dz:L;


% Extract and define all the calculated unknown parameters for each step
    yCO2_A    = a(:,1:N+2)        ; yH2O_A    = a(:,N+3:2*N+4)        ; 
    Press_A   = a(:,2*N+5:3*N+6)  ; Temp_A    = a(:, 5*N+11:6*N+12)   ;
    qCO2_A    = a(:,3*N+7:4*N+8)  ; qH2O_A    = a(:,4*N+9:5*N+10)     ;
    yN2_A     = 1-(yCO2_A+yH2O_A); 
    PCO2_A    = Press_A.*yCO2_A  ; PH2O_A    = Press_A.*yH2O_A  ; 
   
    yCO2_B    = b(:,1:N+2)        ; yH2O_B    = b(:,N+3:2*N+4)        ; 
    Press_B   = b(:,2*N+5:3*N+6)  ; Temp_B    = b(:, 5*N+11:6*N+12)   ;
    qCO2_B    = b(:,3*N+7:4*N+8)  ; qH2O_B    = b(:,4*N+9:5*N+10)     ;
    yN2_B     = 1-(yCO2_B+yH2O_B);
    PCO2_B    = Press_B.*yCO2_B  ; PH2O_B    = Press_B.*yH2O_B  ; 

    yCO2_C    = c(:,1:N+2)        ; yH2O_C    = c(:,N+3:2*N+4)        ; 
    Press_C   = c(:,2*N+5:3*N+6)  ; Temp_C    = c(:, 5*N+11:6*N+12)   ;
    qCO2_C    = c(:,3*N+7:4*N+8)  ; qH2O_C    = c(:,4*N+9:5*N+10)     ;
    yN2_C     = 1-(yCO2_C+yH2O_C) ;
    PCO2_C    = Press_C.*yCO2_C   ; PH2O_C    = Press_C.*yH2O_C  ;
    
    yCO2_D    = d(:,1:N+2)        ; yH2O_D    = d(:,N+3:2*N+4)        ; 
    Press_D   = d(:,2*N+5:3*N+6)  ; Temp_D    = d(:, 5*N+11:6*N+12)   ;
    qCO2_D    = d(:,3*N+7:4*N+8)  ; qH2O_D    = d(:,4*N+9:5*N+10)     ;
    yN2_D     = 1-(yCO2_D+yH2O_D);
    PCO2_D    = Press_D.*yCO2_D   ; PH2O_D    = Press_D.*yH2O_D  ;
  
    yCO2_E    = e(:,1:N+2)        ; yH2O_E    = e(:,N+3:2*N+4)        ; 
    Press_E   = e(:,2*N+5:3*N+6)  ; Temp_E    = e(:, 5*N+11:6*N+12)   ;
    qCO2_E    = e(:,3*N+7:4*N+8)  ; qH2O_E    = e(:,4*N+9:5*N+10)     ;
    yN2_E     = 1-(yCO2_E+yH2O_E);
    PCO2_E    = Press_E.*yCO2_E   ; PH2O_E    = Press_E.*yH2O_E  ;


% Initialize cumulative variables for time and parameters
cumulative_time = [];%zeros(1,length(t1));
cumulative_qCO2 = [];
cumulative_qH2O = [];
cumulative_yCO2 = [];
cumulative_yH2O = [];
cumulative_Temp = [];
cumulative_Press = [];
cumulative_PCO2  = [];
cumulative_PH2O  = [];

% Define all steps and corresponding variables
steps = {
    'Pressurization', t1, qCO2_A, qH2O_A, yCO2_A, yH2O_A, Temp_A, Press_A, PCO2_A, PH2O_A;
    'CO_2 Sorption', t2, qCO2_B, qH2O_B, yCO2_B, yH2O_B, Temp_B, Press_B, PCO2_B, PH2O_B;
    'Air Evacuation', t3, qCO2_C, qH2O_C, yCO2_C, yH2O_C, Temp_C, Press_C, PCO2_C, PH2O_C;    
    'Vapor Stripping', t4, qCO2_D, qH2O_D, yCO2_D, yH2O_D, Temp_D, Press_D, PCO2_D, PH2O_D;
    'Final Desorption', t5, qCO2_E, qH2O_E, yCO2_E, yH2O_E, Temp_E, Press_E, PCO2_E, PH2O_E;
};

% Loop through each step
for step_idx = 1:size(steps, 1)
    % Extract variables for the current step
    step_name = steps{step_idx, 1};
    t = steps{step_idx, 2};
    qCO2_data = steps{step_idx, 3};
    qH2O_data = steps{step_idx, 4};
    yCO2_data = steps{step_idx, 5};
    yH2O_data = steps{step_idx, 6};
    Temp_data = steps{step_idx, 7};
    Press_data = steps{step_idx, 8};
    PCO2_data = steps{step_idx, 9};
    PH2O_data = steps{step_idx, 10};

    % Compute the averaged values
    qCO2 = arrayfun(@(i) trapz(z, qCO2_data(i, :)) / L, 1:length(t));
    qH2O = arrayfun(@(i) trapz(z, qH2O_data(i, :)) / L, 1:length(t));
    yCO2 = arrayfun(@(i) trapz(z, yCO2_data(i, :)) / L, 1:length(t));
    yH2O = arrayfun(@(i) trapz(z, yH2O_data(i, :)) / L, 1:length(t));
    Temp = arrayfun(@(i) trapz(z, Temp_data(i, :)) / L, 1:length(t));
    Press = arrayfun(@(i) trapz(z, Press_data(i, :)) / L, 1:length(t));
    PCO2 = arrayfun(@(i) trapz(z, PCO2_data(i, :)) / L, 1:length(t));
    PH2O = arrayfun(@(i) trapz(z, PH2O_data(i, :)) / L, 1:length(t));

    % Append to cumulative variables
    if isempty(cumulative_time)
        cumulative_time = t';
    else
        cumulative_time = [cumulative_time, cumulative_time(end) + t'];
    end
    cumulative_qCO2 = [cumulative_qCO2, qCO2]; 
    cumulative_qH2O = [cumulative_qH2O, qH2O];
    cumulative_yCO2 = [cumulative_yCO2, yCO2];
    cumulative_yH2O = [cumulative_yH2O, yH2O];
    cumulative_Temp = [cumulative_Temp, Temp];
    cumulative_Press = [cumulative_Press, Press];
    cumulative_PCO2 = [cumulative_PCO2, PCO2];
    cumulative_PH2O = [cumulative_PH2O, PH2O];
end

% Define colors for each step
colors = [
    0.8500, 0.3250, 0.0980;  % Step 1: Red
    0.4660, 0.6740, 0.1880;  % Step 2: Green
    0.0000, 0.4470, 0.7410;  % Step 3: Blue

    0.4940, 0.1840, 0.5560;  % Step 5: Purple
    0.6350, 0.0780, 0.1840   % Step 6: Dark Red
];

% Plot each parameter with unique colors
y_labels = {'qCO_2 [mmol/g]', 'qH_2O [mmol/g]', 'yCO_2 [-]', 'yH_2O [-]', ...
    'Temperature [K]', 'Pressure [kPa]', 'PCO_2 [kPa]', 'PH_2O [kPa]'};
parameters = {cumulative_qCO2/1000, cumulative_qH2O/1000, cumulative_yCO2, ...
    cumulative_yH2O, cumulative_Temp, cumulative_Press*1e-3, cumulative_PCO2*1e-3, cumulative_PH2O*1e-3};

% Loop through each parameter
for param_idx = 1:length(parameters)
    figure(param_idx);
    hold on;
    for step_idx = 1:size(steps, 1)
        % Extract time range for this step
        if step_idx == 1
            start_idx = 1;
        else
            start_idx = sum(cellfun(@length, steps(1:step_idx-1, 2))) + 1;
        end
        end_idx = start_idx + length(steps{step_idx, 2}) - 1;
        
        % Extract data for this step
        time_segment = cumulative_time(start_idx:end_idx);
        param_segment = parameters{param_idx}(start_idx:end_idx);
        
        % Plot with unique color
        plot(time_segment/60, param_segment, 'LineWidth', 3.0, ...
            'Color', colors(step_idx, :));        
        
    end
    hold off;
    
    xlabel('Time [min]', 'FontSize', 20);
    ylabel(y_labels{param_idx}, 'FontSize', 20);
    legend(steps(:, 1), 'Location', 'Best','FontSize', 14);
    legend boxoff    
    set(gca, 'LineWidth', 1.5); 
    
    ax = gca; ax.YAxis.Exponent = 0; ax.XAxis.Exponent = 0;
    ax.FontSize = 20;

end

