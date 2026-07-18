


function PH2O_sat = PressureVaporSat(T)
%----------Saturation Pressure of DI Water------------

% This function is used to find the saturation pressure of water at any given
% temperature between the range 255.9K-373K. The equation used is the 
% Antoine equation with the coefficients adopted from Stull 1947 
% (https://doi.org/10.1021/ie50448a022)

% Antoine Eq: log10(P) = A-(B/(T+C), where: 
% P : Pressure in (bar)
% T : Temperature in (K)
% A : 4.6543, B: 1435.264, C: -64.848


    % log10P    = 4.6543 - (1435.264/(T + (-64.848)))  ; % Antoine Equation, where P is in bar. 
    % PH2O_sat  = (10.^(log10P)).*1e5                  ; % Convert log10 pressure to a number in Pascals from bar. 
    
    % A = 18.678-((T-273.15)/234.5)     ;
    % B = (T-273.15)./(T-16.01)        ;
    % PH2O_sat  = 611.21.*exp(A.*B)     ;

%----------Saturation Pressure of Sea Water (Brackish Water)------------
    PH2O_sat  = (3.8962.*(T.^2)) - (2137.4.*T) + 294066    ;

end
