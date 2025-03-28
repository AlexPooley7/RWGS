% CE30243 - Individual Design Project
% Description - Models the transient behaviour of the reverse water-gas
% shift reaction (RWGS)
% Last edited: 28/03/2025
% Last commit: 28/03/2025
% Author: Alexander Pooley
% -------------------------------------------------------------
% Functions:
% 
% -------------------------------------------------------------
% Structures:
% 
% -------------------------------------------------------------
% Engineering variables:
% 
% -------------------------------------------------------------

clc
clear

% Define constants:
params = struct(); % Initialise the params structure

% Misc inlet parameters
params.inlet.H2O = 0*1000/(60*60*24); % kmol/day -> mol/s
params.inlet.CH4 = 610.43*1000/(60*60*24); % kmol/day -> mol/s
params.inlet.gases = 1295.90*1000/(60*60*24); % kmol/day -> mol/s

% Arrhenius constants
params.arr.preExpFactor = 1.7*10^3;     % m3 kg-1 s-1
params.arr.activationEnergy = 68.5;     % kJ mol-1
params.arr.gasConst = 8.314;            % J/(mol·K)

% Energy balance constants
params.eb.enthalpyReaction = 41.2; % kJ/mol
    
% Carbon Dioxide
params.eb.CO2.Fin = (725.81+6532.28)*1000/(60*60*24); % kmol/day -> mol/s
params.eb.CO2.A = 25.0;
params.eb.CO2.B = 55.2;
params.eb.CO2.C = -33.7;
params.eb.CO2.D = 7.9;
params.eb.CO2.E = -0.1;

% Hydrogen 
params.eb.H2.Fin = 22951.16*1000/(60*60*24); % kmol/day -> mol/s
% params.eb.H2.A = 18.6; % 1000-2500
% params.eb.H2.B = 12.3;
% params.eb.H2.C = -2.9;
% params.eb.H2.D = 0.3;
% params.eb.H2.E = 2.0;

params.eb.H2.A = 33.066; % these are 298-100K
params.eb.H2.B = -11.36;
params.eb.H2.C = 11.433;
params.eb.H2.D = -2.77;
params.eb.H2.E = -0.159;

% Carbon monoxide
params.eb.CO.Fin = 7381.52*1000/(60*60*24); % mol/s
params.eb.CO.A = 25.6;
params.eb.CO.B = 6.1;
params.eb.CO.C = 4.1;
params.eb.CO.D = -2.7;
params.eb.CO.E = 0.1;

% Methane
params.eb.CH4.Fin = params.inlet.CH4; % mol/s
params.eb.CH4.A = -0.7; 
params.eb.CH4.B = 108.5;
params.eb.CH4.C = -42.5;
params.eb.CH4.D = 5.9;
params.eb.CH4.E = 0.7;

% Define initial conditions 
FA0 = params.eb.CO2.Fin; 
FB0 = params.eb.H2.Fin; 
FC0 = params.eb.CO.Fin;
FD0 = 0; 
T0 = 1173;

% Assign initial conditions to vector
Y0 = [FA0, FB0, FC0, FD0, T0];

Wspan = [0 100];

% Pass params to odeSolver using an anonymous function
[w,Y] = ode15s(@(w,Y) odeSolver(w,Y,params), Wspan, Y0); 

T = Y(:,1); FA = Y(:,2); FB = Y(:,3); FC = Y(:,4); FD = Y(:,5);

% Plot all variables
figure;

subplot(2,1,1);
plot(w, FA, 'r', w, FB, 'b', w, FC, 'g', w, FD, 'm', 'LineWidth', 1.5);
xlabel('Weight of Catalyst (kg)');
ylabel('Molar Flow Rate (mol/s)');
legend('CO2 (FA)', 'H2 (FB)', 'CO (FC)', 'H2O (FD)');
title('Molar Flow Rates vs. Catalyst Weight');
grid on;

subplot(2,1,2);
plot(w, T, 'k', 'LineWidth', 1.5);
xlabel('Weight of Catalyst (kg)');
ylabel('Temperature (K)');
title('Temperature vs. Catalyst Weight');
grid on;

checkCapacity(params); % Check the heat capacity formula


% ODE Solver Function
function dYdt = odeSolver(w,Y,params) %#ok<INUSD> 

% Extract state variables
FA = Y(1); FB = Y(2); FC = Y(3); FD = Y(4); T = Y(5);

% Rate constant calculations using the Arrhenius equation
k = params.arr.preExpFactor * exp(-params.arr.activationEnergy / (params.arr.gasConst * T));

% Mole fraction calculations 
totalMol = FA + FB + FC + FD + params.inlet.CH4 + params.inlet.gases;
molFractionCO2 = FA / totalMol;
molFractionH2 = FB / totalMol;

% Assume pressure P = 1 atm (if needed, define it properly)
P = 22; % Placeholder, modify if needed

% Rate of reaction calculations
rRWGS = k * molFractionCO2 * molFractionH2 * ((P^2) / (params.arr.gasConst^2) * (T^2)); 

% Mole balance ODEs
dFA_dw = -rRWGS;
dFB_dw = -rRWGS;
dFC_dw = rRWGS;
dFD_dw = rRWGS;

cpCO2 = schomate(params, T / 1000, 'CO2'); % Convert to kK
cpH2 = schomate(params, T / 1000, 'H2');
cpCO = schomate(params, T / 1000, 'CO');
cpCH4 = schomate(params, T / 1000, 'CH4');

sumNcp = cpCO2*params.eb.CO2.Fin + cpH2*params.eb.H2.Fin + cpCO*params.eb.CO.Fin + cpCH4*params.eb.CH4.Fin;

dT_dw = (params.eb.enthalpyReaction * rRWGS) / (sumNcp);

% Output vector for ODE solver
dYdt = [dFA_dw; dFB_dw; dFC_dw; dFD_dw; dT_dw];

end

%% Schomate Equation Function
function cp = schomate(params,T,component)

switch component
    case 'CO2'
        A = params.eb.CO2.A;
        B = params.eb.CO2.B;
        C = params.eb.CO2.C;
        D = params.eb.CO2.D;
        E = params.eb.CO2.E;
    case 'H2'
        A = params.eb.H2.A;
        B = params.eb.H2.B;
        C = params.eb.H2.C;
        D = params.eb.H2.D;
        E = params.eb.H2.E;
    case 'CO'
        A = params.eb.CO.A;
        B = params.eb.CO.B;
        C = params.eb.CO.C;
        D = params.eb.CO.D;
        E = params.eb.CO.E;
    case 'CH4'
        A = params.eb.CH4.A;
        B = params.eb.CH4.B;
        C = params.eb.CH4.C;
        D = params.eb.CH4.D;
        E = params.eb.CH4.E;
end

% Schomate equation for specific heat capacity (cp)
cp = A + B*T + C*T^2 + D*T^3 + E/T^2 ;
end

%%
function checkCapacity(params)

% Define a range of temperatures for plotting
T_range = linspace(300, 1500, 100);  % Temperature range from 300 K to 1500 K

% Calculate cp for each component over the temperature range
cp_CO2_values = arrayfun(@(T) schomate(params, T / 1000, 'CO2'), T_range);  % Convert T to kK
cp_H2_values = arrayfun(@(T) schomate(params, T / 1000, 'H2'), T_range);
cp_CO_values = arrayfun(@(T) schomate(params, T / 1000, 'CO'), T_range);
cp_CH4_values = arrayfun(@(T) schomate(params, T / 1000, 'CH4'), T_range);

% Plot cp values as a function of temperature
figure;
hold on;
plot(T_range, cp_CO2_values, 'r', 'LineWidth', 1.5);
plot(T_range, cp_H2_values, 'b', 'LineWidth', 1.5);
plot(T_range, cp_CO_values, 'g', 'LineWidth', 1.5);
plot(T_range, cp_CH4_values, 'm', 'LineWidth', 1.5);
hold off;

xlabel('Temperature (K)');
ylabel('Specific Heat Capacity (cp) [J/mol·K]');
legend('CO2', 'H2', 'CO', 'CH4');
title('Specific Heat Capacity vs Temperature');
grid on;

end