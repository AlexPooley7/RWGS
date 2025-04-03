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
params.inlet.temp = 1000; % K
params.inlet.pres = 22*100000; % bar -> Pa

% Arrhenius constants
params.arr.preExpFactor = 1.7/100;     % m3 kg-1 s-1, actually (1.7 * 10^ 3)
params.arr.activationEnergy = 68.5*1000;     % kJ mol-1 -> J mol-1
params.arr.gasConst = 8.314;            % J/(mol·K)

% Energy balance constants
params.eb.enthalpyReaction = 41.2*1000; % kJ/mol -> J/mol
    
% Carbon Dioxide
params.eb.CO2.Fin = ((725.81+6532.28)*1000)/(60*60*24); % kmol/day -> mol/s
params.eb.CO2.A = 25.0;
params.eb.CO2.B = 55.2;
params.eb.CO2.C = -33.7;
params.eb.CO2.D = 7.9;
params.eb.CO2.E = -0.1;
params.CO2.Hf = -393.51*1000; % KJ/mol -> J/mol
params.CO2.S = 213.79;

% Hydrogen 
params.eb.H2.Fin = 22951.16*1000/(60*60*24); % kmol/day -> mol/s
params.eb.H2.A = 18.6; % 1000-2500
params.eb.H2.B = 12.3;
params.eb.H2.C = -2.9;
params.eb.H2.D = 0.3;
params.eb.H2.E = 2.0;
% params.eb.H2.A = 33.066; % these are 298-1000K
% params.eb.H2.B = -11.36;
% params.eb.H2.C = 11.433;
% params.eb.H2.D = -2.77;
% params.eb.H2.E = -0.159;
params.H2.Hf = 0;
params.H2.S = 130.68;

% Carbon monoxide
params.eb.CO.Fin = 7381.52*1000/(60*60*24); % mol/s
params.eb.CO.A = 25.6;
params.eb.CO.B = 6.1;
params.eb.CO.C = 4.1;
params.eb.CO.D = -2.7;
params.eb.CO.E = 0.1;
params.CO.Hf = -110.53*1000; % kJ/mol -> J/mol
params.CO.S = 197.66;

% Methane
params.eb.CH4.Fin = params.inlet.CH4; % mol/s
params.eb.CH4.A = -0.7; 
params.eb.CH4.B = 108.5;
params.eb.CH4.C = -42.5;
params.eb.CH4.D = 5.9;
params.eb.CH4.E = 0.7;

% Water
params.H2O.Hf = -241.83*1000; % kJ/mol -> J/mol
params.H2O.S = 188.84; 
%%
% Ergun equation
params.reactor.diameter = 0.75; % m
params.ergun.bulkDensity = 1200; % kg/m3
params.ergun.particleDensity = 1910; % kg/m3
params.ergun.voidage = (params.ergun.particleDensity-params.ergun.bulkDensity)/params.ergun.particleDensity; % -
params.ergun.particleDiamater = 0.006; % m
params.ergun.csArea = (pi*params.reactor.diameter^2)/4; % m2
params.ergun.initialTotalMolarFlow = params.eb.CO2.Fin+params.eb.H2.Fin + params.eb.CO.Fin + params.inlet.CH4 + params.inlet.gases; % (mol/s)

% Initial Density Calculator
%   Mole fractions
params.inlet.molFrac.CO2 = params.eb.CO2.Fin/params.ergun.initialTotalMolarFlow;
params.inlet.molFrac.H2 = params.eb.H2.Fin/params.ergun.initialTotalMolarFlow;
params.inlet.molFrac.CO = params.eb.CO.Fin/params.ergun.initialTotalMolarFlow;
params.inlet.molFrac.CH4 = params.inlet.CH4/params.ergun.initialTotalMolarFlow;
params.inlet.molFrac.gases = params.inlet.gases/params.ergun.initialTotalMolarFlow;

%   Molar Masses gmol-1
params.molMass.CO2 = 44;
params.molMass.H2 = 2.016;
params.molMass.CO = 28.01;
params.molMass.CH4 = 16.04;
params.molMass.gases = 35;

params.ergun.inletDensity = densityCalculation(params);

% Volumetric flowrate calculation
params.inlet.totalMassFlowrate = (340236.1+287420.26)/(60*60*24); % kgday-1 -> kg s-1
params.inlet.totalVolFlowrate = params.inlet.totalMassFlowrate/params.ergun.inletDensity; % m3 s-1
params.ergun.supVel = params.inlet.totalVolFlowrate/params.ergun.csArea; % m s-1
disp(params.ergun.supVel)
% Gas Flux
params.ergun.gasFlux = params.ergun.inletDensity*params.ergun.supVel; % kg m-2 s-1

% Viscocity
params.ergun.mixtureViscocity = viscocityCalculation(params);
disp(params.ergun.mixtureViscocity)

% Display results in a table
ParamNames = {'Reactor Diameter'; 'Particle Diameter'; 'Bed Voidage'; 'Cross-Sectional Area'; ...
              'Initial Total Molar Flow'; 'Inlet Density'; 'Total Mass Flowrate'; ...
              'Total Volumetric Flowrate'; 'Superficial Velocity'; 'Gas Flux'};

Values = [params.reactor.diameter; params.ergun.particleDiamater; params.ergun.voidage; params.ergun.csArea; ...
          params.ergun.initialTotalMolarFlow; params.ergun.inletDensity; params.inlet.totalMassFlowrate; ...
          params.inlet.totalVolFlowrate; params.ergun.supVel; params.ergun.gasFlux];

Units = {'m'; 'm'; '-'; 'm^2'; 'mol/s'; 'kg/m^3'; 'kg/s'; 'm^3/s'; 'm/s'; 'kg/m^2·s'};

ResultsTable = table(ParamNames, Values, Units);
disp(ResultsTable);
%%
params.cpCO2 = schomate(params, params.inlet.temp / 1000, 'CO2'); % Convert to kK
params.cpH2 = schomate(params, params.inlet.temp / 1000, 'H2');
params.cpCO = schomate(params, params.inlet.temp / 1000, 'CO');
params.cpCH4 = schomate(params, params.inlet.temp / 1000, 'CH4');

% Define initial conditions 
FA0 = params.eb.CO2.Fin; 
FB0 = params.eb.H2.Fin; 
FC0 = params.eb.CO.Fin;
FD0 = 0; 
T0 = params.inlet.temp; % K
P0 = params.inlet.pres; % Pa

% Assign initial conditions to vector
Y0 = [FA0, FB0, FC0, FD0, T0, P0];

Wspan = [0 5000];

% Pass params to odeSolver using an anonymous function
[w,Y] = ode45(@(w,Y) odeSolver(w,Y,params), Wspan, Y0); 
FA = Y(:,1); FB = Y(:,2); FC = Y(:,3); FD = Y(:,4); T = Y(:,5); P = Y(:,6);

% Conversion CO2-> in-out/in
% Calculate conversion values
conversionCO2 = zeros(length(FA), 1);
for i = 1:length(FA)
    conversionCO2(i) = (params.eb.CO2.Fin - FA(i)) / (params.eb.CO2.Fin);
end

% Calculate reactor length
reactorLength = zeros(length(w), 1);
for i = 1:length(w)
    reactorLength(i) = (4*w(i)/(params.ergun.bulkDensity*pi*(params.reactor.diameter^3)));
end 

% Plot CO2 conversion vs. reactor length
figure;
subplot(3,2,1)
plot(reactorLength, conversionCO2, 'b', 'LineWidth', 1.5);
xlabel('Reactor Length (m)');
ylabel('CO_2 Conversion');
title('CO_2 Conversion vs. Reactor Length (m)');
grid on;

% Plot CO2 conversion vs. Weight
subplot(3,2,2);
plot(w, conversionCO2, 'b', 'LineWidth', 1.5);
xlabel('Catalyst Weight (kg)');
ylabel('CO_2 Conversion');
title('CO_2 Conversion vs. Temperature');
grid on;

% Plot all variables
subplot(3,2,3);
plot(w, FA, 'r', w, FB, 'b', w, FC, 'g', w, FD, 'm', 'LineWidth', 1.5);
xlabel('Weight of Catalyst (kg)');
ylabel('Molar Flow Rate (mol/s)');
legend('CO2 (FA)', 'H2 (FB)', 'CO (FC)', 'H2O (FD)');
title('Molar Flow Rates vs. Catalyst Weight');
grid on;

subplot(3,2,4);
plot(w, T, 'k', 'LineWidth', 1.5);
xlabel('Weight of Catalyst (kg)');
ylabel('Temperature (K)');
title('Temperature vs. Catalyst Weight');
grid on;

figure
plot(w, P, 'k', 'LineWidth', 1.5);
xlabel('Weight of Catalyst (kg)');
ylabel('Pressure (Pa)');
title('Pressure vs. Catalyst Weight');
grid on;

%%
% ODE Solver Function
function dYdt = odeSolver(w,Y,params) %#ok<INUSD> 

% Extract state variables
FA = Y(1); FB = Y(2); FC = Y(3); FD = Y(4); T = Y(5); P = Y(6);

% Rate constant calculations using the Arrhenius equation
k = params.arr.preExpFactor * exp(-params.arr.activationEnergy / (params.arr.gasConst * T));

% Equilibrium calculations
deltaHf = (params.CO.Hf+params.H2O.Hf)-(params.CO2.Hf+params.H2.Hf);
deltaS = (params.CO.S+params.H2O.S)-(params.CO2.S+params.H2.S);
deltaG = deltaHf -(T*deltaS);
Keq = exp(-deltaG/(params.arr.gasConst*T));

% Mole fraction calculations 
totalMol = FA + FB + FC + FD + params.inlet.CH4 + params.inlet.gases;
molFractionCO2 = FA / totalMol;
molFractionH2 = FB / totalMol;
molFractionCO = FC / totalMol;
molFractionH2O = FD / totalMol;

% Assume pressure P = 1 atm (if needed, define it properly)
% P = 22*100000; % Placeholder before pressure equation determined

% Rate of reaction calculations
% rRWGS = k * molFractionCO2 * molFractionH2 * ((P^2) / (params.arr.gasConst^2) * (T^2)); 
rRWGS = (k*P^2/(params.arr.gasConst^2*T^2))*((molFractionCO2*molFractionH2)-((molFractionCO*molFractionH2O)/Keq));

sumNcp = params.cpCO2*params.eb.CO2.Fin + params.cpH2*params.eb.H2.Fin + params.cpCO*params.eb.CO.Fin + params.cpCH4*params.eb.CH4.Fin;

beta = ((params.ergun.gasFlux*(1-params.ergun.voidage))/(params.ergun.inletDensity*params.ergun.particleDiamater*params.ergun.voidage^3))*((150*(1-params.ergun.voidage)*params.ergun.mixtureViscocity)/params.ergun.particleDiamater)+(1.75*params.ergun.gasFlux);

% Mole balance ODEs
dFA_dw = -rRWGS;
dFB_dw = -rRWGS;
dFC_dw = rRWGS;
dFD_dw = rRWGS;

dT_dw = (params.eb.enthalpyReaction * rRWGS) / (sumNcp);

dP_dw = (-beta/(params.ergun.csArea*(1-params.ergun.voidage)*params.ergun.particleDensity))*(params.inlet.pres/P)*(T/params.inlet.temp)*(totalMol/params.ergun.initialTotalMolarFlow);

% dP_dw = (1/(params.ergun.bulkDensity*params.ergun.csArea))*(((150*((1-params.ergun.voidage)^2)*params.ergun.mixtureViscocity*params.ergun.supVel)/((params.ergun.voidage^3)*(params.ergun.particleDiamater^2)))+((1.75*(1-params.ergun.voidage)*params.ergun.inletDensity*(params.ergun.supVel^2))/((params.ergun.voidage^3)*params.ergun.particleDiamater)));

% Page 169 elements Fogler

% Output vector for ODE solver
dYdt = [dFA_dw; dFB_dw; dFC_dw; dFD_dw; dT_dw; dP_dw];

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
%%
function inletDensity = densityCalculation(params)

    % Weighted average molar mass
    averageInletMolarMass = ((params.inlet.molFrac.CO2*params.molMass.CO2)+...
                            (params.inlet.molFrac.H2*params.molMass.H2) + ...
                            (params.inlet.molFrac.CO*params.molMass.CO) + ...
                            (params.inlet.molFrac.CH4*params.molMass.CH4) + ...
                            (params.inlet.molFrac.gases*params.molMass.gases))/1000; % gmol-1 -> kgmol-1

    % Gas inlet density calculation
    inletDensity = (params.inlet.pres*averageInletMolarMass)/(params.arr.gasConst*params.inlet.temp);

end
%%
function mixtureViscocity = viscocityCalculation(params)
    
    y = [params.inlet.molFrac.CO2, params.inlet.molFrac.H2, params.inlet.molFrac.CO, params.inlet.molFrac.CH4, params.inlet.molFrac.gases];
    
    % Dynamic viscosities (Pa.s) at 1000 K (approximate values)
    mu = [41.25e-6, 18.77e-6, 40.4e-6, 38.6e-6, 3.5e-5];
        %   https://www.engineeringtoolbox.com/carbon-dioxide-dynamic-kinematic-viscosity-temperature-pressure-d_2074.html#:~:text=Online%20calculator%2C%20figures%20and%20table%20showing%20dynamic%20and,deformation%20by%20shear%20stress%20or%20tensile%20stress%20.,
        %   https://www.lmnoeng.com/Flow/GasViscosity.php
        %   https://www.lmnoeng.com/Flow/GasViscosity.php
        %   https://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html?
        %   The Viscosity and Thermal Conductivity of Ethane in the Limit of Zero Density 
    
    % Molecular weights (kg/mol*1000)
    M = [params.molMass.CO2*1000, params.molMass.H2*1000, params.molMass.CO*1000, params.molMass.CH4*1000, params.molMass.gases*1000];
    
    n = length(y);
    phi = zeros(n, n);
    
    % Compute interaction coefficients
    for i = 1:n
        for j = 1:n
            if i ~= j
                phi(i,j) = ((1 + (mu(i)/mu(j))^0.5 * (M(j)/M(i))^0.25)^2) / ...
                           (sqrt(8) * (1 + M(i)/M(j)));
            else
                phi(i,j) = 1;
            end
        end
    end
    
    % Compute mixture viscosity using Wilke's equation
    denominator = sum(y .* phi, 2);
    mu_m = sum((y .* mu) ./ denominator);
    
    % Element-wise multiplication of mu_m with y
    mu_m_y = mu_m .* y;
    
    % Sum all the values of mu_m_y into a new variable
    mixtureViscocity = sum(mu_m_y);

end