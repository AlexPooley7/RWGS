% CE30243 - Individual Design Project
% Description - Models the transient behaviour of the reverse water-gas
% shift reaction (RWGS)
% Last edited: 28/03/2025
% Last commit: 04/04/2025
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
global thieleModLog effFactorLog %#ok<GVMIS> 
% Define constants:
params = struct(); % Initialise the params structure to hold all constants

% Misc inlet parameters
params.inlet.H2O = 0*1000/(60*60*24); % kmol/day -> mol/s
params.inlet.CH4 = 610.43*1000/(60*60*24); % kmol/day -> mol/s
params.inlet.gases = 1295.90*1000/(60*60*24); % kmol/day -> mol/s
params.inlet.temp = 1000; % K
params.inlet.pres = 22*100000; % bar -> Pa

% Arrhenius constants
params.arr.preExpFactor = 1.7 * 10^ -2; % m3 kg-1 s-1, actually (1.7 * 10^ 3)
params.arr.activationEnergy = 68.5*1000; % kJ mol-1 -> J mol^-1
params.arr.gasConst = 8.314; % J/(mol K)

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
params.CO2.S = 213.79; % J/molK

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
params.H2.S = 130.68; % J/molK

% Carbon monoxide
params.eb.CO.Fin = 7381.52*1000/(60*60*24); % mol/s
params.eb.CO.A = 25.6;
params.eb.CO.B = 6.1;
params.eb.CO.C = 4.1;
params.eb.CO.D = -2.7;
params.eb.CO.E = 0.1;
params.CO.Hf = -110.53*1000; % kJ/mol -> J/mol
params.CO.S = 197.66; % J/molK

% Methane
params.eb.CH4.Fin = params.inlet.CH4; % mol/s
params.eb.CH4.A = -0.7; 
params.eb.CH4.B = 108.5;
params.eb.CH4.C = -42.5;
params.eb.CH4.D = 5.9;
params.eb.CH4.E = 0.7;

% Water
params.H2O.Hf = -241.83*1000; % kJ/mol -> J/mol
params.H2O.S = 188.84; % J/molK

% Heat capacities
params.cpCO2 = schomate(params, params.inlet.temp / 1000, 'CO2'); % Convert to kK
params.cpH2 = schomate(params, params.inlet.temp / 1000, 'H2');
params.cpCO = schomate(params, params.inlet.temp / 1000, 'CO');
params.cpCH4 = schomate(params, params.inlet.temp / 1000, 'CH4');

%%
% Ergun equation
params.reactor.diameter = 0.75; % m
params.ergun.bulkDensity = 1200; % kg/m3
params.ergun.particleDensity = 1910; % kg/m3
params.ergun.voidage = (params.ergun.particleDensity-params.ergun.bulkDensity)/params.ergun.particleDensity; % -
params.ergun.particleDiameter = 0.006; % m
params.ergun.csArea = (pi*params.reactor.diameter^2)/4; % m2
params.ergun.initialTotalMolarFlow = params.eb.CO2.Fin+params.eb.H2.Fin + params.eb.CO.Fin + params.inlet.CH4 + params.inlet.gases; % mol/s

% Initial Density Calculator
% Mole fractions
params.inlet.molFrac.CO2 = params.eb.CO2.Fin/params.ergun.initialTotalMolarFlow; % -
params.inlet.molFrac.H2 = params.eb.H2.Fin/params.ergun.initialTotalMolarFlow; % -
params.inlet.molFrac.CO = params.eb.CO.Fin/params.ergun.initialTotalMolarFlow; % -
params.inlet.molFrac.CH4 = params.inlet.CH4/params.ergun.initialTotalMolarFlow; % -
params.inlet.molFrac.gases = params.inlet.gases/params.ergun.initialTotalMolarFlow; % -
% Molar Masses 
params.molMass.CO2 = 44; % g/mol
params.molMass.H2 = 2.016; % g/mol
params.molMass.CO = 28.01; % g/mol
params.molMass.CH4 = 16.04; % g/mol
params.molMass.gases = 35; % g/mol
% Density
params.ergun.inletDensity = densityCalculation(params); % kg/m3

% Volumetric flowrate calculation
params.inlet.totalMassFlowrate = (340236.1+287420.26)/(60*60*24); % kgday-1 -> kg s-1
params.inlet.totalVolFlowrate = params.inlet.totalMassFlowrate/params.ergun.inletDensity; % m3 s-1
params.ergun.supVel = params.inlet.totalVolFlowrate/params.ergun.csArea; % m s-1
disp(params.ergun.supVel)

% Gas Flux
params.ergun.gasFlux = params.ergun.inletDensity*params.ergun.supVel; % kg m-2 s-1

% Viscocity
params.ergun.mixtureViscocity = viscocityCalculation(params);

%%
% Thiele Modulus
params.theile.charLength = params.ergun.particleDiameter/6;
params.thiele.porosityParticle = 0.4; % Table B-1
params.thiele.tortuosityParticle = 2.8; % Tbale B-1
params.thiele.molDiff = ; % From DATABASE FIND THIS!!!!
params.thiele.poreDiameter = ; % FIND
params.thiele.knudsenDiff = (params.thiele.poreDiameter/3)*sqrt((8*params.arr.gasConst*T)/(pi*params.molMass.CO2)); % Same unit of molecular diff

thieleModLog = []; 
effFactorLog = []; 
%%
% Define initial conditions 
init.FA0 = params.eb.CO2.Fin; 
init.FB0 = params.eb.H2.Fin; 
init.FC0 = params.eb.CO.Fin;
init.FD0 = 0; 
init.T0 = params.inlet.temp; % K
init.P0 = params.inlet.pres; % Pa

% Assign initial conditions to vector
Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, init.T0, init.P0];

% Set span of catalyst mass to integrate over
Wspan = [0 5000];

% Pass params to odeSolver using odeSolver function
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

% Plot data and call optimisation functions
% plotOriginal(reactorLength,conversionCO2,w,FA,FB,FC,FD,T,P)
% displayTable(params)
% plotThieleEff(thieleModLog, effFactorLog, T)
% optimiseTemp(init,params) % Temperature optimisation function
% optimisePressure(init, params) % Pressure optimisation function
% optimiseParticleDiamater(init, params) % Particle diameter optimisation function
optimiseSuperficialVelocity(init, params)

%%
function dYdt = odeSolver(w,Y,params) %#ok<INUSD> 
    % ODE Solver Function

    global thieleModLog effFactorLog %#ok<GVMIS> 

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

% Rate of reaction calculations
rRWGS = (k*P^2/(params.arr.gasConst^2*T^2))*((molFractionCO2*molFractionH2)-((molFractionCO*molFractionH2O)/Keq));

% Eb denominator
sumNcp = params.cpCO2*params.eb.CO2.Fin + params.cpH2*params.eb.H2.Fin + params.cpCO*params.eb.CO.Fin + params.cpCH4*params.eb.CH4.Fin;

% Beta constant for Ergun equation
beta = ((params.ergun.gasFlux*(1-params.ergun.voidage))/(params.ergun.inletDensity*params.ergun.particleDiameter*params.ergun.voidage^3))*((150*(1-params.ergun.voidage)*params.ergun.mixtureViscocity)/params.ergun.particleDiameter)+(1.75*params.ergun.gasFlux);

% Call the thiele modulus function
[thieleMod,effFactor] = thieleModulus(params,T,k);
thieleModLog(end+1,:) = [w, thieleMod];
effFactorLog(end+1,:) = [w,effFactor];

% Mole balance ODEs
dFA_dw = -rRWGS;
dFB_dw = -rRWGS;
dFC_dw = rRWGS;
dFD_dw = rRWGS;

% Temperature ODE
dT_dw = (params.eb.enthalpyReaction * rRWGS) / (sumNcp);

% Pressure ODE
dP_dw = (-beta/(params.ergun.csArea*(1-params.ergun.voidage)*params.ergun.particleDensity))*(params.inlet.pres/P)*(T/params.inlet.temp)*(totalMol/params.ergun.initialTotalMolarFlow);

% Output vector for ODE solver
dYdt = [dFA_dw; dFB_dw; dFC_dw; dFD_dw; dT_dw; dP_dw];

end

%% Schomate Equation Function
function cp = schomate(params,T,component)
    % Solves Cp for a given temperature, depending on case given

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
    % Plot all heat capacities, to ensure equation works. Only used for
    % debugging

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
    % Calucaltes the weighted inlet density

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
    % Calculates the viscocity of the gas, using Wilkes equation
    
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
%%

function plotOriginal(reactorLength,conversionCO2,w,FA,FB,FC,FD,T,P)
    % Plot the pre optimisation reactor

    % Plot CO2 conversion vs. reactor length
    figure(1);
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
    
    figure(2);
    plot(w, P, 'k', 'LineWidth', 1.5);
    xlabel('Weight of Catalyst (kg)');
    ylabel('Pressure (Pa)');
    title('Pressure vs. Catalyst Weight');
    grid on;
    hold off
end
%%
function [thieleMod,effFactor] = thieleModulus(params,T,k)
    % Estimates the thiele modulus for a given set of reaction rates, taken
    % when the ODE solver is run. 

    knudsenDiff = (params.thiele.poreDiameter/3)*sqrt((8*params.arr.gasConst*T)/(pi*params.molMass.CO2)); % Same unit of molecular diff
    poreDiff = ((1/params.thiele.molDiff)+(1/knudsenDiff))^-1;
    effectiveDiff = (params.thiele.porosityParticle.params.thiele.tortuosityParticle)*poreDiff; 

    thieleMod = params.thiele.charLength*sqrt(k/effectiveDiff);

    effFactor = tanh(thieleMod)/thieleMod;
end

%%
function plotThieleEff(thieleModLog, effFactorLog, T)
    % Plot t against theile mod and effectiveness factor
end

%%
function optimiseTemp(init, params)
    % This function runs the ODE with 10 different temperatures, 
    % 5 above and 5 below 1000K, in 50K increments.
    
    % Define temperature range (5 below and 5 above 1000K in 50K increments)
    temps = 1000 + (-5:5) * 50;

    % Store results
    W_all = cell(1, length(temps)); % Store W values
    T_all = cell(1, length(temps)); % Store Temperature profiles
    FA_all = cell(1, length(temps)); % Store FA values
    CO2ConversionAll = cell(1, length(temps)); % Store CO2 conversion

    % Prepare figures before loop
    figure(3); % For CO2 conversion profiles

    % Loop over each temperature
    for i = 1:length(temps)
        T0 = temps(i);
        
        % Assign initial conditions to vector
        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, T0, init.P0];
        
        % Define span for integration
        Wspan = [0 5000];

        % Solve ODE
        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params), Wspan, Y0);

        % Store data
        W_all{i} = W; 
        T_all{i} = Y(:,5); % Extract Temperature
        FA_all{i} = Y(:,1); % Extract FA
        
        % Compute CO2 conversion
        CO2ConversionAll{i} = (params.eb.CO2.Fin - FA_all{i}) / params.eb.CO2.Fin;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(3);
        plot(W, CO2ConversionAll{i}, 'DisplayName', sprintf('%dK', temps(i)));
        hold on;
    end 

    % Formatting for CO2 Conversion Profile (after loop)
    figure(3);
    xlabel('W (kg catalyst)');
    ylabel('CO2 Conversion');
    title('CO2 Conversion vs. Catalyst Weight for Different Initial Temperatures');
    legend show;
    grid on;
    hold off;
end

%%
function optimisePressure(init, params)
    % This function runs the ODE with 10 different pressure, 
    % 5 above and 5 below 22bar, in 1bar increments - (In pascals)
    
    % Define temperature range (5 below and 5 above 1000K in 50K increments)
    pressure = 22*100000 + (-5:5) * 100000;

    % Store results
    W_all = cell(1, length(pressure)); % Store W values
    P_all = cell(1, length(pressure)); % Store Temperature profiles
    FA_all = cell(1, length(pressure)); % Store FA values
    CO2ConversionAll = cell(1, length(pressure)); % Store CO2 conversion

    % Prepare figures before loop
    figure(4); % For CO2 conversion profiles
    figure(5);

    % Loop over each temperature
    for i = 1:length(pressure)
        P0 = pressure(i);
        
        % Assign initial conditions to vector
        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, init.T0, P0];
        
        % Define span for integration
        Wspan = [0 5000];

        % Solve ODE
        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params), Wspan, Y0);

        % Store data
        W_all{i} = W; 
        P_all{i} = Y(:,6); % Extract Pressure
        FA_all{i} = Y(:,1); % Extract FA
        
        % Compute CO2 conversion
        CO2ConversionAll{i} = (params.eb.CO2.Fin - FA_all{i}) / params.eb.CO2.Fin;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(4);
        plot(W, CO2ConversionAll{i}, 'DisplayName', sprintf('%dPa', pressure(i)));
        hold on;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(5);
        plot(W, P_all{i}, 'DisplayName', sprintf('%dPa', pressure(i)));
        hold on;
    end 

    % Formatting for CO2 Conversion Profile (after loop)
    figure(4);
    xlabel('W (kg catalyst)');
    ylabel('CO2 Conversion');
    title('CO2 Conversion vs. Catalyst Weight for Different Initial Pressures');
    legend show;
    grid on;
    hold off;
    
    % Formatting for CO2 Conversion Profile (after loop)
    figure(5);
    xlabel('W (kg catalyst)');
    ylabel('Pressure (Pa)');
    title('Pressure vs. Catalyst Weight for Different Initial Pressures');
    legend show;
    grid on;
    hold off;
end

%%
function optimiseParticleDiamater(init, params)
    % This function runs the ODE with 10 different catalyst diamaters.
    
    diameters = params.ergun.particleDiameter + (-5:5) * 0.0005;

    % Store results
    W_all = cell(1, length(diameters)); % Store W values
    P_all = cell(1, length(diameters)); % Store Temperature profiles
    FA_all = cell(1, length(diameters)); % Store FA values
    CO2ConversionAll = cell(1, length(diameters)); % Store CO2 conversion

    % Prepare figures before loop
    figure(6); % For CO2 conversion profiles
    figure(7);

    % Loop over each temperature
    for i = 1:length(diameters)
        params.ergun.particleDiameter = diameters(i);
        
        % Assign initial conditions to vector
        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, init.T0, init.P0];
        
        % Define span for integration
        Wspan = [0 5000];

        % Solve ODE
        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params), Wspan, Y0);

        % Store data
        W_all{i} = W; 
        P_all{i} = Y(:,6); % Extract Pressure
        FA_all{i} = Y(:,1); % Extract FA
        
        % Compute CO2 conversion
        CO2ConversionAll{i} = (params.eb.CO2.Fin - FA_all{i}) / params.eb.CO2.Fin;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(6);
        plot(W, CO2ConversionAll{i}, 'DisplayName', sprintf('%dm', diameters(i)));
        hold on;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(7);
        plot(W, P_all{i}, 'DisplayName', sprintf('%dm', diameters(i)));
        hold on;
    end 

    % Formatting for CO2 Conversion Profile (after loop)
    figure(6);
    xlabel('W (kg catalyst)');
    ylabel('CO2 Conversion');
    title('CO2 Conversion vs. Catalyst Weight for Different Catalyst Diameters');
    legend show;
    grid on;
    hold off;
    
    % Formatting for CO2 Conversion Profile (after loop)
    figure(7);
    xlabel('W (kg catalyst)');
    ylabel('Pressure (Pa)');
    title('Pressure vs. Catalyst Weight for Different Catalyst Diameters');
    legend show;
    grid on;
    hold off;
end

%%
function optimiseSuperficialVelocity(init, params)
    % This function runs the ODE with 10 different catalyst diamaters.
    
    supVels = 4 + (-5:5) * 0.5;

    % Store results
    W_all = cell(1, length(supVels)); % Store W values
    P_all = cell(1, length(supVels)); % Store Temperature profiles
    FA_all = cell(1, length(supVels)); % Store FA values
    CO2ConversionAll = cell(1, length(supVels)); % Store CO2 conversion

    % Prepare figures before loop
    figure(8); % For CO2 conversion profiles
    figure(9);

    % Loop over each temperature
    for i = 1:length(supVels)
        params.ergun.supVel = supVels(i);
        params.ergun.csArea = (params.inlet.totalVolFlowrate/params.ergun.supVel);
        params.reactor.diameter = sqrt((4*params.ergun.csArea)/pi);
        params.ergun.gasFlux = params.ergun.inletDensity*params.ergun.supVel;
        disp(params.ergun.gasFlux)
        % Assign initial conditions to vector
        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, init.T0, init.P0];
        
        % Define span for integration
        Wspan = [0 2000];

        % Solve ODE
        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params), Wspan, Y0);

        % Store data
        W_all{i} = W; 
        P_all{i} = Y(:,6); % Extract Pressure
        FA_all{i} = Y(:,1); % Extract FA
        
        % Compute CO2 conversion
        CO2ConversionAll{i} = (params.eb.CO2.Fin - FA_all{i}) / params.eb.CO2.Fin;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(8);
        plot(W, CO2ConversionAll{i}, 'DisplayName', sprintf('%dm/s', supVels(i)));
        hold on;

        % Plot W vs. CO2 Conversion (only once outside the loop for figure 3)
        figure(9);
        plot(W, P_all{i}, 'DisplayName', sprintf('%dm/s', supVels(i)));
        hold on;
    end 

    % Formatting for CO2 Conversion Profile (after loop)
    figure(8);
    xlabel('W (kg catalyst)');
    ylabel('CO2 Conversion');
    title('CO2 Conversion vs. Catalyst Weight for Different Superficial Velocities');
    legend show;
    grid on;
    hold off;
    
    % Formatting for CO2 Conversion Profile (after loop)
    figure(9);
    xlabel('W (kg catalyst)');
    ylabel('Pressure (Pa)');
    title('Pressure vs. Catalyst Weight for Different Superficial Velocities');
    legend show;
    grid on;
    hold off;
end

%%
function displayTable(params)
    % Display results in a table
    ParamNames = {'Reactor Diameter'; 'Particle Diameter'; 'Bed Voidage'; 'Cross-Sectional Area'; ...
                  'Initial Total Molar Flow'; 'Inlet Density'; 'Total Mass Flowrate'; ...
                  'Total Volumetric Flowrate'; 'Superficial Velocity'; 'Gas Flux'};
    
    Values = [params.reactor.diameter; params.ergun.particleDiameter; params.ergun.voidage; params.ergun.csArea; ...
              params.ergun.initialTotalMolarFlow; params.ergun.inletDensity; params.inlet.totalMassFlowrate; ...
              params.inlet.totalVolFlowrate; params.ergun.supVel; params.ergun.gasFlux];
    
    Units = {'m'; 'm'; '-'; 'm^2'; 'mol/s'; 'kg/m^3'; 'kg/s'; 'm^3/s'; 'm/s'; 'kg/m^2·s'};
    
    ResultsTable = table(ParamNames, Values, Units);
    disp(ResultsTable);
end

