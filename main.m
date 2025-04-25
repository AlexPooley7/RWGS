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
close all
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
params.arr.preExpFactor = 1.7 * 10^3; % m3 kg-1 s-1
params.arr.activationEnergy = 68.5*1000; % kJ mol-1 -> J mol^-1
params.arr.gasConst = 8.314; % J/molK

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
params.eb.H2O.A = 30.09; 
params.eb.H2O.B = 6.8325;
params.eb.H2O.C = 6.7934;
params.eb.H2O.D = -2.5345;
params.eb.H2O.E = 0.0821;
params.H2O.Hf = -241.83*1000; % kJ/mol -> J/mol
params.H2O.S = 188.84; % J/molK

% Langmuir Hinshelwood
params.LH.adsH2 = 6.12*10^-4 * 101325; % Pa-1
params.LH.adsCO2 = 8.23*10^-5 * 101325; % Same as CO Pa-1
params.LH.adsCO = 8.23*10^-5 * 101325; % Pa-1
params.LH.adsH2O = 1.77*10^-5 * 101325; % Pa-1

% Heat capacities
params.cpCO2 = schomate(params, params.inlet.temp / 1000, 'CO2'); % Convert to kK
params.cpH2 = schomate(params, params.inlet.temp / 1000, 'H2');
params.cpCO = schomate(params, params.inlet.temp / 1000, 'CO');
params.cpCH4 = schomate(params, params.inlet.temp / 1000, 'CH4');

% Ergun equation
params.reactor.diameter = 1; % Choose
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

% Gas Flux
params.ergun.gasFlux = params.ergun.inletDensity*params.ergun.supVel; % kg m-2 s-1

% Viscocity
params.ergun.mixtureViscocity = viscocityCalculation(params); % Pa s

% Thiele Modulus and Effectiveness Factor
params.thiele.charLength = params.ergun.particleDiameter/2; % m
params.thiele.porosityParticle = 0.4; 
params.thiele.tortuosityParticle = 1.6; 
params.thiele.molDiff = 1.83*10^-4; % m2/s
params.thiele.poreDiameter = 2.3*10^-7; % m

% Initialise the thieleMod and Eff factor structures
thieleModLog = []; 
effFactorLog = []; 

% Initialise the shortcut design parameters
params.shortcut.pressureInit = params.inlet.pres;
params.shortcut.temperatureInit = params.inlet.temp;

%% === Shortcut ===
% Initialise the ode initial conditions
init.FA0 = params.eb.CO2.Fin; 
init.FB0 = params.eb.H2.Fin; 
init.FC0 = params.eb.CO.Fin;
init.FD0 = 0; 

% Assign initial conditions to vector
Y0 = [init.FA0, init.FB0, init.FC0, init.FD0];

% Set span of catalyst mass to integrate over
Wspan = [0 3000];

% Pass params to odeSolver using odeSolver function
[wS,Y] = ode45(@(w,Y) odeSolverShortcut(w,Y,params), Wspan, Y0); 
FAS = Y(:,1); FBS = Y(:,2); FCS = Y(:,3); FDS = Y(:,4);

% Calculate conversion of CO2
conversionCO2S = conversion(params,FAS);

% Calculate reactor length
reactorLengthS = reactorLengthFunc(params,wS);

% Plot data
% plotOriginalShortcut(reactorLengthS,conversionCO2S,wS,FAS,FBS,FCS,FDS)

%% === Rigorous ===
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
Wspan = [0 10000];

% Set the condition for Power or LH kinetics
%% Power Law
condition = 'P'; % Power law

% Pass params to odeSolver using odeSolver function
[wp,Y] = ode45(@(w,Y) odeSolver(w,Y,params, condition), Wspan, Y0); 
FAp = Y(:,1); FBp = Y(:,2); FCp = Y(:,3); FDp = Y(:,4); Tp = Y(:,5); Pp = Y(:,6);

% Calculate conversion
conversionCO2p = conversion(params,FAp);

% Calculate reactor length
reactorLengthp = reactorLengthFunc(params,wp);

%% Langmuir-Hinshelwood
condition = 'L-H'; % Langmuir-Hinshelwood

% Pass params to odeSolver using odeSolver function
[wLH,Y] = ode45(@(w,Y) odeSolver(w,Y,params, condition), Wspan, Y0); 
FALH = Y(:,1); FBLH = Y(:,2); FCLH = Y(:,3); FDLH = Y(:,4); TLH = Y(:,5); PLH = Y(:,6);

% Calculate conversion
conversionCO2LH = conversion(params,FALH);

% Calculate reactor length
reactorLengthLH = reactorLengthFunc(params,wLH);

%%  Plot data and call optimisation functions
% plotOriginal(reactorLengthp,conversionCO2p,wp,FAp,FBp,FCp,FDp,Tp,Pp,reactorLengthLH,conversionCO2LH,wLH,FALH,FBLH,FCLH,FDLH,TLH,PLH)
% displayTable(params)
% plotThieleEff(thieleModLog, effFactorLog, wLH)
optimiseTemp(init,params) % Temperature optimisation function
% optimisePressure(init, params) % Pressure optimisation function
% optimiseParticleDiamater(init, params) % Particle diameter optimisation function
% optimiseSuperficialVelocity(init, params)
% optimiseParticleAndBedDiameter(init, params)
% plotThieleEff(thieleModLog, effFactorLog, wLH)

%%
function dYdt = odeSolverShortcut(w,Y,params) 
    % ODE Solver Function

    % Extract state variables
    FA = Y(1); FB = Y(2); FC = Y(3); FD = Y(4);
    
    T = params.shortcut.temperatureInit;
    % Rate constant calculations using the Arrhenius equation
    k = (params.arr.preExpFactor*10^-5)*exp(-(params.arr.activationEnergy/(params.arr.gasConst*T)));
    
    % Mole fraction calculations 
    totalMol = FA + FB + FC + FD + params.inlet.CH4 + params.inlet.gases;
    molFractionCO2 = FA / totalMol;
    molFractionH2 = FB / totalMol;
    molFractionCO = FC / totalMol;
    molFractionH2O = FD / totalMol;
    
    P = params.shortcut.pressureInit;

    % Rate of reaction calculations
    rRWGS = ((k*(P^2))/((params.arr.gasConst^2)*(T^2)))*(molFractionCO2*molFractionH2);
               
    % Mole balance ODEs
    dFA_dw = -rRWGS;
    dFB_dw = -rRWGS;
    dFC_dw = rRWGS;
    dFD_dw = rRWGS;
       
    % Output vector for ODE solver
    dYdt = [dFA_dw; dFB_dw; dFC_dw; dFD_dw];

end
%%
function dYdt = odeSolver(w,Y,params, condition) 
    % Rigorous ODE Solver Function

    % Global variable statement
    global thieleModLog effFactorLog %#ok<GVMIS> 

    % Extract state variables
    FA = Y(1); FB = Y(2); FC = Y(3); FD = Y(4); T = Y(5); P = Y(6);
    
    % Rate constant calculations using the Arrhenius equation
    k = (params.arr.preExpFactor*10^-5)*exp(-(params.arr.activationEnergy/(params.arr.gasConst*T)));
    
    % Equilibrium calculations
    deltaHf = (params.CO.Hf+params.H2O.Hf)-(params.CO2.Hf+params.H2.Hf);
    deltaS = (params.CO.S+params.H2O.S)-(params.CO2.S+params.H2.S);
    deltaG = deltaHf -(1000*deltaS);
    Keq = exp(-deltaG/(params.arr.gasConst*T));
    
    % Mole fraction calculations 
    totalMol = FA + FB + FC + FD + params.inlet.CH4 + params.inlet.gases;
    molFractionCO2 = FA / totalMol;
    molFractionH2 = FB / totalMol;
    molFractionCO = FC / totalMol;
    molFractionH2O = FD / totalMol;

    % Partial Pressures
    ppCO2 = P*molFractionCO2;
    ppH2 = P*molFractionH2;
    ppCO = P*molFractionCO;
    ppH2O = P*molFractionH2O;

    % Keq conversion
    Kr = Keq/8.314*T;
    cT = 10^2.5;

    if condition == 'P'
        rRWGS = ((k*(P^2))/((params.arr.gasConst^2)*(T^2)))*((molFractionCO2*molFractionH2)-((molFractionCO*molFractionH2O)/Keq));
    else
        kineticTerm = k*params.LH.adsCO2*params.LH.adsH2*cT^2;
        drivingForce = (ppCO2*ppH2)-((ppH2O*ppCO)/Kr);
        adsorptionTerm = (1 + params.LH.adsCO2*ppCO2 + params.LH.adsH2*ppH2 + ppCO/params.LH.adsCO + ppH2O/params.LH.adsH2O)^2;
        rRWGS = (kineticTerm*drivingForce)/adsorptionTerm;
    end
   
    % Eb denominator
    params.cpCO2 = schomate(params, T / 1000, 'CO2'); % Convert to kK
    params.cpH2 = schomate(params, T / 1000, 'H2');
    params.cpCO = schomate(params, T / 1000, 'CO');
    params.cpH2O = schomate(params,T / 1000, 'H2O');
    params.cpCH4 = schomate(params, T / 1000, 'CH4');
    
    % sumNcp = params.cpCO2*params.eb.CO2.Fin + params.cpH2*params.eb.H2.Fin + params.cpCO*params.eb.CO.Fin + params.cpCH4*params.eb.CH4.Fin;
    sumNcp = params.cpCO2*FA + params.cpH2*FB + params.cpCO*FC + params.cpH2O*FD + params.cpCH4*params.eb.CH4.Fin;

    % Beta constant for Ergun equation
    beta = ((params.ergun.gasFlux*(1-params.ergun.voidage))/(params.ergun.inletDensity*params.ergun.particleDiameter*params.ergun.voidage^3))*((150*(1-params.ergun.voidage)*params.ergun.mixtureViscocity)/params.ergun.particleDiameter)+(1.75*params.ergun.gasFlux);
    
    % Thiele Modulus and Effectiveness factor calculations
    knudsenDiff = (params.thiele.poreDiameter/3)*sqrt((8*params.arr.gasConst*T)/(pi*(params.molMass.CO2/1000))); % m2/s
    
    poreDiff = ((1/params.thiele.molDiff)+(1/knudsenDiff))^-1; % m2/s
    effectiveDiff = (params.thiele.porosityParticle*params.thiele.tortuosityParticle)*poreDiff; % m2/s 
    
    thieleMod = sqrt((k*1000*w*(params.ergun.particleDiameter/2)^2)/effectiveDiff);
    
    if thieleMod == 0 || isnan(thieleMod)
        effFactor = 1;  
    else
        effFactor = (1/thieleMod)*((1/tanh(3*thieleMod))-(1/(3*thieleMod)));
    end
    
    % Call the thiele modulus function
    % [thieleMod,effFactor] = thieleModulus(params,T,k,w);
    thieleModLog(end+1) = thieleMod;
    effFactorLog(end+1) = effFactor;
    
    % Mole balance ODEs
    dFA_dw = -rRWGS*effFactor;
    dFB_dw = -rRWGS*effFactor;
    dFC_dw = rRWGS*effFactor;
    dFD_dw = rRWGS*effFactor;
    
    % Temperature ODE
    dT_dw = (-params.eb.enthalpyReaction * rRWGS*effFactor)/(sumNcp);
    
    % Pressure ODE
    dP_dw = (-beta/(params.ergun.csArea*(1-params.ergun.voidage)*params.ergun.particleDensity))*(params.inlet.pres/P)*(T/params.inlet.temp)*(totalMol/params.ergun.initialTotalMolarFlow);
    
    % Output vector for ODE solver
    dYdt = [dFA_dw; dFB_dw; dFC_dw; dFD_dw; dT_dw; dP_dw];

end


%%
function plotOriginal(reactorLengthp,conversionCO2p,wp,FAp,FBp,FCp,FDp,Tp,Pp,reactorLengthLH,conversionCO2LH,wLH,FALH,FBLH,FCLH,FDLH,TLH,PLH)
    % Plot the pre optimisation reactor
    
    % Get exact mass when convrsion is 0.376
    indexp = find(conversionCO2p > 0.376, 1, 'first');
    wpValue = wp(indexp);
    disp(wpValue)
    
    % Langmuir-Hinshelwood
    indexLH = find(conversionCO2LH > 0.376, 1, 'first');
    wLHValue = wLH(indexLH);
    disp(wLHValue)

    % Get exact length when convrsion is 0.376
    lpValue = reactorLengthp(indexp);
    disp(lpValue)
    
    lLHValue = reactorLengthLH(indexLH);
    disp(lLHValue)
    
    % Temperature
    tpValue = Tp(indexp);
    disp(tpValue)

    tLHValue = TLH(indexLH);
    disp(tLHValue)

    % Pressure
    ppValue =Pp(indexp);
    disp(ppValue)

    pLHValue = PLH(indexLH);
    disp(pLHValue)

    % Plot CO2 conversion vs. catalyst weight
    figure (1)
    plot(wp, conversionCO2p, 'g', 'LineWidth', 1.5)
    hold on
    plot(wLH, conversionCO2LH, 'b', 'LineWidth', 1.5)
    xlabel('Catalyst Mass (kg)');
    ylabel('CO2 Conversion');
    yline(0.376, '--r', 'LineWidth', 1.2);  % Add horizontal line
    % Plot vertical line from y=0 to y=0.376 at wpValue
    line([wpValue wpValue], [0 0.376], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    % Plot vertical line from y=0 to y=0.376 at wLHValue
    line([wLHValue wLHValue], [0 0.376], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    grid off
    legend('Power Law', 'Langmuir-Hinshelwood', 'Desired Conversion','Required Catalyst Mass')
    
    % Plot CO2 conversion vs. reactor length
    figure (2)
    plot(reactorLengthp, conversionCO2p, 'g', 'LineWidth', 1.5)
    hold on
    plot(reactorLengthLH, conversionCO2LH, 'b', 'LineWidth', 1.5)
    xlabel('Reactor Length (m)');
    ylabel('CO2 Conversion');
    yline(0.376, '--r', 'LineWidth', 1.2);  % Add horizontal line
    grid off
     % Plot vertical line from y=0 to y=0.376 at wpValue
    line([lpValue lpValue], [0 0.376], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    % Plot vertical line from y=0 to y=0.376 at wLHValue
    line([lLHValue lLHValue], [0 0.376], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    legend('Power Law', 'Langmuir-Hinshelwood', 'Desired Conversion','Required ReactorLength')
    
    % Plot temperature vs. reactor length
    figure (3)
    plot(reactorLengthp,Tp, 'g', 'LineWidth', 1.5)
    hold on
    plot(reactorLengthLH,TLH, 'b', 'LineWidth', 1.5)
    xlabel('Reactor Length (m)');
    ylabel('Temperature (K)');
    grid off
    xlim([0 7])
    % Plot vertical line from y=0 to T at lpValue
    line([lpValue lpValue], [880 tpValue], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    % Plot vertical line from y=0 to y=0.376 at lLHValue
    line([lLHValue lLHValue], [880 tLHValue], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    legend('Power Law', 'Langmuir-Hinshelwood', 'Required Reactor Length')

    % Plot pressure vs. reactor length 
    figure (4)
    plot(reactorLengthp,Pp, 'g', 'LineWidth', 1.5)
    hold on
    plot(reactorLengthLH,PLH, 'b', 'LineWidth', 1.5)
    xlabel('Reactor Length (m)');
    ylabel('Pressure (Pa)');
    grid off
    % Plot vertical line from y=0 to T at lpValue
    line([lpValue lpValue], [2.15e6 ppValue], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    % Plot vertical line from y=0 to y=0.376 at lLHValue
    line([lLHValue lLHValue], [2.15e6 pLHValue], 'Color', 'y', 'LineStyle', '--', 'LineWidth', 1.2);
    legend('Power Law', 'Langmuir-Hinshelwood', 'Required Reactor Length')
    xlim([0 7])

    % Plot all flowrates of respective components power law
    figure(5)
    plot(reactorLengthp, FAp, 'g', reactorLengthp, FBp, 'b', reactorLengthp, FCp, 'y', reactorLengthp, FDp, 'r', 'LineWidth',1.5)
    xlabel('Reactor Length (m)');
    ylabel('Molar Flowrate (mol/s)');
    grid off
    xline([lpValue lpValue], '--k', 'LineWidth', 1.2)
    legend('CO2', 'H2', 'CO', 'H20', 'Required Reactor Length')
    xlim([0 7])
    
    % Plot all flowrates of respective components LH
    figure(6)
    plot(reactorLengthLH, FALH, 'g', reactorLengthLH, FBLH, 'b', reactorLengthLH, FCLH, 'y', reactorLengthLH, FDLH, 'r', 'LineWidth',1.5)
    xlabel('Reactor Length (m)');
    ylabel('Molar Flowrate (mol/s)');
    grid off
    xline([lLHValue lLHValue], '--k', 'LineWidth', 1.2)
    legend('CO2', 'H2', 'CO', 'H20', 'Required Reactor Length')
    xlim([0 7])


end

%% 
function plotThieleEff(thieleModLog, effFactorLog, wLH)
    % Plot Thiele modulus and effectiveness factor against weight of catalyst
    
    % Interpolate Thiele modulus data
    thieleModInterp = interp1(1:length(thieleModLog), thieleModLog, linspace(1, length(thieleModLog), length(wLH)));
    
    % Interpolate effectiveness factor data
    effInterp = interp1(1:length(effFactorLog), effFactorLog, linspace(1, length(effFactorLog), length(wLH)));

    % Perform polynomial regression on the interpolated Thiele modulus data
    pThiele = polyfit(wLH, thieleModInterp, 2);
    
    % Evaluate the polynomial at the original temperature points for Thiele modulus
    thieleModRegressed = polyval(pThiele, wLH);
    
    % Plot the Thiele modulus with regression line
    figure(2);
    plot(wLH, thieleModRegressed, 'r--', 'LineWidth', 2);  % Regression line
    xlabel('Weight [kg]');
    ylabel('Thiele Modulus');
    title('Thiele Modulus vs Weight');
    xlim([0 7000])
    grid off;
    
    % Plot the effectiveness factor vs. Thiele modulus
    figure(3);
    loglog(thieleModLog, effFactorLog, 'g-', 'LineWidth', 2);  % Regression line
    xlabel('Thiele Modulus');
    ylabel('Effectiveness Factor');
    xticks([0.1 0.2 0.4 0.6 1 2 4 6 10 20 30 40]);
    xticklabels({'0.1','0.2','0.4','0.6','1','2','4','6'});
    yticks([0.1 0.2 0.4 0.6 0.8 1]);
    yticklabels({'0.1','0.2','0.4','0.6','0.8','1'});
    xlim([0.1 7]);
    ylim([0 1])
    grid off;
 
end

%%
function plotOriginalShortcut(reactorLength,conversionCO2,w,FA,FB,FC,FD)
    % Plot the pre optimisation reactor

    % Plot CO2 conversion vs. catalyst weight
    figure (1)
    plot(w, conversionCO2, 'b', 'LineWidth', 1.5)
    xlabel('Catalyst Weight (kg)');
    ylabel('CO2 Conversion');
    yline(0.376, '--r', 'LineWidth', 1.2);  % Add horizontal line
    grid off
    legend('Shortcut','Desired Conversion')

    % Plot CO2 conversion vs. reactor length
    figure (2)
    plot(reactorLength, conversionCO2, 'b', 'LineWidth', 1.5)
    xlabel('Reactor Length (m)');
    ylabel('CO2 Conversion');
    yline(0.376, '--r', 'LineWidth', 1.2);  % Add horizontal line
    grid off
    legend('Shortcut', 'Desired Conversion')
    
    % Plot all flowrates of respective components
    figure(3)
    plot(reactorLength, FA, 'g', reactorLength, FB, 'b', reactorLength, FC, 'y', reactorLength, FD, 'r', 'LineWidth',1.5)
    xlabel('Reactor Length (m)');
    ylabel('Molar Flowrate (mol/s)');
    grid off
    legend('CO2', 'H2', 'CO', 'H20')
    xlim([0 3])

end
%%
function optimiseTemp(init, params)
    % Temperature range
    temps = 1000 + (-5:5) * 50;

    % Store results
    wAll = cell(1, length(temps));
    faAll = cell(1, length(temps));
    co2ConversionAll = cell(1, length(temps));
    wAtTargetConversion = nan(1, length(temps));

    targetConversion = 0.376;

    figure(1); % CO2 conversion profiles

    % Loop over each temperature
    for i = 1:length(temps)
        t0 = temps(i);
        y0 = [init.FA0, init.FB0, init.FC0, init.FD0, t0, init.P0];
        wSpan = [0 10000];
        condition = 'LH';

        % Solve ODE
        [w, y] = ode45(@(w, y) odeSolver(w, y, params, condition), wSpan, y0);

        % Store and process
        wAll{i} = w;
        faAll{i} = y(:,1);
        co2Conversion = (params.eb.CO2.Fin - y(:,1)) / params.eb.CO2.Fin;
        co2ConversionAll{i} = co2Conversion;

        % Find catalyst weight for target conversion
        idx = find(co2Conversion >= targetConversion, 1);
        if ~isempty(idx)
            wAtTargetConversion(i) = w(idx);
        end

        % Plot CO2 conversion
        figure(1);
        plot(w, co2Conversion, 'DisplayName', sprintf('%dK', t0));
        hold on;
    end

    % Final formatting for CO2 conversion plot
    figure(1);
    xlabel('Catalyst Mass [kg]');
    ylabel('CO2 Conversion');
    legend('show', 'Location', 'best');
    yline(targetConversion , '--r', 'DisplayName', 'Desired Conversion', 'LineWidth', 1.2);
    grid off;
    hold off;

    % Extract valid data for regression and cost analysis
    validIdx = ~isnan(wAtTargetConversion);
    tValid = temps(validIdx);
    wValid = wAtTargetConversion(validIdx);

    % Plot Temperature vs Catalyst Weight
    figure(2);
    xlabel('Temperature [K]');
    ylabel('Catalyst Weight for 0.376 Conversion [kg]');
    title('Catalyst Weight Required vs Temperature');
    grid off;
    hold on;
    % Polynomial regression
    p = polyfit(tValid, wValid, 2);
    weightFcn = @(T) polyval(p, T);
    Tplot = linspace(min(tValid), max(tValid), 100);
    plot(Tplot, weightFcn(Tplot), 'b-', 'DisplayName', 'Polynomial Fit');
    legend('show', 'Location', 'best');
    hold off;

    % Costing parameters
    Cc = 22.46;     % $/kg, catalyst cost
    Ch = 0.1026/(60*1000);    % $/W, heating cost per watt
    Tamb = 300;     % K, ambient temperature

    % Define cost functions
    catalystCostFcn = @(T) Cc * weightFcn(T);
    heatingPowerFcn = @(T) ...
        (schomate(params,T/1000, 'CO2') * params.eb.CO2.Fin + ...
         schomate(params,T/1000, 'H2')  * params.eb.H2.Fin  + ...
         schomate(params,T/1000, 'CO')  * params.eb.CO.Fin  + ...
         schomate(params,T/1000, 'CH4') * params.eb.CH4.Fin) * (T - Tamb); % Watts
    heatingCostFcn = @(T) Ch * heatingPowerFcn(T);
    totalCostFcn = @(T) catalystCostFcn(T) + heatingCostFcn(T);

    % Optimise temperature to minimise total cost
    Topt = fminbnd(totalCostFcn, min(tValid), max(tValid));

    % Plot cost vs temperature
    totalCost = arrayfun(totalCostFcn, Tplot);
    figure(3);
    plot(Tplot, totalCost, 'r-', 'LineWidth', 2, 'DisplayName', 'Total Cost');
    hold on;
    plot(Topt, totalCostFcn(Topt), 'ro', 'MarkerSize', 8, 'DisplayName', ...
    sprintf('Optimal T = %.1f K', Topt));
    xlabel('Temperature [K]');
    ylabel('Total Cost [$]');
    legend('show', 'Location', 'best');
    grid off;
    hold off;

    % Evaluate catalyst weight at the optimal temperature
    optimalCatalystWeight = weightFcn(Topt);

    % Display the result
    fprintf('Optimal Temperature: %.2f K\n', Topt);
    fprintf('Catalyst Weight at Optimal Temperature: %.2f kg\n', optimalCatalystWeight);
end


%%
function optimisePressure(init, params)
    % This function runs the ODE with 10 different pressures, 
    % 5 above and 5 below 22 bar, in 1 bar increments - (In pascals)
    
    % Define pressure range (5 below and 5 above 22 bar in 1 bar increments)
    pressure = 22 * 100000 + (-5:5) * 100000;

    % Store results
    W_all = cell(1, length(pressure)); % Store W values
    P_all = cell(1, length(pressure)); % Store Pressure profiles
    FA_all = cell(1, length(pressure)); % Store FA values
    CO2ConversionAll = cell(1, length(pressure)); % Store CO2 conversion

    % Prepare figures before loop
    figure(4); % For CO2 conversion profiles
    figure(5); % For pressure profiles
    figure(6); % For pressure drop profiles

    % Loop over each pressure
    for i = 1:length(pressure)
        P0 = pressure(i);
        
        % Assign initial conditions to vector
        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, 1100, P0];
        
        % Define span for integration
        Wspan = [0 5000];
        condition = 'LH';

        % Solve ODE
        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params, condition), Wspan, Y0);

        % Store data
        W_all{i} = W; 
        P_all{i} = Y(:,6); % Extract Pressure
        FA_all{i} = Y(:,1); % Extract FA
        
        % Compute CO2 conversion
        CO2ConversionAll{i} = (params.eb.CO2.Fin - FA_all{i}) / params.eb.CO2.Fin;

        % Plot W vs. CO2 Conversion
        figure(4);
        plot(W, CO2ConversionAll{i}, 'DisplayName', sprintf('%dbar', pressure(i)/100000));
        hold on;

        % Plot W vs. Pressure
        figure(5);
        plot(W, P_all{i}, 'DisplayName', sprintf('%dbar', pressure(i)/100000));
        hold on;

        % Calculate pressure drop
        pressureDrop = P0 - P_all{i};

        % Plot W vs. Pressure Drop
        figure(6);
        plot(W, pressureDrop, 'DisplayName', sprintf('%dbar', P0 / 100000));
        hold on;
    end 

    % Formatting for CO2 Conversion Profile
    figure(4);
    xlabel('Catalyst Mass (kg)');
    ylabel('CO2 Conversion');
    yline(0.376 , '--r', 'DisplayName', 'Desired Conversion', 'LineWidth', 1.2);
    legend show;
    grid off;
    hold off;

    % Formatting for Pressure Profile
    figure(5);
    xlabel('Catalyst Mass (kg)');
    ylabel('Pressure (Pa)');
    legend show;
    grid off;
    hold off;

    % Formatting for Pressure Drop Profile
    figure(6);
    xlabel('Catalyst Mass (kg)');
    ylabel('Pressure Drop (Pa)');
    legend show;
    grid off;
    hold off;
end

%%
function optimiseParticleDiamater(init, params)
    
    diameters = params.ergun.particleDiameter + (-5:5) * 0.0005;

    W_all = cell(1, length(diameters));      
    P_all = cell(1, length(diameters));      
    FA_all = cell(1, length(diameters));     
    CO2ConversionAll = cell(1, length(diameters));  

    P0 = 23 * 100000;  % Initial pressure in Pa
    targetConversion = 0.376;

    figure(6); % CO2 conversion profiles
    figure(7); % Pressure profiles
    figure(8); % Pressure drop profiles

    for i = 1:length(diameters)
        params.ergun.particleDiameter = diameters(i);

        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, 1100, P0];
        Wspan = [0 3000];
        condition = 'LH';

        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params, condition), Wspan, Y0);

        W_all{i} = W;
        P_all{i} = Y(:,6);
        FA_all{i} = Y(:,1);

        CO2ConversionAll{i} = (params.eb.CO2.Fin - FA_all{i}) / params.eb.CO2.Fin;

        % Plot CO2 Conversion
        figure(6);
        plot(W, CO2ConversionAll{i}, 'DisplayName', sprintf('%.4fm', diameters(i)));
        hold on;

        % Plot Pressure
        figure(7);
        plot(W, P_all{i}, 'DisplayName', sprintf('%.4fm', diameters(i)));
        hold on;

        % Plot Pressure Drop
        pressureDrop = P0 - P_all{i};
        figure(8);
        plot(W, pressureDrop, 'DisplayName', sprintf('%.4fm', diameters(i)));
        hold on;
    end

    figure(6);
    xlabel('Catalyst Mass (kg)');
    ylabel('CO2 Conversion');
    yline(0.376 , '--r', 'DisplayName', 'Desired Conversion', 'LineWidth', 1.2);
    legend show;
    grid off;
    hold off;

    figure(7);
    xlabel('Catalyst Mass (kg)');
    ylabel('Pressure (Pa)');
    legend show;
    grid off;
    hold off;

    figure(8);
    xlabel('Catalyst Mass (kg)');
    ylabel('Pressure Drop (Pa)');
    legend show;
    grid off;
    hold off;

    % Now find catalyst mass and pressure drop at target conversion
    objectives = zeros(1, length(diameters));

    for i = 1:length(diameters)
        conv = CO2ConversionAll{i};
        W = W_all{i};
        P = P_all{i};

        % Find first index where conversion >= targetConversion
        idx = find(conv >= targetConversion, 1, 'first');

        if isempty(idx)
            % Conversion not reached within Wspan, assign large value or NaN
            objectives(i) = NaN;
            continue
        end

        W_target = W(idx);
        deltaP = P0 - P(idx);

        if deltaP == 0
            deltaP = 1e-6; % avoid div zero
        end

        % Objective: catalyst mass / pressure drop (lower is better)
        objectives(i) = W_target / deltaP;
    end

    % Plot objective vs diameter (ignore NaNs)
    figure;
    validIdx = ~isnan(objectives);
    plot(diameters(validIdx), objectives(validIdx), '-o', 'LineWidth', 1.5);
    xlabel('Particle Diameter (m)');
    ylabel('Catalyst Mass / Pressure Drop');
    grid off;
    hold on;

    [minVal, minIdx_rel] = min(objectives(validIdx));
    minIdx = find(validIdx);
    minIdx = minIdx(minIdx_rel);

    plot(diameters(minIdx), minVal, 'ro', 'MarkerSize', 12);
    text(diameters(minIdx), minVal, sprintf('  Min at %.4fm', diameters(minIdx)));
    hold off;
end


function optimiseSuperficialVelocity(init, params)
    % Function to evaluate reactor length needed for a target conversion
    % and plot conversion, pressure, pressure drop, superficial velocity,
    % and length vs diameter

    params.ergun.particleDiameter = 0.0035;
    targetConversion = 0.376;
    diameters = 1 + (-5:5) * 0.1; % [0.5 to 1.5 m]
    lengths = zeros(1, length(diameters));
    supVels = zeros(1, length(diameters)); % Superficial velocity array

    % Prepare figures without grid lines and titles
    figure(8); clf; hold on;
    xlabel('Catalyst Mass (kg)');
    ylabel('CO2 Conversion');

    figure(9); clf; hold on;
    xlabel('Catalyst Mass (kg)');
    ylabel('Pressure (Pa)');

    figure(12); clf; hold on;
    xlabel('Catalyst Mass (kg)');
    ylabel('Pressure Drop (Pa)');

    for i = 1:length(diameters)
        D = diameters(i);
        params.reactor.diameter = D;

        % Update Ergun parameters
        params.ergun.csArea = (pi / 4) * D^2;
        params.ergun.supVel = params.inlet.totalVolFlowrate / params.ergun.csArea;
        params.ergun.gasFlux = params.ergun.inletDensity * params.ergun.supVel;
        supVels(i) = params.ergun.supVel;

        % Initial conditions
        Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, 1100, 23e5];
        Wspan = [0 3000];
        condition = 'LH';

        % Solve ODE
        [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params, condition), Wspan, Y0);

        FA = Y(:,1);
        P = Y(:,6);
        conversion = (params.eb.CO2.Fin - FA) / params.eb.CO2.Fin;

        % Plot conversion vs. W
        figure(8);
        plot(W, conversion, 'DisplayName', sprintf('%.2f m', D));

        % Plot pressure vs. W
        figure(9);
        plot(W, P, 'DisplayName', sprintf('%.2f m', D));

        % Calculate pressure drop and plot
        deltaP = P(1) - P;
        figure(12);
        plot(W, deltaP, 'DisplayName', sprintf('%.2f m', D));

        % Find W for target conversion and compute length
        idx = find(conversion >= targetConversion, 1);
        if ~isempty(idx)
            requiredW = W(idx);
            bulkDensity = params.ergun.bulkDensity;
            lengths(i) = (4 * requiredW) / (bulkDensity * pi * D^3);
        else
            warning(['Target conversion not reached for diameter = ', num2str(D)]);
            lengths(i) = NaN;
        end
    end

    figure(8); legend show;
    figure(9); legend show;
    figure(12); legend show;

%     % Reactor Length vs Diameter
%     figure(10); clf;
%     plot(diameters, lengths, '-', 'LineWidth', 1.5);
%     xlabel('Reactor Diameter (m)');
%     ylabel('Reactor Length (m)');

    % Superficial Velocity vs Diameter
    figure(11); clf;
    plot(diameters, supVels, '-', 'LineWidth', 1.5);
    xlabel('Reactor Diameter (m)');
    ylabel('Superficial Velocity (m/s)');
end


function optimiseParticleAndBedDiameter(init, params)
    % Define ranges
    d_particles = 0.0035 : 0.0005 : 0.0085;% particle diameters (m)
    d_beds = 0.5:0.1:1.5;              % bed diameters (m)

    % Constants
    targetConversion = 0.376;
    P0 = 23e5; % Initial pressure

    % Preallocate objective matrix
    objVals = NaN(length(d_particles), length(d_beds));

    % Result storage for analysis
    bestConfig = struct('dp', NaN, 'D', NaN, 'W', NaN, 'deltaP', NaN, 'L', NaN, 'i', NaN, 'j', NaN);
    minObj = Inf;

    % Loop over all combinations
    for i = 1:length(d_particles)
        for j = 1:length(d_beds)
            dp = d_particles(i);
            D = d_beds(j);

            % Update parameters
            params.ergun.particleDiameter = dp;
            params.reactor.diameter = D;
            params.ergun.csArea = (pi/4) * D^2;
            params.ergun.supVel = params.inlet.totalVolFlowrate / params.ergun.csArea;
            params.ergun.gasFlux = params.ergun.inletDensity * params.ergun.supVel;

            % Initial conditions and solve ODE
            Y0 = [init.FA0, init.FB0, init.FC0, init.FD0, 1100, P0];
            Wspan = [0 3000];
            condition = 'LH';

            try
                [W, Y] = ode45(@(w, Y) odeSolver(w, Y, params, condition), Wspan, Y0);
            catch
                warning('ODE solver failed at dp=%.4f, D=%.2f', dp, D);
                continue;
            end

            FA = Y(:,1);
            P = Y(:,6);
            conversion = (params.eb.CO2.Fin - FA) / params.eb.CO2.Fin;

            idx = find(conversion >= targetConversion, 1);
            if isempty(idx)
                continue;
            end

            Wtarget = W(idx);
            deltaP = P0 - P(idx);

            % Reactor length
            bulkDensity = params.ergun.bulkDensity;
            L = (4 * Wtarget) / (bulkDensity * pi * D^2);

            % Apply geometry constraint: D must be <= 0.33 * L
            if D > 0.33 * L
                continue;
            end

            % Objective function: catalyst mass × pressure drop
            obj = Wtarget / deltaP;
            objVals(i, j) = obj;

            % Store best result
            if obj < minObj
                minObj = obj;
                bestConfig.dp = dp;
                bestConfig.D = D;
                bestConfig.W = Wtarget;
                bestConfig.deltaP = deltaP;
                bestConfig.L = L;
                bestConfig.i = i;
                bestConfig.j = j;
            end
        end
    end

    % Plot the objective function as a surface
    [D_grid, dp_grid] = meshgrid(d_beds, d_particles);
    figure;
    surf(D_grid, dp_grid, objVals);
    xlabel('Bed Diameter (m)');
    ylabel('Particle Diameter (m)');
    zlabel('Objective: W \times \DeltaP');
    grid off;
    colorbar;
    hold on;

    % Add red marker at optimum point
    plot3(bestConfig.D, bestConfig.dp, objVals(bestConfig.i, bestConfig.j), ...
        'ro', 'MarkerSize', 10, 'LineWidth', 2);
    legend('Objective Surface', 'Optimal Configuration', 'Location', 'northeast');

    % Output best configuration
    fprintf('Best config:\n');
    fprintf('Particle Diameter = %.4f m\n', bestConfig.dp);
    fprintf('Bed Diameter      = %.2f m\n', bestConfig.D);
    fprintf('Catalyst Mass     = %.2f kg\n', bestConfig.W);
    fprintf('Pressure Drop     = %.2f Pa\n', bestConfig.deltaP);
    fprintf('Reactor Length    = %.2f m\n', bestConfig.L);
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
    case 'H2O'
        A = params.eb.H2O.A;
        B = params.eb.H2O.B;
        C = params.eb.H2O.C;
        D = params.eb.H2O.D;
        E = params.eb.H2O.E;
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

function conversionCO2s = conversion(params, FAs)
    % Calculate conversion values
    
    conversionCO2s = zeros(length(FAs), 1);
    for i = 1:length(FAs)
        conversionCO2s(i) = (params.eb.CO2.Fin - FAs(i)) / (params.eb.CO2.Fin);
    end

end

function reactorLengths = reactorLengthFunc(params,ws)

    % Calculate reactor length

    reactorLengths = zeros(length(ws), 1);
    for i = 1:length(ws)
        reactorLengths(i) = (4*ws(i)/(params.ergun.bulkDensity*pi*(params.reactor.diameter^3)));
    end 

end