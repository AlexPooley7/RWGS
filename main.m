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

% Arrhenius constants

% Rate equation constants


Wspan = [0 1000];
[w,Y] = ode15s(@questiontwo_func, Wspan, Y0, options); % Calls the ode
T = Y(:,1); CA = Y(:,2); CB = Y(:,3); CS = Y(:,4); CD = Y(:,5); P = Y(:,6);

    %%
function dYdt = odeSolver(w,Y)

% Translate components of input vector into variables from the model
FA = Y(2); FB = Y(3); FC = Y(4); FD = Y(5); T = Y(1); P = Y(6);

% Rate constant calculations using the Arrhenius equation
k = A*exp(-Ea_RWGS./(R*T));

% Rate of reaction calculations
rRWGS = ; 

%%   
% Mole balances ODEs
dFA_dw = -rRWGS;
dFB_dw = -rRWGS;
dFC_dw = rRWGS;
dFD_dw = rRWGS;

% Energy balance ode
dT_dw;

% Pressure ODE
dP_dw;

% Forms output column vector for ode solver
dYdt = [dFA_dw; dFB_dw; dFC_dw; dFD_dw; dT_dw; dP_dw];


end