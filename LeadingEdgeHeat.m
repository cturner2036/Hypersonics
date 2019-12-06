function [qt, qw] = LeadingEdgeHeat(M, a,Tinf,Pinf,Rhoinf, Tw)
% This function computes the area specific Leading Edge heat.
% Parameters computed at total conditions
% Inputs: 
%   M: Mach number in duct []
%   a: angle of attack [deg]
%   Tinf: Static Temp [K]
%   Pinf: Static Pressure [Pa]
%   Rhoinf: Static Density
%   Tw: allowable wall temperature [K]

rn = 1/100; %radius of leading edge lip (assume 1 cm)

Gamma = 1.4;  
[~, Pr, Rho, ~] = normCond(M, Gamma);
P_e = Pinf*Pr;
Rho_e = Rhoinf*Rho;

rho_w = 0.3627; %kg/m/s
mu_w = 4.111e-5; %kg/m/s

Tt = Tinf*(1 + (Gamma-1)/2*M^2);
mu_t = 6.63e-5; %kg/m/s
rho_t = 0.1553; %kg/m^3 
Cp_t = 1264; %J/kg/K
dVdX = sqrt(2*(P_e-Pinf)/Rho_e)/rn; %

Kl = cosd(a); % Multiplication Factor?



Ht = Cp_t*Tt; Hw = Cp_t*Tw; %specific entalpy at recovery and wall temps

qt = 0.763*9.81*Kl*(rho_t*mu_t)^0.4*(rho_w*mu_w)^0.1*dVdX^0.5*Ht;
qw = 0.763*9.81*Kl*(rho_t*mu_t)^0.4*(rho_w*mu_w)^0.1*dVdX^0.5*Hw;
end

