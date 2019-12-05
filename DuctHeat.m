function [qr, qw] = DuctHeat(M, T, x, Tw)
%This function computes the area specific heat transfer for duct flow
% reference temperature is conventionally 25 degrees Celsius 
% Inputs: 
%   M: Mach number in duct []
%   T: Static Temp in duct [K]
%   x: distance from flow leading edge [m]
%   Tw: allowable wall temperature [K]

Gamma = 1.4; R = 287; 

mu_ref = 1.849e-5; %kg/m/s
rho_ref = 1.184; %kg/m^3
Cp_ref = 1007; %J/kg/K
k_ref = 0.02551; %W/m/K
Ve = M*sqrt(Gamma*R*T); %duct flow velocity defined by Mach number [m/s]


Re_ref = rho_ref*Ve*x/mu_ref; %Reynolds number at ref conditions
Pr_ref = mu_ref*Cp_ref/k_ref; %Prandtl's number at ref cond
r = Pr_ref^(1/3); %recovery factor
Tr = T*(1 + (Gamma-1)/2*M^2*r); %recovery temperature

Hr = Cp_ref*Tr; Hw = Cp_ref*Tw; %specific entalpy at recovery and wall temps

qr = 0.0296*Re_ref^-0.2*Pr_ref^-0.67*rho_ref*Ve*Hr; %heat [W/m^2]
qw = 0.0296*Re_ref^-0.2*Pr_ref^-0.67*rho_ref*Ve*Hw;
end

