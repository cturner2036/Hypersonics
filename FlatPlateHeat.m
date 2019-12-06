function [qr, qw] = FlatPlateHeat(M,a,T,x,Tw)
%This function computes the area specific heat transfer for a flat plate
% Parameters computed at recovery temp.  Assumed to be 2,000 C (roughly
% midpoint of Mach 4 to Mach 8 revocery temp.
% Inputs: 
%   M: Mach number in duct []
%   a: angle of attack [deg]
%   T: Static Temp [K]
%   x: distance from flow leading edge [m]
%   Tw: allowable wall temperature [K]

Gamma = 1.4; R = 287; 

mu_r = 6.63e-5; %kg/m/s
rho_r = 0.1553; %kg/m^3
Cp_r = 1264; %J/kg/K
k_r = 0.11113; %W/m/K
Ve = M*sqrt(Gamma*R*T); %duct flow velocity defined by Mach number [m/s]

Kt = cosd(a); % Multiplication Factor?
Re_r = rho_r*Ve*x/mu_r; %Reynolds number at ref conditions
Pr_r = mu_r*Cp_r/k_r; %Prandtl's number at ref cond
r = Pr_r^(1/3); %recovery factor
Tr = T*(1 + (Gamma-1)/2*M^2*r); %recovery temperature

Hr = Cp_r*Tr; Hw = Cp_r*Tw; %specific entalpy at recovery and wall temps

qr = 0.185*9.81*Kt/Pr_r*mu_r*rho_r*Ve*Hr*2.584/log10(Re_r + 3000);
qw = 0.185*9.81*Kt/Pr_r*mu_r*rho_r*Ve*Hw*2.584/log10(Re_r + 3000);
end

