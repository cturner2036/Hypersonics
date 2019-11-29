function nu = PME_nu(gamma,Mach)
% Simple Prandtl-Meyers Function
%   Input gamma value and Mach value to get prandtl-meyer value
lam = sqrt((gamma-1)/(gamma+1));
B = sqrt((Mach^2)-1);
nu = (1/lam)*atand(lam*B) - atand(B);
end
