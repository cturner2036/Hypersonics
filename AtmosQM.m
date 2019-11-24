function [T, P, rho,hgeom] = AtmosQM(Q,M)
% Geometric altitude vs. Mach
% Q and M must be same size
% Q and P in /bf/ft^2
% T is in deg R
% rho is slug/ft^3
hgeom=cumsum(ones(2992,1)*100);
[T_atmos,P_atmos,rho_atmos,~]=atmosphere4(hgeom,1);
gamma=1.4;
P = (2.*Q)./(gamma.*M.^2);
h=interp1(P_atmos,hgeom,P,'pchip');
[T,~,rho,hgeom]=atmosphere4(h,1);
end 
