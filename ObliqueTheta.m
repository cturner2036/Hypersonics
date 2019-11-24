%% Theta Function for known Beta by Joshua E. McMillin
function t=ObliqueTheta(B1, M, gamma)
    t = B1 - atan((2+(gamma-1)*(M*sin(B1))^2)/((gamma+1)*sin(B1)*cos(B1)*M^2));
    % Elements of Gas Turbine Propulsion, J.D. Mattingly, 1996
end
