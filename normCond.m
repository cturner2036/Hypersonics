function [Ptr, Pr, Rho, M1] = normCond(M, g)
% compute condition ratios across a normal shock for given Mach and Gamma
% -- ratios are conditions after shock/conditions before
Ptr = (((g+1)*M^2)/((g-1)*M^2+2))^(g/(g-1))*((g+1)/(2*g*M^2-(g-1)))^(1/(g-1));
M1 = sqrt(((g-1)*M^2+2)/(2*g*M^2-(g-1)));
Pr = (2*g*M^2 - (g-1))/(g+1);
Rho = ((g+1)*M^2)/((g-1)*M^2 + 2);
end

