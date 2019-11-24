function MFP=FlowMFP(M,gamma,R,g_c)
%%FLOWMFP Mass Flow Parameter Function By Joshua McMillin
% Function returns the Mass Flow Parameter for given Flow Parameters
% M=Mach#
% gamma= ratio of specific heats
%R=gas constant
%g_c=gravitational constant
MFP = sqrt(gamma*g_c/R).*M./(1+0.5.*(gamma-1).*M.^2).^((gamma+1)/(2*(gamma-1)));
end
