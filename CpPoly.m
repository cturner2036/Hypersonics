function Cp = CpPoly(T)
    T = T.*1.8; %Convert to Rankine
    Cp = (3.2240709 + 2.5453604e-4.*T + 2.5222882e-7.*T.^2 -1.2013920e-10.*T.^3 + 1.7003460e-14.*T.^4).*1716.5;
    Cp = Cp.*0.1672;%J/kg/K Conversion
end