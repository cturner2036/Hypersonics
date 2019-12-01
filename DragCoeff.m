function [D] = DragCoeff(Mach,alpha,phi,NFC,AFC,q,weight)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variables:
% - Mach: Mach number of the vehicle
% - alpha: Angle of attack of the vehicle
% - phi: Equivalence Ratio of RJE (May need to set permaently to account
% for difference in engine between GHV
% - NFC: Normal Force Coefficient Table from GHV data {run Cd_Import.m} [N]
% - AFC: Axial Force Coefficient Table from GHV data {run Cd_Import.m to
% obtain} [N]
% - q: dynamic pressure atmosphere your flying through 
%%%
%Outputs:
% - D: Drag Value [N]
% - L: Lift Value [N]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the data before placing function in loop to save on speed

%Normal Force Coeff.
%if Phi
if phi == .5
    TruncNFC = NFC(1:6,:);
elseif phi == .6
    TruncNFC = NFC(7:12,:);
elseif phi == .7
    TruncNFC = NFC(13:18,:);
elseif phi == .8
    TruncNFC = NFC(19:24,:);
elseif phi == .9
    TruncNFC = NFC(25:30,:);
elseif phi == 1
    TruncNFC = NFC(31:36,:);
elseif phi == 1.1
    TruncNFC = NFC(37:42,:);
elseif phi == 1.2
    TruncNFC = NFC(43:48,:);
end
%Angle of Attack
if alpha == -4
    TruncNFC = TruncNFC(1,:);
elseif alpha == -2
    TruncNFC = TruncNFC(2,:);
elseif alpha == 0
    TruncNFC = TruncNFC(3,:);
elseif alpha == 2
    TruncNFC = TruncNFC(4,:);
elseif alpha == 4
    TruncNFC = TruncNFC(5,:);
elseif alpha == 5
    TruncNFC = TruncNFC(6,:);
end 
 %Mach
if Mach == 4
    NCd = TruncNFC.M4;
elseif Mach > 4 && Mach < 4.5
    NCd = TruncNFC.M4 + (((Mach - 4)/(4.5-4))*(TruncNFC.M4_5-TruncNFC.M4)); 
elseif Mach == 4.5
    NCd = TruncNFC.M4_5;
elseif Mach > 4.5 && Mach < 5
    NCd = TruncNFC.M4_5 + (((Mach - 4.5)/(5-4.5))*(TruncNFC.M5-TruncNFC.M4_5)); 
elseif Mach == 5
    NCd = TruncNFC.M5;
elseif Mach > 5 && Mach < 5.5
    NCd = TruncNFC.M5 + (((Mach - 5)/(5.5-5))*(TruncNFC.M5_5-TruncNFC.M5)); 
elseif Mach == 5.5
    NCd = TruncNFC.M5_5;
elseif Mach > 5.5 && Mach < 6
    NCd = TruncNFC.M5_5 + (((Mach - 5.5)/(6-5.5))*(TruncNFC.M6-TruncNFC.M5_5)); 
elseif Mach == 6
    NCd = TruncNFC.M6;
elseif Mach > 6 && Mach < 6.5
    NCd = TruncNFC.M6 + (((Mach - 6)/(6.5-6))*(TruncNFC.M6-TruncNFC.M6_5)); 
elseif Mach == 6.5
    NCd = TruncNFC.M6_5;
elseif Mach > 6.5 && Mach < 7
    NCd = TruncNFC.M6_5 + (((Mach - 6.5)/(7-6.5))*(TruncNFC.M7-TruncNFC.M6_5)); 
elseif Mach == 7
    NCd = TruncNFC.M7;
end

%Axial Force Coeff.
%if Phi
if phi == .5
   TruncAFC = AFC(1:6,:);
elseif phi == .6
    TruncAFC = AFC(7:12,:);
elseif phi == .7
    TruncAFC = AFC(13:18,:);
elseif phi == .8
    TruncAFC = AFC(19:24,:);
elseif phi == .9
    TruncAFC = AFC(25:30,:);
elseif phi == 1
    TruncAFC = AFC(31:36,:);
elseif phi == 1.1
    TruncAFC = AFC(37:42,:);
elseif phi == 1.2
    TruncAFC = AFC(43:48,:);
end
%Angle of Attack
if alpha == -4
    TruncAFC = TruncAFC(1,:);
elseif alpha == -2
    TruncAFC = TruncAFC(2,:);
elseif alpha == 0
    TruncAFC = TruncAFC(3,:);
elseif alpha == 2
    TruncAFC = TruncAFC(4,:);
elseif alpha == 4
    TruncAFC = TruncAFC(5,:);
elseif alpha == 5
    TruncAFC = TruncAFC(6,:);
end 
 %Mach
 if Mach == 4
     ACd = TruncAFC.M4;
 elseif Mach > 4 && Mach < 4.5
     ACd = TruncAFC.M4 + (((Mach - 4)/(4.5-4))*(TruncAFC.M4_5-TruncAFC.M4)); 
 elseif Mach == 4.5
     ACd = TruncAFC.M4_5;
 elseif Mach > 4.5 && Mach < 5
     ACd = TruncAFC.M4_5 + (((Mach - 4.5)/(5-4.5))*(TruncAFC.M5-TruncAFC.M4_5)); 
 elseif Mach == 5
     ACd = TruncAFC.M5;
 elseif Mach > 5 && Mach < 5.5
     ACd = TruncAFC.M5 + (((Mach - 5)/(5.5-5))*(TruncAFC.M5_5-TruncAFC.M5)); 
 elseif Mach == 5.5
     ACd = TruncAFC.M5_5;
 elseif Mach > 5.5 && Mach < 6
     ACd = TruncAFC.M5_5 + (((Mach - 5.5)/(6-5.5))*(TruncAFC.M6-TruncAFC.M5_5)); 
 elseif Mach == 6
     ACd = TruncAFC.M6;
 elseif Mach > 6 && Mach < 6.5
     ACd = TruncAFC.M6 + (((Mach - 6)/(6.5-6))*(TruncAFC.M6-TruncAFC.M6_5)); 
 elseif Mach == 6.5
     ACd = TruncAFC.M6_5;
 elseif Mach > 6.5 && Mach < 7
     ACd = TruncAFC.M6_5 + (((Mach - 6.5)/(7-6.5))*(TruncAFC.M7-TruncAFC.M6_5)); 
 elseif Mach == 7
     ACd = TruncAFC.M7;
 end

 

%Calculate Drag
Cd = ACd;
Cl = NCd;
Lift_to_Drag = Cl/Cd;
D = (weight/Lift_to_Drag)/1000; %[kN]

end