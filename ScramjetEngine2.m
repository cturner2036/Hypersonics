function EnginePerf = ScramjetEngine(InletMap,DynamicPressure,FreestreamMach,AngleofAttack,Dry_Weight,Fuel_Weight,SI_Flag,AFC,NFC,FP_Angle)
    %Constants
    Total_Weight = Dry_Weight + Fuel_Weight;
    gamma = 1.4;
    %R_BTU_lbmolR = 1.986;
    %R_ftlbf_lbmolR = 1545;
    R_J_kmolK = 8314; 
    MW_air = 28.965;
    %Cp_air_BTUlbmR = (R_BTU_lbmolR/MW_air)*gamma/(gamma-1);
    CP_air_j_kmolK = (R_J_kmolK/MW_air)*gamma/(gamma-1);
    grav = 9.81; %[m/s^2]
    drop_height = 12192; %[m]
    [T_amb, P_amb, rho] = atmosphere4(drop_height,0);
    xx(1) = 0;
    yy(1) = 12192;
    MM(1) = FreestreamMach;
    Isp(1) = 3000;
    timestep = 1; %[s]
    cnt = 1;
    i = 2;
    while Total_Weight >= Dry_Weight
        %Assume run 1 second iterations, depending on mf_dot will determine
        
        % get inlet performance data
        Inlet = getGHV(InletMap,FreestreamMach,AngleofAttack,DynamicPressure,SI_Flag);

        % Station 0 - Freestream
        Station0=Station;
        Station0.Mach = FreestreamMach;
        Station0.Velocity_ms = Inlet.FreestreamVelocity_ms;
        Station0.Pressure_Pa = P_amb;
        Station0.TotalPressure_Pa = Inlet.FreestreamTotalPressure_Pa;
        Station0.Temperature_K = Inlet.FreestreamTemperature_K;
        Station0.TotalTemperature_K = Inlet.FreeStreamTotalTemperature_K;
        Station0.MassFlowRate_kgs = Inlet.MassFlowRate_kgs;
        Station0.Area_m2 = Inlet.EffectiveStreamtubeCapture_m2;
        EnginePerf.Station0 = Station0;

        % Station 1 - Isolator Entrance
        Station1=Station;
        Station1.Mach = Inlet.ThroatMachNumber;
        Station1.Velocity_ms = Inlet.ThroatVelocity_ms;
        Station1.Pressure_Pa = Inlet.ThroatPressure_Pa;
        Station1.TotalPressure_Pa = Inlet.ThroatTotalPressure_Pa;
        Station1.Temperature_K = Inlet.ThroatTemperature_K;
        Station1.TotalTemperature_K = Inlet.ThroatTotalTemperature_K;
        Station1.MassFlowRate_kgs = Inlet.MassFlowRate_kgs;
        Station1.Area_m2 = Inlet.ThroatArea_m2;
        EnginePerf.Station1 = Station1;

        %Station 2 - Isolator Exit / Combustor Entrance
        Station2 = Station;
        Station2.Mach = Inlet.IsolatorExitMach;
        Station2.Velocity_ms = Inlet.IsolatorExitMach*sqrt(gamma*(R_J_kmolK/MW_air)*Inlet.IsolatorExitTemperature_K);
        Station2.Pressure_Pa = Inlet.IsolatorExitPressure_Pa;
        Station2.TotalPressure_Pa = Inlet.IsolatorExitTotalPressure_Pa;
        Station2.Temperature_K = Inlet.IsolatorExitTemperature_K;
        Station2.TotalTemperature_K = Inlet.IsolatorExitTotalTemperature_K;
        Station2.MassFlowRate_kgs = Inlet.MassFlowRate_kgs;
        Station2.Area_m2 = Inlet.IsolatorEffectiveExitArea_m2;
        EnginePerf.Station2 = Station2;

        %get combustor data
        Station4 = getCombustorOutlet(Station2);
        EnginePerf.Station4 = Station4;

        % get Nozzle data
        Station9 = getNozzle(Station4,Station0);
        EnginePerf.Station9 = Station9;
        %% Adding Scripting for Nozzle%%

        %Calculate Total Thrust from Engine
        %   Set Nozzle Geom.
        if cnt == 1
            %   Ideal Nozzle
            %[throat_height, throat_angle, cowl_height, body_width, step_size] = Plug_Nozzle(throat_angle, throat_height, cowl_height, body_width, step_size);
            [throat_height,throat_angle,cowl_height,body_width,step_size] = Plug_Nozzle(67,0.2,1,1,100);
            %   Truncated Nozzle
            %[x,y] = Plug_Nozzle_Style2(AR,eta_b,throat_height,step_size);
            [x,y,local_turn] = Plug_Nozzle_Style2(3.3839,0.05,throat_height,step_size);
            spd_snd = sqrt(gamma*R_J_kmolK*Station0.Temperature_K);
            veloc = FreestreamMach*spd_snd;
            height = drop_height;
        end
        %   Run Flow_Properties to calculate Thrust Values
        %[Engine_Thrust, Engine_Lift] = Flow_Properties(step_size,local_turn,P_amb,T_exit,Pt_exit(Station4.TotalPressure_Pa),throat_angle,throat_height,body_width,M_throat,y,x,alpha,Q,mdot,gamma);
        [Engine_Thrust, Engine_Lift, StagnationTemp] = Flow_Properties(step_size,local_turn,Station0.Pressure_Pa,Station4.Temperature_K,Station4.TotalPressure_Pa,throat_angle,throat_height,body_width,1.15,y,x,AngleofAttack,71820,Station4.MassFlowRate_kgs,1.4);
  
        % Calculate Aero Drag Thrust
        %Run Cd_Import Before execution
        %D = DragCoeff(FreestreamMach, AngleofAttack, phi, NFC, AFC, DynamicPressure, TotalWeight);
        [D, L] = DragCoeff(FreestreamMach, AngleofAttack, 1.0, NFC, AFC, DynamicPressure, Total_Weight);
        Drag(i) = D;
        
        %Update Vehicle Weight and Update Flight Values
        %[total_weight] = actualWeight(full_weight,mdot_f(***need from combustor***),timestep[s])
        mdot_ff(i) = (Station4.MassFlowRate_kgs - Station2.MassFlowRate_kgs);
        [Total_Weight] = actualWeight(Total_Weight,mdot_ff(i),timestep);
        Vehicle_Weight(i) = Total_Weight;
        
        %Rocket EQ
        %F = ma > a = F/m   %%% Change AngleofAttack to Pitch
        Engine_Total_Thrust = (Engine_Thrust*1000 + Inlet.InletAxialForce_N - Station2.MassFlowRate_kgs*Station2.Velocity_ms);
        %Isp should be around 1000-2000
        Isp(i) = Engine_Total_Thrust/(mdot_ff(i)*grav);
        Total_Thrust = (Engine_Total_Thrust/1000) - D - (Total_Weight*grav*sind(FP_Angle))/1000;
        accel = (Total_Thrust*1000)/Total_Weight;
        accel_y = (L - (Total_Weight*grav*cosd(FP_Angle)))/Total_Weight;
        veloc = veloc + accel*timestep;
        distance = veloc*timestep;
        xx(i) = distance*cosd(FP_Angle) + xx(cnt);
        yy(i) = distance*sind(FP_Angle) + yy(cnt);
        MM(i) = FreestreamMach;
        
        
        %calculate new ambient conditions with new height
        height = height + (yy(i)-yy(cnt));
        height_ft = height*3.281;
        [T_amb, P_amb, rho] = atmosphere4(height_ft,0);
        P_amb = P_amb*47.88;
        P_amb_plot(i) = P_amb;
        St4_P(i) = Station4.TotalPressure_Pa;
        St2_P(i) = Station2.TotalPressure_Pa;
        accel_p(i) = accel;
        accely_p(i) = accel_y;
        spd_snd = sqrt(gamma*R_J_kmolK*Station0.Temperature_K);
        FreestreamMach = veloc/spd_snd;
        cnt = cnt + 1;
        i = i + 1
        
        %%%%%% Plotting Stuff %%%%%%%%%%%%%%%
        dd = linspace(0,i,i-1);
        tiledlayout(3,4);
        %h = figure;
        
        %Top Left Plot
        ax1y1 = nexttile;
        plot(ax1y1,xx,yy)
        title(ax1y1,["Distance (m), Flight Angle:" + num2str(FP_Angle)])
        xlabel('Distance (m)')
        ylabel("height (m)")
        %xlim([0 6*10^6])
        %ylim([12000 50000])
        grid on
        
        %Top Right Plot
        ax2y1 = nexttile(2);
        plot(ax2y1,dd,MM)
        title("Mach # vs time")
        ylabel('Mach #')
        xlabel('iterations')
        %xlim([0 1200])
        %ylim([4 8])
        
        %Middle Left Plot
        ax1y2 = nexttile(3);
        plot(ax1y2,dd,Isp)
        title("Isp vs. Iterations")
        xlabel("iterations")
        ylabel("Specific Thrust (Isp)")
        %xlim([0 1200])
        %ylim([2500 8100])
        
        %Middle Right Plot
        ax2y2 = nexttile(4);
        plot(ax2y2,dd,mdot_ff)
        title("fuel mass flow vs. iterations")
        xlabel("iterations")
        ylabel("fuel mass flow")
        %xlim([0 1200])
        %ylim([0 2])
        
        %Bottom Left Plot
        ax1y3 = nexttile(5);
        plot(ax1y3,dd,Vehicle_Weight)
        title("Vehicle Weight vs. Iterations")
        xlabel("iterations")
        ylabel("weight [kgs]")
        %xlim([0 1200])
        %ylim([2300 3900])
        
        %Bottom Right Plot
        ax2y3 = nexttile(6);
        plot(ax2y3,MM,Drag)
        title("Drag Coeff vs. M#")
        xlabel("Mach #")
        ylabel("Drag Coeff")
        
        %Bottom Right Plot
        ax1y4 = nexttile(7);
        plot(ax1y4,dd,P_amb_plot)
        title("Ambient Pressure over time")
        xlabel("iter")
        ylabel("Amb_Pressure (Pa)")
                
        %Bottom Right Plot
        ax2y4 = nexttile(8);
        plot(ax2y4,dd,St4_P)
        title("Station4 Pressure over time")
        xlabel("iter")
        ylabel("Station 4 (Pa)")
        
        %Bottom Right Plot
        ax2y5 = nexttile(9);
        plot(ax2y5,dd,St2_P)
        title("Station2 Pressure over time")
        xlabel("iter")
        ylabel("Station 4 (Pa)")
        
        %Bottom Right Plot
        ax2y6 = nexttile(10);
        plot(ax2y6,dd,accel_p)
        title("Acceleration Curve")
        xlabel("iter")
        ylabel("Acceleration (m/s^2)")
        
        %Bottom Right Plot
        ax1y3 = nexttile(11);
        plot(ax1y3,dd,accel_y)
        title("Acceleration Curve Y - dir")
        xlabel("iter")
        ylabel("Acceleration (m/s^2)")
        
        %saveas(h,sprintf('Fig%d.png',i));
    end

   %plot(xx,yy)
   %plot(i,MM)
   %plot(i,Isp)
   %plot mf_dot over i
    
end


function S4 = getCombustorOutlet(S2)
    gamma = 1.4;
    R_J_kmolK = 8314; 
    MW_air = 28.965;
    CP_air_J_kgK = (R_J_kmolK/MW_air)*gamma/(gamma-1);
    H_f_J_kg = 43.5e6; %heating value of Fuel
    T_f_K = 298; %initial Temperature of Fuel
    H_vap_J_kg = 446e3; % Heat of vaporization
    Cp_f_J_kgK = 480*4.1868; % Specific heat from cal/kgK to J/kgK

    S4 = Station;
    S4.Mach = 1; % adaptive geometry and fuel flow ensures flow is choked
    Tt4 = 2600;
    S4.TotalTemperature_K = Tt4; % adaptive fuel flow will maximize stagnation temperature
    S4.Temperature_K = S4.TotalTemperature_K / (1+0.5*(gamma-1)*S4.Mach^2); % solve for exit Temperature
    S4.Velocity_ms = S4.Mach *sqrt(gamma*(R_J_kmolK/MW_air)*S4.Temperature_K); % Solve for exit velocity
    
    % Solve for Fuel Flow Rate to hit Tt4
    Tt2 = S2.TotalTemperature_K;
    mdot_air = S2.MassFlowRate_kgs;
    dH_air = mdot_air*CP_air_J_kgK*(Tt4-Tt2); % Joules/sec
    mdot_f = dH_air /(H_f_J_kg - Cp_f_J_kgK*(Tt4-T_f_K)- H_vap_J_kg); %kg/sec
    S4.MassFlowRate_kgs = mdot_air + mdot_f;
    
    %Solve for Pressure at exit - This is an estimate... Cal's model is
    %real model
    S4.Pressure_Pa = S2.Pressure_Pa*((1+0.5*(gamma-1)*S2.Mach^2)/(1+0.5*(gamma-1))); % Estimate of pressure gain for supersonic combustion
    S4.TotalPressure_Pa = S4.Pressure_Pa*(1+0.5*(gamma-1)*S4.Mach^2)^(gamma/(gamma-1)); % Results in stagnation Pressure loss
    S4.TotalPressure_Pa = S4.TotalPressure_Pa*(0.4);
    %Solve for variable throat Area
    S4.Area_m2 = S4.MassFlowRate_kgs * (R_J_kmolK/MW_air) * S4.Temperature_K / (S4.Pressure_Pa*S4.Velocity_ms); % m^2
end

function S9 = getNozzle(S4,S0)
%Constants
    gamma = 1.4;
    R_J_kmolK = 8314; 
    MW_air = 28.965;
% conservation of mass
    S9.MassFlowRate_kgs = S4.MassFlowRate_kgs;
    
% For now, assume perfect expansion - update with Marks Nozzle model later
    Pt9 = S4.TotalPressure_Pa;
    P9 = S0.Pressure_Pa;
    S9.Pressure_Pa = P9;
    S9.TotalPressure_Pa = Pt9;
   
    Tt9 = S4.TotalTemperature_K;
    S9.TotalTemperature_K = Tt9;
    T9 = Tt9 * (P9/Pt9)^((gamma-1)/gamma);
    S9.Temperature_K = T9;
    
    M9 = (((Pt9/P9)^((gamma-1)/gamma) -1)*2/(gamma-1))^0.5;
    S9.Mach = M9;
    S9.Velocity_ms = S9.Mach *sqrt(gamma*(R_J_kmolK/MW_air)*S9.Temperature_K); % Solve for exit velocity

    %Solve for variable final flow Area
    S9.Area_m2 = S9.MassFlowRate_kgs * (R_J_kmolK/MW_air) * S9.Temperature_K / (S9.Pressure_Pa*S9.Velocity_ms); % m^2

end
