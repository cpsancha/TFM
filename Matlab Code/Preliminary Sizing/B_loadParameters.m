%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                 	   LOAD PARAMETERS FROM ROSKAM                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CONVERSION FACTORS
% Movidos al apartado A para usarlos al importar los aviones semejantes


%% Load engine data base

switch ME.MissionType
    case 5
    case 11
        run('importTurbopropDatabase.m')
end
%% GENERATE A BUNCH OF DISTINGUISHABLE COLORS
Parameters.Colors = distinguishable_colors(20);


%% DYNAMIC PRESSURE RELATIONS
%Relaciones de presión dinámica segun torenbeek
    Parameters.q1_qinf = 1.00; 
    Parameters.q2_qinf = 0.85; %Torenbeek E-10.1
    

%% Table 2.2: Suggested Values For L/D, c_j, n_p and c_p for several mission phases:
switch ME.MissionType
    case 5
        %5.Business Jets
        % ROSKAM VALUES --> Old, so not used, instead obtained from actual similar planes
        Parameters.Cruise.L_D = [10, 12];
        Parameters.Cruise.c_j = [0.5, 0.9]*CF.TSFC2SI;  %lbm/lbf/hr to kg/(N·s)
        Parameters.Cruise.c_p = NaN;                    %Bussiness JETS can't be propeller driven
        Parameters.Cruise.n_p = NaN;                    %Bussiness JETS can't be propeller driven
        
        Parameters.Loiter.L_D = [12, 14];
        Parameters.Loiter.c_j = [0.4, 0.6]*CF.TSFC2SI;  %lbm/lbf/hr to kg/N/s
        Parameters.Loiter.c_p = NaN;                    %Bussiness JETS can't be propeller driven
        Parameters.Loiter.n_p = NaN;                    %Bussiness JETS can't be propeller driven
        
        %CUSTOM VALUES FROM SIMILAR PLANES
        Parameters.Cruise.L_D    = DP.CruiseEfficiency;         %[-]
        Parameters.Cruise.c_j    = DP.CruiseTSFC*CF.TSFC2SI;    %[kg/(N·s)]
        Parameters.Loiter.L_D    = DP.LoiterEfficiency;         %[-]
        Parameters.Loiter.c_j    = DP.LoiterTSFC*CF.TSFC2SI;    %[kg/(N·s)]
        Parameters.LowHeight.L_D = DP.LowHeightEfficiency;      %[-]
        Parameters.LowHeight.c_j = DP.LowHeightTSFC*CF.TSFC2SI; %[kg/(N·s)]
        
    case 11
        %11. Flying boats, amphibious, float airplanes
        Parameters.Cruise.L_D = [10, 12];
        Parameters.Cruise.c_j = [0.5, 0.9]*CF.TSFC2SI; %lbm/lbf/hr to kg/N/s
        Parameters.Cruise.c_p = [0.5, 0.7]*CF.c_p2SI;  %lbs/hp/hr to N/Watts/s
        Parameters.Cruise.n_p = 0.82;
        
        Parameters.Loiter.L_D = [13, 15];
        Parameters.Loiter.c_j = [0.4, 0.6]*CF.TSFC2SI; %lbm/lbf/hr to kg/N/s
        Parameters.Loiter.c_p = [0.5, 0.7]*CF.c_p2SI;  %lbs/hp/hr to N/Watts/s
        Parameters.Loiter.n_p = 0.77;
        
        %Custom Values
        Parameters.Cruise.L_D = 14;
        Parameters.Cruise.c_j = 0.5*CF.TSFC2SI; %lbm/lbf/hr to kg/N/s
        Parameters.Cruise.c_p = 0.46*CF.c_p2SI;  %lbs/hp/hr to N/Watts/s
        Parameters.Cruise.n_p = 0.82;
        
        Parameters.Loiter.L_D = 17;
        Parameters.Loiter.c_j = 0.4*CF.TSFC2SI; %lbm/lbf/hr to kg/N/s
        Parameters.Loiter.c_p = 0.46*CF.c_p2SI;  %lbs/hp/hr to N/Watts/s
        Parameters.Loiter.n_p = 0.77;
        
    otherwise
        warning('Unexpected mission type.')
end



%% Table 2.1: Suggested Fuel-Fractions for Several Mission Phases
switch ME.MissionType
    case 5
        %5.Business Jets
        Parameters.fuelFraction(1).phase = 'W1/W_TO: Engine start, warm-up';
        Parameters.fuelFraction(1).value = 0.990; %Roskam --> 0.990
        
        Parameters.fuelFraction(2).phase = 'W2/W1: Taxi';
        Parameters.fuelFraction(2).value = 0.995; %Roskam --> 0.995
        
        Parameters.fuelFraction(3).phase = 'W3/W2: Take-off';
        Parameters.fuelFraction(3).value = 0.995; %Roskam --> 0.995
        
        Parameters.fuelFraction(4).phase = 'W4/W3: Climb';
        Parameters.fuelFraction(4).value = 0.988; %Roskam --> 0.980
        %Can also be calculated from Breguet's equation
        %Hay que descontar del alcance ya que se recorre un cacho majo.
        %En el ejemplo habla de una velocidad de ascenso de 2500ft/min (12.7m/s)
        %mientras se avanza con velocidad constante de 275kts (510km/h), por lo que 
        %en el tiempo de subir a 12.500m (16 min) el avión recorrería 160 km
        ME.Climb.E_cl = ME.Cruise.Altitude/ME.Climb.Rate_cl; %aprox 900s (15min)
        
        Parameters.fuelFraction(5).phase = 'W5/W4: Cruise';
        Parameters.fuelFraction(5).value = NaN;  %From Breguet's equation
        
        Parameters.fuelFraction(6).phase = 'W6/W5: Loiter';
        Parameters.fuelFraction(6).value = NaN;  %From Breguet's equation
        
        Parameters.fuelFraction(7).phase = 'W7/W6: Descent';
        Parameters.fuelFraction(7).value = 0.990; %Roskam --> 0.990
        
        Parameters.fuelFraction(8).phase = 'W8/W7: Fly to alternate and descend';
        Parameters.fuelFraction(8).value = NaN;
        
        Parameters.fuelFraction(9).phase = 'W9/W8: Landing, taxi and shut-down';
        Parameters.fuelFraction(9).value = 0.992; %Roskam --> 0.992
        
    case 11
        %11. Flying boats, amphibious, float airplanes
        Parameters.fuelFraction(1).phase = 'W1/W_TO: Engine start, warm-up';
        Parameters.fuelFraction(1).value = 0.992;
        
        Parameters.fuelFraction(2).phase = 'W2/W1: Taxi';
        Parameters.fuelFraction(2).value = 0.990;
        
        Parameters.fuelFraction(3).phase = 'W3/W2: Take-off';
        Parameters.fuelFraction(3).value = 0.996;
        
        Parameters.fuelFraction(4).phase = 'W4/W3: Climb';
        Parameters.fuelFraction(4).value = 0.985;  %Can be calculated from Breguet's equation
        
        Parameters.fuelFraction(5).phase = 'W5/W4: Cruise';
        Parameters.fuelFraction(5).value = NaN;  %From Breguet's equation
        
        Parameters.fuelFraction(6).phase = 'W6/W5: Loiter';
        Parameters.fuelFraction(6).value = NaN;  %From Breguet's equation
        
        Parameters.fuelFraction(7).phase = 'W7/W6: Descent';
        Parameters.fuelFraction(7).value = 0.990;
        
        Parameters.fuelFraction(8).phase = 'W8/W7: Landing, taxi and shut-down';
        Parameters.fuelFraction(8).value = 0.990;
        
    otherwise
        warning('Unexpected mission type.')
end



%% Breguet´s equation to get fuel fraction during cruise and loiter:
%Updating structure 'parameters' with computed fuel fractions from Breguet's equation
switch ME.MissionType
    case 5
        %5.Business Jets
        Parameters = getFuelFraction('W5/W4: Cruise',ME,Parameters,CST,CF);
        Parameters = getFuelFraction('W6/W5: Loiter',ME,Parameters,CST,CF);
        Parameters = getFuelFraction('W8/W7: Fly to alternate and descend',ME,Parameters,CST,CF);
%         Mff = Parameters.fuelFraction(1).value * Parameters.fuelFraction(2).value * Parameters.fuelFraction(3).value * Parameters.fuelFraction(4).value * ...
%               Parameters.fuelFraction(5).value * Parameters.fuelFraction(6).value * Parameters.fuelFraction(7).value * Parameters.fuelFraction(8).value
    case 11
        %11. Flying boats, amphibious, float airplanes    
        Parameters = getFuelFraction('W5/W4: Cruise',ME,Parameters,CST,CF);
        Parameters = getFuelFraction('W6/W5: Loiter',ME,Parameters,CST,CF);
end



%% Table 2.15: Regresion constants A and B of equation 2.16:
%Regresion for MTOW and EW in lbf (IMPORTANTE)
switch ME.MissionType
    case 5
        %5.Business Jets
        Parameters.Table_2_15.a = 0.2678;
        Parameters.Table_2_15.b = 0.9979; 
        
    case 11
        %11. Flying boats, amphibious, float airplanes
        Parameters.Table_2_15.a = 0.1703;
        Parameters.Table_2_15.b = 1.0083;
end



%% Update Regresion Constants A and B with data from similar planes
switch ME.MissionType
    case 5
        %5.Business Jets
        [Parameters.Table_2_15.a,Parameters.Table_2_15.b] = getWTORegression(ME, SP, CST, CF, Parameters, DP );
    
    case 11
        %11. Flying boats, amphibious, float airplanes
        [Parameters.Table_2_15.a,Parameters.Table_2_15.b] = getWTORegression(ME, SP, CST, CF, Parameters, DP );
end



%% Table 3.1 Typical Values for Maximum Lift Coefficient
% Sobreescribir con semejantes?
switch ME.MissionType
    case 5
        % 5.Business Jets
        Parameters.CL_max    = [1.4, 1.8];
        Parameters.CL_max_TO = [1.6, 2.2];
        Parameters.CL_max_L  = [1.6, 2.6];
        
    case 11
        %11. Flying boats, amphibious, float airplanes
        Parameters.CL_max    = [1.2, 1.8];
        Parameters.CL_max_TO = [1.6, 2.2];
        Parameters.CL_max_L  = [1.8, 3.4];

end





%% Figures 3.21 Equivalent Skin Friction (cf)
switch ME.MissionType
    case 5
        % 5.Business Jets
        Parameters.cf = 0.0035; %fig 3.21c
    case 11
        %11. Flying boats, amphibious, float airplanes
        Parameters.cf = 0.005; %fig 3.21b
end



%% Table 3.5 Correlation coefficients for parasite area versus wetted area
% Performs interpolation in table 3.5 in function of cf.
cf = [0.009, 0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002];
a  = [-2.0458, - 2.0969, - 2.1549, -2.2218, - 2.3010, -2.3979, -2.5229, -2.6990];

Parameters.Table_3_4.a = interp1(cf,a,Parameters.cf);
Parameters.Table_3_4.b = 1;
clear cf a



%% Table 3.5 Regression line coefficients for take-off weight (in lbf) versus wetted area (in ft^2)
switch ME.MissionType
    case 5
        % 5.Business Jets
        Parameters.Table_3_5.c = 0.2263;
        Parameters.Table_3_5.d = 0.6977;
    case 11
        %11. Flying boats, amphibious, float airplanes
        Parameters.Table_3_5.c = 0.6295;
        Parameters.Table_3_5.d = 0.6708;
end



%% Table 3.6 First estimates for deltaCD0 and  'e' with flaps and gear down 
% Highly dependent on the used flaps and gear type. Split flaps are more draggy than Fowler flaps. Full span flaps are more draggy than partial span
% flaps. Wing mounted landing gears on high wing airplanesa are more draggy than those on low wing airplanes. More information in Roskam Part VI
switch ME.MissionType
    case 5
        %delta_CD0
        Parameters.Table_3_6.deltaC_D0.clean = 0;
        Parameters.Table_3_6.deltaC_D0.take_off_flaps = 0.015; %Roskam --> [0.010 , 0.020]
        Parameters.Table_3_6.deltaC_D0.landing_flaps  = 0.065; %Roskam --> [0.055 , 0.075]
        Parameters.Table_3_6.deltaC_D0.landing_gear   = 0.020; %Roskam --> [0.015 , 0.025]
        %e 
        Parameters.Table_3_6.e.clean          = 0.85; %Roskam --> [0.8 , 0.85]
        Parameters.Table_3_6.e.take_off_flaps =  0.8; %Roskam --> [0.75 , 0.8]
        Parameters.Table_3_6.e.landing_flaps  = 0.75; %Roskam --> [0.7 , 0.75]
        Parameters.Table_3_6.e.landing_gear   =  NaN; %Roskam --> No effect 
    case 11
        %delta_CD0
        Parameters.Table_3_6.deltaC_D0.clean = 0;
        Parameters.Table_3_6.deltaC_D0.take_off_flaps = 0.01; %[0.010 , 0.020];
        Parameters.Table_3_6.deltaC_D0.landing_flaps  = 0.06; %[0.055 , 0.075];
        Parameters.Table_3_6.deltaC_D0.landing_gear   = 0.015;%[0.015 , 0.025];
        %e
        Parameters.Table_3_6.e.clean          = 0.85; %[0.8 , 0.85];
        Parameters.Table_3_6.e.take_off_flaps = 0.8; %[0.75 , 0.8];
        Parameters.Table_3_6.e.landing_flaps  = 0.75; %[0.7 , 0.75];
        Parameters.Table_3_6.e.landing_gear   = NaN; %No effect
end



%% Weight Reduction due to the use of new composite materials and techniques
%Weight reduction of the empty weight as being fully manufatured in composite materials
%From: http://www.compositesworld.com/news/revolutionary-fuselage-concept-unveiled-by-mtorres --> Fuselage weight
%reduction estimated between 10% and 30%, and from Roskam Part V, Chapter 2 page 11, average FuselageWeight/EW=0.2
%so, total EW reduction estimated between 2% and 6%
switch ME.MissionType
    case 5
        Parameters.EWnew_EWold = DP.EWnew_EWold;
    case 11
        Parameters.EWnew_EWold = 0.94; 
end



%% Engine options
switch ME.MissionType
    case 5
        Parameters.EngineOptions = engines();
        % [ ] 2*Rolls-Royce AE 3007A1E --> Embraer Legacy family
        Parameters.EngineOptions(1).Model    = 'Rolls-Royce AE 3007A1E';
        Parameters.EngineOptions(1).Thrust   =    41; %[kN] at sea level
        Parameters.EngineOptions(1).Weight   = 751.6; %[kg] dry weight
        Parameters.EngineOptions(1).TSFC     = 0.625; %[lb/(lbf·h)] at cruise
        Parameters.EngineOptions(1).Diameter = 1.080; %[m]
        Parameters.EngineOptions(1).Length   = 1.900; %[m]
        
        % [ ] 2*Snecma Silvercrest --> Cessna Citation Hemisphere
        Parameters.EngineOptions(2).Model    = 'Snecma Silvercrest 2C';
        Parameters.EngineOptions(2).Thrust   =  53.4; %[kN] at sea level
        Parameters.EngineOptions(2).Weight   =  1040; %[kg] dry weight
        Parameters.EngineOptions(2).TSFC     = 0.597; %[lb/(lbf·h)] at cruise 0.95*0.628
        Parameters.EngineOptions(2).Diameter = 0.978; %[m]
        Parameters.EngineOptions(2).Length   = 2.705; %[m]
        
        % [X] 2*Snecma Silvercrest --> Dassault Falcon 5X
        Parameters.EngineOptions(3).Model    = 'Snecma Silvercrest 2D';
        Parameters.EngineOptions(3).Thrust   =  45.1; %[kN] at sea level (50.9kN)
        Parameters.EngineOptions(3).Weight   =  1040; %[kg] dry weight
        Parameters.EngineOptions(3).TSFC     = 0.597; %[lb/(lbf·h)] at cruise 0.95*0.628
        Parameters.EngineOptions(3).Diameter = 1.080; %[m]
        Parameters.EngineOptions(3).Length   = 1.900; %[m]
        
        % [ ] 2*Rolls-Royce BR710A2-20 --> Bombardier Global 5000
        Parameters.EngineOptions(4).Model  = 'Rolls-Royce BR710A2-20';
        Parameters.EngineOptions(4).Thrust = 65.6; %[kN] at sea level
        Parameters.EngineOptions(4).Weight = 1633; %[kg] dry weight
        Parameters.EngineOptions(4).TSFC   = 0.63; %[lb/(lbf·h)] at cruise
        
        % [ ] 2*P&W Canada PW800 --> Gulfstream G500/G600 (67.36/69.75 kN)
        Parameters.EngineOptions(5).Model  = 'P&W Canada PW814GA';
        Parameters.EngineOptions(5).Thrust = 67.36; %[kN] at sea level
        Parameters.EngineOptions(5).Weight =  1422; %[kg] dry weight
        Parameters.EngineOptions(5).TSFC   =    []; %[lb/(lbf·h)] at cruise
        
%        % [ ] 2*P&W PW1215G --> Mitsubishi Regional Jet
%        Parameters.EngineOptions(6).Model  = 'P&W PW1215G';
%        Parameters.EngineOptions(6).Thrust = 66.7; %[kN] at sea level
%        Parameters.EngineOptions(6).Weight =   []; %[kg] dry weight
%        Parameters.EngineOptions(6).TSFC   =   []; %[lb/(lbf·h)] at cruise
        
%         % [ ] 2*GE Passport 20–17BB1A (44 to 89 kN)       -->  ~90kN ~kg ~lb/(lbf·h)--> Bombardier Global 7000/8000 (73.4 kN)
%         Parameters.EngineOptions(7).Model   = 'GE Passport 20–17BB1A';
%         Parameters.EngineOptions(7).Thrust = 74.8; %[kN] at sea level
%         Parameters.EngineOptions(7).Weight = []; %[kg] dry weight
%         Parameters.EngineOptions(7).TSFC   = []; %[lb/(lbf·h)] at cruise
        
    case 11
end


%% AUXILIAR FUNCTIONS DEFINITION:
function [Parameters] = getFuelFraction( phase, ME, Parameters, CST, CF )
%DOC: Computes fuel fraction using Breguet's equation for range/endurance
%according to the mission phase, mission especification and propulsion paramaters.
switch phase
    case 'W5/W4: Cruise'
        %Cruise
        switch ME.Powerplant.Type
            case 'jet'
                cj  = Parameters.Cruise.c_j;    %Custom selected
                L_D = Parameters.Cruise.L_D;    %Custom selected
                ClimbRange = ME.Climb.E_cl*ME.Climb.V_cl; %The horizontal distance flight while climbing, is relevant [in meters]
                
                Parameters.fuelFraction(5).value=inv(exp((ME.Cruise.Range*1e3-ClimbRange)*CST.GravitySI*cj/(ME.Cruise.Speed*L_D)));
                
            case 'propeller'
                R   = ME.Cruise.Range;
                cp  = Parameters.Cruise.c_p(end);  %Worst case: maximo consumo
                np  = Parameters.Cruise.n_p(1);    %n_p is not an array, just in case
                L_D = Parameters.Cruise.L_D(1);    %Worst case: minima eficiencia
                
                Parameters.fuelFraction(5).value=1/exp(R*(cp/np)/L_D);
        end
        
    case 'W6/W5: Loiter'
        %Loiter
        switch ME.Powerplant.Type
            case 'jet'
                cj  = Parameters.Loiter.c_j;    %Custom selected
                L_D = Parameters.Loiter.L_D;    %Custom selected
                
                Parameters.fuelFraction(6).value=inv(exp(ME.Loiter.E_ltr*cj*CST.GravitySI/L_D));
                
            case 'propeller'
                cp  = Parameters.Loiter.c_p(end);  %Worst case: maximo consumo
                np  = Parameters.Loiter.n_p(1) ;   %n_p is not an array, just in case
                L_D = Parameters.Loiter.L_D(1);    %Worst case: minima eficiencia
                
                Parameters.fuelFraction(6).value=1/exp(ME.Loiter.E_ltr*ME.Loiter.V_ltr*(cp/np)/L_D);
        end
        
    case 'W8/W7: Fly to alternate and descend'
        %Fly to alternate and descend
        %NBAA IFR reserves are quoted as the amount of fuel required for the following profile:
        % (1) - a 5 minute approach at sea level
        % (2) - climb to 5000 feet
        % (3) - a 5 minute hold at 5000ft
        % (4) - climb to cruise altitude for the diversion to alternate
        % (5) - cruise at long range power
        % (6) - descend to sea level
        % (7) - land with 30 mins of holding fuel at 5000ft.
        switch ME.Powerplant.Type
            case 'jet'
                Cruise_cj  = Parameters.Cruise.c_j;          %Custom selected
                Cruise_L_D = Parameters.Cruise.L_D;          %Custom selected
                Loiter_cj  = Parameters.Loiter.c_j;          %Custom selected
                Loiter_L_D = Parameters.Loiter.L_D;          %Custom selected
                LowHeight_cj  = Parameters.LowHeight.c_j;    %Custom selected
                LowHeight_L_D = Parameters.LowHeight.L_D;    %Custom selected
                
                % (1) - a 5 minute approach at sea level               
                W10 = inv(exp(5*CF.min2sec*LowHeight_cj*CST.GravitySI/LowHeight_L_D));
                % (2) - climb to 5000 feet
                E_Climb5000ft = 5000*CF.ft2m/ME.Climb.Rate_cl;
                W21 = inv(exp(E_Climb5000ft*LowHeight_cj*CST.GravitySI/LowHeight_L_D));
                % (3) - a 5 minute hold at 5000ft
                W32 = inv(exp(5*CF.min2sec*Loiter_cj*CST.GravitySI/Loiter_L_D));
                % (4) - climb to cruise altitude for the diversion to alternate
                E_ClimbCruise = (ME.Cruise.Altitude-5000*CF.ft2m)/ME.Climb.Rate_cl;
                W43 = inv(exp(E_ClimbCruise*Cruise_cj*CST.GravitySI/Cruise_L_D)); %According Roskam, should be aprox 0.980
                % (5) - cruise at long range power
                ClimbRange = E_ClimbCruise*ME.Climb.V_cl; %The horizontal distance flight while climbing, is relevant [in meters]
                W54 = inv(exp((ME.Alternate.R_alt*1e3-ClimbRange)*CST.GravitySI*Cruise_cj/(ME.Cruise.Speed*Cruise_L_D)));
                % (6) - descend to sea level      
                W65 = 0.990; % Typical value from Roskam
                % (7) - land with 30 mins of holding fuel at 5000ft.
                W76 = inv(exp(30*CF.min2sec*Loiter_cj*CST.GravitySI/Loiter_L_D));
                %NBAA IFR reserves
                W70 = W76*W65*W54*W43*W32*W21*W10; %#ok<MINV> %Total weight ralation during the attempt landing, flight to alternate and touch ground with remaining fuel
                
                Parameters.fuelFraction(8).value = W70;
                
            case 'propeller'
                
        end

    otherwise
        warning('Wrong mission phase specifiec')
end
end

function [ A, B ] = getWTORegression(ME, SP, CST, CF, Parameters, DP )
%GETWTOGUESS Summary of this function goes here
%   Detailed explanation goes here

%Calculate regresion
for i=1:length(SP)
    if (~isempty(SP{i}.Weight.OEW))&&(~strcmp(SP{i}.Model,'Phenom 300'))
        Model(i) = SP{i}.Model;
        MTOW(i)  = SP{i}.Weight.MTOW; %#ok<*AGROW>
        EW(i)    = SP{i}.Weight.EW;
        index(i) = true;
    else
        Model(i) = '';
        MTOW(i)  = NaN;
        EW(i)    = NaN;
        index(i) = false;
    end
end


%Polyfit with weights in lbf, to obtain regresion constants according with Roskam
fitLbf = polyfit(log10(EW(index)*CST.GravitySI*CF.N2lbf),log10(MTOW(index)*CST.GravitySI*CF.N2lbf),1);
A = fitLbf(2);
B = fitLbf(1);


%Polyfit with weights in kg, for plotting
[fit, rsquared] = polyfitR2(log10(EW(index)),log10(MTOW(index)),1);


%Image customization
switch ME.MissionType
    case 5 %Business jet
        if DP.ShowReportFigures
        figure()
        loglog(EW(index),MTOW(index),'*','LineWidth',1,'Color',Parameters.Colors(1,:)); hold on;
        plot(linspace(min(EW)-0.5e3,max(EW)+1e3,2),10.^polyval(fit,linspace(log10(min(EW)-0.5e3),log10(max(EW)+1e3),2)),'LineWidth',1.25,'Color',Parameters.Colors(2,:));
        xlabel('EW$$\ [kg]$$','Interpreter','latex')
        ylabel('MTOW$$\ [kg]$$','Interpreter','latex')
        title('Logarithmic correlation between $$($$EW$$)$$ and $$($$MTOW$$)$$','Interpreter','latex')
        %Formating
        grafWidth   = 16;
        grafAR      = 0.6;
%         grid on
        xlim([9e3,25e3]);
        ylim([15e3,50e3]);
        set(gcf,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'FontSize',10,'FontName','Times new Roman','box','on')

        %Show equation in the graph
        x = min(EW) + 0.3*(max(EW)-min(EW));
        y = 10.^polyval(fit,log10(x));
        txt0 = '$$\ \ \leftarrow$$ $$\log_{10}($$MTOW$$)=A+B\log_{10}($$EW$$)$$';
        txt1 = strcat('$$\ \ \ \ \ \ A=',num2str(fit(2)),'$$');
        txt2 = strcat('$$\ \ \ \ \ \ B=',num2str(fit(1)),'$$');
        txt3 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(rsquared),'$$');
        text(x,y,txt0,'Interpreter','latex','FontSize',11)
        text(x,y-1.5e3,txt1,'Interpreter','latex','FontSize',11)
        text(x,y-2.8e3,txt2,'Interpreter','latex','FontSize',11)
        text(x,y-4.0e3,txt3,'Interpreter','latex','FontSize',11)
        %Image positioning
        %     set(gcf,'pos',[500   300   800   450])
        %         xlim([3.9,4.5])
        %         ylim([4.15,4.7])
        %Show models names
        for i=1:length(Model)
%             if index(i)
%                 if strcmp(Model(i),'Falcon 7X')
%                     text(EW(i)+0.3e3,MTOW(i),char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Falcon 900LX')
%                     text(EW(i)+0.2e3,MTOW(i),char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Falcon 2000LXS')
%                     text(EW(i)+0.15e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Global 6000')
%                     text(EW(i)-3.3e3,MTOW(i),char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Global 5000')
%                     text(EW(i)-3.5e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Challenger 350')
%                     text(EW(i)-1.9e3,MTOW(i)+0.2e3,char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Gulfstream G280')
%                     text(EW(i)+0.1e3,MTOW(i)-0.1e3,char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Legacy 450')
%                     text(EW(i)+0.1e3,MTOW(i)-0.15e3,char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Legacy 500')
%                     text(EW(i)-1.45e3,MTOW(i)+0.1e3,char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Phenom 300')
%                     text(EW(i),MTOW(i),char(Model(i)),'FontSize',9)
%                 elseif strcmp(Model(i),'Citation X+')
%                     text(EW(i)+0.1e3,MTOW(i),char(Model(i)),'FontSize',9)
%                 else
%                     text(EW(i),MTOW(i),char(Model(i)),'FontSize',9)
%                 end
%             end
            if index(i)
                if strcmp(Model(i),'Falcon 7X')
                    text(EW(i)+0.3e3,MTOW(i),char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Falcon 900LX')
                    text(EW(i)+0.2e3,MTOW(i),char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Falcon 2000LXS')
                    text(EW(i)+0.15e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Global 6000')
                    text(EW(i)-2.9e3,MTOW(i),char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Global 5000')
                    text(EW(i)-2.9e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Challenger 350')
                    text(EW(i)-1.7e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Gulfstream G280')
                    text(EW(i)+0.1e3,MTOW(i)-0.1e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Legacy 450')
                    text(EW(i)+0.1e3,MTOW(i)-0.15e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Legacy 500')
                    text(EW(i)-1.30e3,MTOW(i)+0.1e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Phenom 300')
                    text(EW(i),MTOW(i),char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Citation X+')
                    text(EW(i)+0.1e3,MTOW(i)+0.05e3,char(Model(i)),'FontSize',8)
                else
                    text(EW(i),MTOW(i),char(Model(i)),'FontSize',8)
                end
            end
        end
        p=findobj(gca,'Type','text');
        for i=1:length(p)
            uistack(p(i), 'top')
        end
        clear p
        %Save Figure
        saveFigure(ME.FiguresFolder,'EW_MTOW_Correlation')
        
        end
        
        
    case 11 %Amphibious
    figure()
    loglog(EW(index),MTOW(index),'*','LineWidth',1,'Color',Parameters.Colors(1,:)); hold on;
    plot(linspace(min(EW),max(EW),2),10.^polyval(fit,linspace(min(log10(EW)),max(log10(EW)),2)),'LineWidth',1.25,'Color',Parameters.Colors(2,:));
    xlabel('EW$$\ [kg]$$','Interpreter','latex')
    ylabel('MTOW$$\ [kg]$$','Interpreter','latex')
    title('Logarithmic correlation between $$($$EW$$)$$ and $$($$MTOW$$)$$','Interpreter','latex')
        
        
        %Formating
        grafWidth   = 16;
        grafAR      = 0.6;
% %         grid on
        xlim([9e3,35e3]);
        ylim([15e3,55e3]);
        set(gcf,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'FontSize',10,'FontName','Times new Roman','box','on')

        %Show equation in the graph
        x = min(EW) + 0.3*(max(EW)-min(EW));
        y = 10.^polyval(fit,log10(x));
        txt0 = '$$\ \ \leftarrow$$ $$\log_{10}($$MTOW$$)=A+B\log_{10}($$EW$$)$$';
        txt1 = strcat('$$\ \ \ \ \ \ A=',num2str(fit(2)),'$$');
        txt2 = strcat('$$\ \ \ \ \ \ B=',num2str(fit(1)),'$$');
        txt3 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(rsquared),'$$');
        text(x,y,txt0,'Interpreter','latex','FontSize',11)
        text(x,y-1.5e3,txt1,'Interpreter','latex','FontSize',11)
        text(x,y-2.8e3,txt2,'Interpreter','latex','FontSize',11)
        text(x,y-4.0e3,txt3,'Interpreter','latex','FontSize',11)
        for i=1:length(Model)

            if index(i)
                if strcmp(Model(i),'Bombardier 415')
                    text(EW(i),MTOW(i)+0.5e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Grumann Albatross')
                    text(EW(i)+0e3,MTOW(i)-0.4e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'ShinMaywa US-1A')
                    text(EW(i)-6e3,MTOW(i)-0.25e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Global 6000')
                    text(EW(i)-2.9e3,MTOW(i),char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Global 5000')
                    text(EW(i)-2.9e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Challenger 350')
                    text(EW(i)-1.7e3,MTOW(i)+0.25e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Gulfstream G280')
                    text(EW(i)+0.1e3,MTOW(i)-0.1e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Legacy 450')
                    text(EW(i)+0.1e3,MTOW(i)-0.15e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Legacy 500')
                    text(EW(i)-1.30e3,MTOW(i)+0.1e3,char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Phenom 300')
                    text(EW(i),MTOW(i),char(Model(i)),'FontSize',8)
                elseif strcmp(Model(i),'Citation X+')
                    text(EW(i)+0.1e3,MTOW(i)+0.05e3,char(Model(i)),'FontSize',8)
                else
                    text(EW(i),MTOW(i),char(Model(i)),'FontSize',8)
                end
            end
        end
        
        
                %Save Figure
        saveFigure(ME.FiguresFolder,'EW_MTOW_Correlation')
end


end


