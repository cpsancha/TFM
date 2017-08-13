%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   MISSION ESPECIFICATIONS (ME)                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEAR CMD WINDOW & WORKSPACE
clear all %#ok<CLALL>
close all
clc


%% MISSION TYPE:
%Choose one:
ME.MissionType = 5 ; % Business Jets
% ME.MissionType = 11; % Flying boats, amphibious, float airplanes


%% IMPORT CONSTANTS
CST = importConstants();


%% DEFINE FIGURES FOLDER AND POSITION
switch ME.MissionType
    case 5
        ME.FiguresFolder = '5_Figures';
    case 11
        ME.FiguresFolder = '11_Figures';
end



%% CREATE ERRORS LIST
ME.errorList=cell(0);

%% DEFINE CONVERSION FACTORS
%length
    CF.nm2m     = 1852;
    CF.m2ft     = 3.28084;
    CF.ft2m     = 0.3048;
    CF.sm2m     = 1609.3; %Static miles to m
%speed
    CF.mph2ms   = 0.44704;
    CF.kts2ms   = 0.514444;
    CF.ms2kts   = 1.94384;
    CF.fps2kts  = 0.592484;
%mass
    CF.kg2lbm   = convmass(1,'kg','lbm');
    CF.lbm2kg   = convmass(1,'lbm','kg');
    CF.slug2kg  = 14.593903;
%force
    CF.lbf2N    = 4.448222;
    CF.N2lbf    = 0.224809;
%time
    CF.hour2sec = 3600;
    CF.min2sec  = 60;
%Pressure
    CF.psf2Pa   = CF.lbf2N/CF.ft2m^2;
%Power
    CF.hp2watts = 745.7;
%Other
    CF.TSFC2SI  = CF.lbm2kg/(CF.lbf2N*CF.hour2sec);
    CF.c_p2SI   = CF.lbm2kg*(1/CF.hp2watts)*(1/3600)*CST.GravitySI;
    
    

%% CREATE DESIGN AIRCRAFT 
AC = aircraft();



%% LOAD CUSTOM DESING PARAMETERS
run DesignParameters.m



%% PAYLOAD
switch ME.MissionType
    case  5 % business jet
        ME.Passengers = 6;
%         ME.CabinWeight = 1000; %Other weight to be added in kg
%         %Mission payload weight in kg:
%         ME.Payload =  ME.Passengers*(CST.PassengerWeightSI+CST.PassBaggWeightSI)+ME.CabinWeight;
        ME.Payload = DP.Payload; %kg --> Se saca en el excel, mejorar el dato cuando tenga todos los semejantes
        
    case 11 % amphibious
        ME.Payload = 3000;
end


%% CREW
switch ME.MissionType
    case  5 % business jet
        ME.Crew = DP.CrewNumber; %Check FAR 91.215 for minimun crew members
        %Mission crew weight in kg:
        ME.CrewWeight = ME.Crew*(CST.CrewWeightSI+CST.CrewBaggWeightSI);
        
    case 11 % amphibious
        ME.Crew = 2; %Check FAR 91.215 for minimun crew members
        AC.Payload.crew = ME.Crew;
        %Mission crew weight in kg:
        ME.CrewWeight = ME.Crew*(CST.CrewWeightSI+CST.CrewBaggWeightSI);
end


%% TAKE-OFF
switch ME.MissionType
    case  5 % business jet
        ME.TakeOff.Altitude = 0;     % Take off altitude in m
        ME.TakeOff.S_TOFL = DP.TOFL; % Take off distance in m
        
    case 11 % amphibious
        ME.TakeOff.Altitude = 0;  % Take off altitude in m
        ME.TakeOff.S_TOFL = 500;  % Take off distance in m
end


%% CLIMB
switch ME.MissionType
    case  5 % business jet
        ME.Climb.E_cl    = NaN;           %Time to climb in seconds --> se sobrescribe luego
        ME.Climb.Rate_cl = DP.ClimbRate;  %Rate of climb speed in m/s (vertical) (2500 ft/min)
        ME.Climb.V_cl    = DP.ClimbSpeed; %Climb speed in m/s (horizontal) (270 kts)
        
    case 11 % amphibious
        ME.Climb.E_cl    = NaN;  %Time to climb in seconds
        ME.Climb.Rate_cl = NaN;  %Rate of climb speed in m/s (vertical)
        ME.Climb.V_cl    = NaN;  %Climb speed in m/s (horizontal)
end


%% CRUISE
switch ME.MissionType
    case  5 % business jet
        ME.Cruise.Range    = DP.Range;          % in km
        ME.Cruise.Altitude = DP.CruiseAltitude; % in m
        ME.Cruise.Speed    = DP.CruiseSpeed;    % m/s

        
    case 11 % amphibious
        ME.Cruise.Range = 2000*1e3;     % in m
        ME.Cruise.Altitude =  6000;     % in m
        ME.Cruise.Speed = 138.8889;     % m/s
end



%% LOITER
switch ME.MissionType
    case  5 % business jet
        ME.Loiter.E_ltr = DP.LoiterTime;    %Loiter time in seconds
        ME.Loiter.V_ltr = ME.Cruise.Speed;  %Loiter speed in m/s
        
    case 11 % amphibious
        ME.Loiter.E_ltr = 30*CF.min2sec;    %Loiter time in seconds
        ME.Loiter.V_ltr = ME.Cruise.Speed;  %Loiter speed in m/s
end


%% ALTERNATE
switch ME.MissionType
    case 5 % business jet
        ME.Alternate.R_alt = DP.AlternateRange; %Range to alternate airport (200 nautic miles --> 370km)
        ME.Alternate.V_alt = 250*CF.kts2ms;     %If flight is below 10.000ft (typical) max speed is 250kts because of FAA regulations
        
    case 11 % amphibious
        ME.Alternate.R_alt = 370*1e3;       %Range to alternate airport (200 nautic miles --> 370km)
        ME.Alternate.V_alt = 250*CF.kts2ms; %If flight is below 10.000ft (typical) max speed is 250kts because of FAA regulations
        
end


%% LANDING
switch ME.MissionType
    case 5
        ME.Landing.S_LFL = DP.LFL; %Length of the landing field [m]
    case 11
        ME.Landing.S_LFL = 1000; %Length of the landing field [m]
end


%% POWERPLANT
%Choose one: 'propeller', 'jet'

switch ME.MissionType
    case  5 % business jet
        ME.Powerplant.Type = 'jet';
        
    case 11 % amphibious
        ME.Powerplant.Type = 'propeller';
end


%% OTHERS
ME.Pressurization  = NaN;
ME.Mission_Profile = NaN;





%% LOAD SIMILAR PLANES
switch ME.MissionType
    case  5 % business jet
        SP = importSimilarPlanes(ME.MissionType,CST,CF);
    case 11 % amphibious
%         SP = importSimilarPlanes(ME.MissionType,CST,CF);
end




%% EXAMPLE OF HOW TO OBTAIN MEAN VALUES FROM SIMILAR PLANES
% for i=1:length(SP)
%     if ~isempty(SP{i}.Actuations.Vcruise)
%         index(i)   = true; %#ok<SAGROW>
%     	vcruise(i) = SP{i}.Actuations.Vcruise; %#ok<SAGROW>
%     else
%         index(i)   = false; %#ok<SAGROW>
%         vcruise(i) = NaN; %#ok<SAGROW>
%     end
% end
% Vcruise = mean(vcruise(index));
% clear index vcruise i


%% CONTINUE...
run B_loadParameters.m
run C05_weightEstimation.m
run D05_airplaneDesignParameters.m
AC.Wing2.deltaCLdeltaE = 0;
run F05_wingConfiguration
globalOptions = optimoptions('fsolve', 'StepTolerance',1e-9, 'Display','none');
X0 = [DP.Incidence_1, DP.Incidence_2, DP.Stagger];
DP.ShowReportFigures  = false;
DP.ShowAircraftLayout = false;
[X,~,exitflag,~] = fsolve(@(X)getWingsIncidence(X, AC, ME, DP, Parameters, CST, CF),X0,globalOptions);
DP.Incidence_1 = X(1);
DP.Incidence_2 = X(2);
DP.Stagger     = X(3);
if ~isequal(exitflag,1)
    error('El solver que calcula las incidencias y el stagger no ha logrado converger correctamente. Se debería revisar el resultado.')
else
    clear exitflag  X0 X
end
X0 = [DP.fuselage_AoA, 0];
AircraftWeight = linspace(AC.Weight.MTOW*prod([Parameters.fuelFraction(1:6).value]),AC.Weight.MTOW*prod([Parameters.fuelFraction(1:4).value]),50);
DP.ShowReportFigures  = false;
DP.ShowAircraftLayout = false;
for i=1:length(AircraftWeight)
    [X,~,exitflag,~] = fsolve(@(X)trimAircraft(X, AircraftWeight(i), AC, ME, DP, Parameters, CST, CF),X0,globalOptions);
    if ~isequal(exitflag,1)
        error('El solver del trimado no ha logrado converger correctamente. Se debería revisar el resultado.')
    else
        clear exitflag
    end
    [~] = trimAircraft(X, AircraftWeight(i), AC, ME, DP, Parameters, CST, CF);
    Polar.Weight(i)       = AircraftWeight(i);
    Polar.Fuselage_AoA(i) = AC.Fuselage.fuselage_AoA;
    Polar.Wing1.CL(i)     = AC.Wing1.CL_wf;
    Polar.Wing1.Cma(i)    = AC.Wing1.Cm_ac_wf;
    Polar.Wing2.CL(i)     = AC.Wing2.CL_wf;
    Polar.Wing2.Cma(i)    = AC.Wing2.Cm_ac_wf;
    Polar.Wing2.deltaCLdeltaE(i) = AC.Wing2.deltaCLdeltaE;
    Polar.CL(i) = Parameters.q1_qinf * AC.Wing1.Sw/AC.Wing.Sw * AC.Wing1.CL_wf + Parameters.q2_qinf * AC.Wing2.Sw/AC.Wing.Sw * AC.Wing2.CL_wf;
    run G05_polarPrediction
    Polar.CD(i) = CD;
    Polar.D{i}  = D;
    Polar.DS{i} = DS;
    clear D DS
end
Polar.PolarFit = polyfit(Polar.CL,Polar.CD,2);
figure(); hold on;
plot(Polar.CL,Polar.CD,'r')
plot(linspace(min(Polar.CL),max(Polar.CL),50), polyval(Parameters.Polar.LongRangeCruise,linspace(min(Polar.CL),max(Polar.CL),50)),'b')
legend('Calculated Polar','Design Polar')
xlabel('C_L')
ylabel('C_D')
% AC.Fuselage.fuselage_AoA = 0;
% AC.Wing2.deltaCLdeltaE = 0;
% run F05_wingConfiguration
% run G05_polarPrediction



% test_values =[-5,-2.5,0,2.5,5];
% figure(10)
% hold on
% for i=1:length(test_values)
%     DP.TipTwist = test_values(i);
%     run B_loadParameters.m
%     run C05_weightEstimation.m
%     run D05_airplaneDesignParameters.m
%     run F05_wingConfiguration
%     figure(10); hold on
%     plot(eta1,Cl1)
% end
% legend('-5º','-2.5º','0º','2.5º','5º')

clear X X0 i exitflag globalOptions AircraftWeight

function [Error] = getWingsIncidence(X, AC, ME, DP, Parameters, CST, CF) %#ok<INUSL,INUSD>
  %INPUTS:  
    %X(1) = Wing1 Incidence [º]
    %X(2) = Wing2 Incidence [º]
    %X(3) = Stagger [m]
  %OUTPUTS:  
    %Error(1) = LiftCoeff - CL0 [-]
    %Error(2) = MomentCoeff [-]
    %Error(3) = Lift wing1 - 0.7Weight [kN]
    
  %Parse inputs
    DP.Incidence_1 = X(1);
    DP.Incidence_2 = X(2);
    DP.Stagger     = X(3);
  
  %Run wing's script
    run F05_wingConfiguration
    
  %Necessary calculation
    designWeight    = AC.Weight.EW + ME.Payload + AC.Weight.MFW/2;
    [rho,~,~,~,~,~] = atmos(DP.CruiseAltitude);
    CL0             = designWeight*CST.GravitySI / (0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw);
    
  %Parse outputs
    Error(1) =  AC.Wing.CL_wf - CL0;
    Error(2) =  AC.Wing.Cm_wf;
    Error(3) =  DP.Wing1_Wing2 - (AC.Wing1.Sw / AC.Wing.Sw * AC.Wing1.CL_w / CL0);
%     Error(3) = (DP.Wing1_Wing2*designWeight*CST.GravitySI - 0.5*rho*DP.CruiseSpeed^2*AC.Wing1.Sw*AC.Wing1.CL_wf)*1e-3;
    
end

function [Error] = trimAircraft(X, AircraftWeight, AC, ME, DP, Parameters, CST, CF)  %#ok<INUSD,INUSL>
  %Parse inputs
    AC.Fuselage.fuselage_AoA = X(1);
    AC.Wing2.deltaCLdeltaE   = X(2);
    
  %Run wing's script
    run F05_wingConfiguration
    
  %Necessary calculation
	[rho,~,~,~,~,~] = atmos(DP.CruiseAltitude);
 	CL0             = AircraftWeight*CST.GravitySI / (0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw); 
    
  %Parse outputs
    Error(1) =  AC.Wing.CL_wf - CL0;
    Error(2) =  AC.Wing.Cm_wf;   
end