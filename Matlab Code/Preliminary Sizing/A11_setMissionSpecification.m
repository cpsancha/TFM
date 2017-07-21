%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                   MISSION ESPECIFICATIONS (ME)                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEAR CMD WINDOW & WORKSPACE
clearvars
close all
clc

%% MISSION TYPE:
%Choose one:
% ME.MissionType = 5 ; % Business Jets
ME.MissionType = 11; % Flying boats, amphibious, float airplanes


%% IMPORT CONSTANTS
CST = importConstants();


%% DEFINE FIGURES FOLDER AND POSITION
switch ME.MissionType
    case 5
        ME.FiguresFolder = '5_Figures';
    case 11
        ME.FiguresFolder = '11_Figures';
end



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
        ME.Crew = 2; %Check FAR 91.215 for minimun crew members
        %Mission crew weight in kg:
        ME.CrewWeight = ME.Crew*(CST.CrewWeightSI+CST.CrewBaggWeightSI);
        
    case 11 % amphibious
        ME.Crew = 2; %Check FAR 91.215 for minimun crew members
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
        ME.TakeOff.S_TOFL = 750;  % Take off distance in m
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
        ME.Cruise.Range    = DP.Range;          % in m
        ME.Cruise.Altitude = DP.CruiseAltitude; % in m
        ME.Cruise.Speed    = DP.CruiseSpeed;    % m/s
        
    case 11 % amphibious
        ME.Cruise.Range = 3500*1e3;     % in m
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
        ME.Landing.Altitude = 0;
end


%% POWERPLANT
%Choose one: 'propeller', 'jet'

switch ME.MissionType
    case  5 % business jet
        ME.Powerplant.Type = 'jet';
        
    case 11 % amphibious
        ME.Powerplant.Type = 'propeller';
        ME.Powerplant.Number = 2;
end


%% OTHERS
ME.Pressurization = NaN;
ME.Mission_Profile = NaN;


%% LOAD SIMILAR PLANES
switch ME.MissionType
    case  5 % business jet
        SP = importSimilarPlanes(ME.MissionType,CST,CF);
    case 11 % amphibious
        [SP, MV] = importSimilarPlanes(ME.MissionType,CST,CF);
end


%% CREATE DESIGN AIRCRAFT 
AC = aircraft();

%% FUSELAGE SHAPE

        AC.Fuselage.fusLength  = 30;
        AC.Fuselage.fusWidth   = 3;
        AC.Fuselage.fusHeight  = 3;


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
run C_weightEstimation.m
run D11_airplaneDesignParameters.m
run E11_configurationDesign.m