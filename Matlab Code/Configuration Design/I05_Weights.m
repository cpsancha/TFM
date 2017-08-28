%************************************************************************************
%*																					*
%*					AIRCRAFT WEIGHTS CALCULATION SCRIPT								*
%*																					*
%************************************************************************************


%% WEIGHT SUBDIVISION:
options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');

%Speed definitions
[~, asound, P, rho] = atmosisa(DP.CruiseAltitude);
V_D     = DP.MaxSpeed;  %Design dive speed
V_D_EAS = correctairspeed(V_D, asound, P, 'TAS', 'EAS');
Vapp_SP   = sqrt((loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*DP.LFL);
VStall_L  = Vapp_SP/1.3; %[m/s]

% Ultimate Load
n_ult = getUltimateLoad(AC, AC.Wing, DP, CST, CF);

%AIRFRAME STRUCTURE:
    %Lead Wing group
    [WE.wing1,~,exitFlag,~]=fsolve(@(x)getWingWeight(x, AC, CST, AC.Wing1, V_D, VStall_L, n_ult, DP.delta_f),1e3,options);
    checkExitFlag(exitFlag);
    WE.wing1 = 0.75 * WE.wing1;
    
    %Rear Wing group
    [WE.wing2,~,exitFlag,~]=fsolve(@(x)getWingWeight(x, AC, CST, AC.Wing2, V_D, VStall_L, n_ult, 0),1e3,options);
    checkExitFlag(exitFlag);
    WE.wing2 = 0.75 * 0.60 * WE.wing2;
    
    % VTP
    WE.VTP = max( [0.64*(n_ult*AC.VTP.Sw^2)^0.75, 0.015*AC.Weight.EW] );

    %Body group
    lt = (AC.Wing2.Root_LE-0.25*AC.Wing2.RootChord) - (AC.Wing1.Root_LE-0.25*AC.Wing1.RootChord); %Distance, fig D-2
    WE.fuselage = 0.23*sqrt(V_D_EAS*lt/(AC.Fuselage.fusWidth+AC.Fuselage.fusHeight))*AC.Fuselage.Swet^1.2;
    WE.fuselage = WE.fuselage*1.08; % 8% Fuselage penality for pressurized cabins
    WE.fuselage = WE.fuselage*1.04; % 4% Fuselage penality for fuselage-mounted engines
    WE.fuselage = WE.fuselage*1.07; % 7% Fuselage penality for main landing gear attached to the fuselage
    WE.fuselage = WE.fuselage*0.85; % 15% Fuselage reduction as indicated by Roskam Part I table 2.16 for composites/metal
    WE.fuselage = WE.fuselage*0.80; % 25% Fuselage reduction as indicated by MTorres
    
    %Alighting gear group
    %[Nose, Main]
    A = [5.4, 15.0];
    B = [0.049, 0.033];
    C = [0, 0.021];
    D = [0, 0];
    WE.undercarriageNose = (1.08 * (A(1) + B(1)*AC.Weight.MTOW^(3/4) + C(1)*AC.Weight.MTOW + D(1)*AC.Weight.MTOW^(3/2)));
    WE.undercarriageMain = (1.08 * (A(2) + B(2)*AC.Weight.MTOW^(3/4) + C(2)*AC.Weight.MTOW + D(2)*AC.Weight.MTOW^(3/2)));
    
    % Surface control group
    WE.cockpit = 50; %50kg
    WE.autopilot = 9*AC.Weight.MTOW^(1/5);
    WE.systems = 0.008 * AC.Weight.MTOW;    %Todos estos valores son funcion del MTOW, realmente habria que iterar
    WE.surfaceControls = WE.cockpit + WE.autopilot + WE.systems;

    
%PROPULSION GROUP:
    %Engines
    k_pg  = 1.15; %Jet transports podded engines
    k_thr = 1.18; %With thrust reversers installed
    WE.engines = k_pg * k_thr * AC.Engine.Number * AC.Engine.Weight;
    
    % Engines installation and operation
    % Fuel systems
    % Thrust reversing probisions
    
    % Nacelles
    SwetNacelles = pi * AC.Engine.Diameter * 0.25*AC.Engine.Length;
    WE.nacelles = AC.Engine.Number * ...
                  (0.405 * sqrt(V_D_EAS) * SwetNacelles^1.3 + ...
                   0.05*AC.Engine.Weight + ...
                  (14.6+1.71+5.51)*SwetNacelles);
    

%AIRFRAME EQUIPMENT AND SERVICES
    % APUs
    k_apu  = 2;    %Installation factor 2-2.5
    k_wapu = 0.65; % Recent APU engines used on wide-body transports have a specific weight of only 65% of this value
    WE.APU = k_apu * k_wapu * 11.7 * (0.5*ME.HighDensityPassengers)^(3/5);

    % Instruments, navigational equipment & electronical groups
    k_ieg = 0.347;
    W_DE  = AC.Weight.EW; %kg Delivery empty weight
    R_D   = 10e3; %km maximum range
    WE.instruments = k_ieg * W_DE^(5/9)*R_D^(1/4);
    
    % Hidraulic
    WE.hidraulic = 0.011*W_DE + 181;
    
    % Electronics
    cabVol = AC.Fuselage.cabLength * pi/4 * AC.Fuselage.cabWidth * AC.Fuselage.cabHeight;
    if cabVol<227
        Pel = 0.565*cabVol;
    else
        Pel = 3.64*cabVol^0.7;
    end
    WE.electronics = 16.3*Pel*(1-0.033*sqrt(Pel));
    WE.electronics = 0.65*WE.electronics;
    
    % Air-conditioning & Anti-Icing
    WE.airConditioning = 14.0*AC.Fuselage.cabLength^1.28;
    
    
%FURNISHING & FIXED EQUIPMENT
    % Fixed ballast
    WE.flightDeck       = 16.5*W_DE^0.285;
    
    % Passenger seats
    WE.passengerSeats   = 22.7*ME.HighDensityPassengers; %pag 96 pdf torenbeek 22.7 kg, asientos vip
    WE.galleyProvisions = 29.5+45.3+113.4; %main meal gallery, snack pantry, coffee bar
    WE.toilet           = 136.0;
    
    % Floor covering
    WE.floorCovering    = 1.25*(AC.Fuselage.cabLength*AC.Fuselage.cabWidth)^1.15; %Scf = AC.Fuselage.cabLength*AC.Fuselage.cabWidth
    
    % Removable walls
    WE.misc             = 6.17*cabVol^1.07;
    
    % Basic emergency equipment
    WE.emergencyOxygen  = 18.1 + 1.09*ME.HighDensityPassengers;
    WE.emergencyFire    = 0.0012 *AC.Weight.MTOW;
    WE.emergencyEscape  = 0.453*(ME.Crew + ME.HighDensityPassengers);
 
    
% OPERATIONAL ITEMS & RENOVABLE EQUIPMENT
    WE.crewProvisions       = ME.CrewWeight; %ME.Crew*93 + 0*68; % 2 pilots,0 azafatas
    WE.passCabSupplies      = (8.62+2.27+2.27+2.27)*ME.HighDensityPassengers; %Long range first class with snacks
    WE.waterToiletChemicals = max([90.7,2.95*ME.HighDensityPassengers]); %Un solo baño
    WE.safetyEquipment      = 3.4*ME.HighDensityPassengers;

    % Fluids in a closed systems
    WE.residualFuel         = 0.151*(AC.Weight.MFW/CST.fuelDensity)^(2/3); 


% PAYLOAD
    WE.Payload = ME.Payload;
    
    
% TOTAL FUEL
WE.wing1Fuel    = AC.Wing1.fuelVolume * CST.fuelDensitySI; %[kg]
WE.wing2Fuel    = AC.Wing2.fuelVolume * CST.fuelDensitySI; %[kg]
WE.depositsFuel = DP.depositsUsage*(AC.Fuselage.cabLength*DP.depositsArea)*CST.fuelDensitySI; %[kg]
WE.tailFuel     = AC.Weight.MFW - WE.wing1Fuel - WE.wing2Fuel - WE.depositsFuel; %[kg]
totalFuelVolume = AC.Weight.MFW/CST.fuelDensitySI;   %[m^3]
fusFuelVolume   = WE.depositsFuel/CST.fuelDensitySI; %[m^3]
tailFuelVolume  = WE.tailFuel/CST.fuelDensitySI;     %[m^3]


% TOTAL WEIGHT
WE.Airframe    = WE.wing1 + WE.wing2 + WE.VTP + WE.fuselage + WE.undercarriageMain + WE.undercarriageNose + WE.surfaceControls;
WE.Propulsion  = WE.engines + WE.nacelles;
WE.Equipment   = WE.APU + WE.instruments + WE.hidraulic + WE.electronics + WE.airConditioning;
WE.Furnishing  = WE.flightDeck + WE.passengerSeats + WE.galleyProvisions + WE.toilet + ...
                 WE.floorCovering + WE.misc + WE.emergencyOxygen + WE.emergencyFire + WE.emergencyEscape;
WE.Operational = WE.crewProvisions + WE.passCabSupplies + WE.waterToiletChemicals + ...
                 WE.safetyEquipment + WE.residualFuel;
WE.OEW         = WE.Airframe + WE.Propulsion + WE.Equipment + WE.Operational;
WE.EW          = WE.OEW-WE.crewProvisions;
WE.MTOW        = WE.OEW + WE.Payload + AC.Weight.MFW;
disp(['  Calculated: MTOW=',num2str(WE.MTOW),'  EW=',num2str(WE.EW)])





%% LONGITUDINAL CENTER OF GRAVITY
%AIRFRAME
XCG.wing1             = AC.Wing1.Root_LE + 0.7*getChord(0.35*(AC.Wing1.WingSpan/2), AC.Wing1.WingSpan, AC.Wing1.TaperRatio, AC.Wing1.Sw);
XCG.wing2             = AC.Wing2.Root_LE + 0.7*getChord(0.35*(AC.Wing2.WingSpan/2), AC.Wing2.WingSpan, AC.Wing2.TaperRatio, AC.Wing2.Sw);
XCG.VTP               = AC.VTP.Root_LE + 0.42*getChord(0.38*(AC.Wing1.WingSpan/2), AC.Wing1.WingSpan, AC.Wing1.TaperRatio, AC.Wing1.Sw);
XCG.fuselage          = 0.47*AC.Fuselage.fusLength;
XCG.undercarriageNose = DP.UnderCarriageNose;
XCG.undercarriageMain = DP.UnderCarriageMain;
XCG.surfacecontrols   = ((WE.autopilot+WE.cockpit) * (AC.Fuselage.ln/2)+ ...
                         (WE.systems/2) * (AC.Wing1.CMA_LE + AC.Wing1.CMA) + ...
                         (WE.systems/2) * (AC.Wing2.CMA_LE + AC.Wing2.CMA))/WE.surfaceControls;

%PROPULSION
XCG.nacelles = DP.xEngine - AC.Engine.Length + 0.4*AC.Engine.Length;
XCG.engines  = DP.xEngine - AC.Engine.Length + 0.5*AC.Engine.Length;


%EQUIPMENT
XCG.Equipment = (WE.APU*AC.Fuselage.fusLength + ...
                (WE.instruments + WE.hidraulic + WE.electronics + WE.airConditioning) * ...
                0.5*AC.Fuselage.fusLength)/WE.Equipment;

%FURNISHING
XCG.Furnishing = (WE.flightDeck*(AC.Fuselage.ln/2) + ...
                  WE.passengerSeats*(AC.Fuselage.ln+AC.Fuselage.galleryLength+AC.Fuselage.passengersLength/2) + ...
                  WE.galleyProvisions*(AC.Fuselage.ln+AC.Fuselage.galleryLength/2) + ...
                  WE.toilet*(AC.Fuselage.ln+AC.Fuselage.cabLength-AC.Fuselage.lavatoryLength/2) + ...
                  WE.floorCovering*(AC.Fuselage.ln+AC.Fuselage.cabLength/2) + ...
                  WE.misc*(AC.Fuselage.ln+AC.Fuselage.cabLength/2) + ...
                 (WE.emergencyOxygen + WE.emergencyFire + WE.emergencyEscape)*(AC.Fuselage.ln+AC.Fuselage.cabLength/2))...
                 /WE.Furnishing;

%OPERATIONAL
XCG.Operational = (WE.crewProvisions*(AC.Fuselage.ln/2) + ...
                   WE.passCabSupplies*(AC.Fuselage.ln+AC.Fuselage.galleryLength/2) + ...
                   WE.waterToiletChemicals*(AC.Fuselage.ln+AC.Fuselage.cabLength-AC.Fuselage.lavatoryLength/2) + ...
                   WE.safetyEquipment*(AC.Fuselage.ln+AC.Fuselage.cabLength-AC.Fuselage.lavatoryLength/2) + ...
                   WE.residualFuel*(AC.Fuselage.fusLength/2))/WE.Operational;

%PAYLOAD
XCG.Payload = AC.Fuselage.ln+AC.Fuselage.galleryLength+AC.Fuselage.passengersLength/2;


%FUEL
XCG.Fuel = (WE.wing1Fuel*(AC.Wing1.CMA_LE+AC.Wing1.CMA/2) + ...
            WE.wing2Fuel*(AC.Wing2.CMA_LE+AC.Wing2.CMA/2) + ...
            WE.depositsFuel*(AC.Fuselage.ln+AC.Fuselage.cabLength/3) + ...
            WE.tailFuel*(AC.Fuselage.ln+AC.Fuselage.cabLength+AC.Fuselage.la/2))/AC.Weight.MFW;

               
%TOTAL
XCG.OEW = (XCG.wing1*WE.wing1 + XCG.wing2*WE.wing2 + XCG.VTP*WE.VTP + XCG.fuselage*WE.fuselage + ...
           XCG.undercarriageNose*WE.undercarriageNose + XCG.undercarriageMain*WE.undercarriageMain + ...
           XCG.surfacecontrols*WE.surfaceControls + XCG.engines*WE.engines + XCG.nacelles*WE.nacelles + ...
           XCG.Equipment*WE.Equipment + XCG.Furnishing*WE.Furnishing + XCG.Operational*WE.Operational)/WE.OEW;


%MTOW
XCG.MTOW = (XCG.OEW*WE.OEW + XCG.Payload*ME.Payload + XCG.Fuel*AC.Weight.MFW)/(WE.MTOW);
            
%ON FLIGHT
XCG.FLIGHT = (XCG.OEW*WE.OEW + XCG.Payload*ME.Payload + XCG.Fuel*AC.Weight.MFW/2)/(WE.OEW+WE.Payload+AC.Weight.MFW/2);
    
disp(['Xn=',num2str(AC.Weight.x_n),'  Xcg_OEW=',num2str(XCG.OEW),'  Xcg_MTOW=',num2str(XCG.MTOW),'  Xcg_FLIGHT=',num2str(XCG.FLIGHT)])

if DP.ShowReportFigures
    plotReportFigures(ME, WE, XCG);
end

clear A B C D options k_pg k_thr lt asound rho P V_D V_D_EAS Vapp_SP VStall_L SwetNacelles
clear SCautopilot SCcockpit SCsystems n_ult exitFlag k_apu k_wapu labels X k_ieg R_D W_DE Pel cabVol

%% USEFUL FUNCTIONS
function [] = checkExitFlag(exitFlag)
    if exitFlag ~= 1
        disp('Problema con la convergencia del solver')
        error('Problema con la convergencia')
    end
end

function [] = plotReportFigures(ME, WE, XCG)
        %AIRFRAME
    X = [WE.wing1, WE.wing2, WE.VTP, WE.fuselage,...
         WE.undercarriageMain+WE.undercarriageNose, WE.surfaceControls,...
         WE.engines, WE.nacelles];
    labels =  {'Ala delantera: ';...
               'Ala trasera: ';...
               'VTP: ';...
               'Fuselaje: ';...
               'Tren de aterrizaje: ';...
               'Sistemas de Control : ';...
               'Motores: ';...
               'Gondolas: '};
    drawCustomPieChart(X,labels,false,false);
    title('Contribuciones al peso estructural y planta propulsora','interpreter','latex')
    saveFigure(ME.FiguresFolder,'AirframePie') 
    
    %EQUIPMENT
    X = [WE.APU, WE.instruments, WE.hidraulic, WE.electronics, WE.airConditioning];
    labels =  {'APU: ';...
               'Equipos electronicos y navegacion: ';...
               'Equipos hidraulicos y neumaticos: ';...
               'Sistemas electricos: ';...
               'Aire acondicionado: '};
    drawCustomPieChart(X,labels,false,false);
    title('Contribuciones al peso de los equipos','interpreter','latex')
    saveFigure(ME.FiguresFolder,'EquipmentPie') 
    
    %FURNISHING
    X = [WE.flightDeck, WE.passengerSeats, WE.galleyProvisions, WE.toilet, WE.floorCovering, WE.misc, ...
         WE.emergencyOxygen + WE.emergencyFire + WE.emergencyEscape];
    labels =  {'Cabina de los pilotos: ';...
               'Asientos de los pasajeros: ';...
               'Cocina (Estructura y consumibles): ';...
               'Lavabo (Estructura y consumibles): ';...
               'Recubrimiento del suelo: ';...
               'Otros: ';...
               'Equipos de emergencia: '};
    drawCustomPieChart(X,labels,false,false);
    title('Contribuciones al peso del mobiliario','interpreter','latex')
    saveFigure(ME.FiguresFolder,'FurnishingPie') 
    
    %OPERATIONAL
    X = [WE.crewProvisions, WE.passCabSupplies, WE.waterToiletChemicals, WE.safetyEquipment, WE.residualFuel;];
    labels =  {'Tripulacion: ';...
               'Suministros de los pasajeros: ';...
               'Agua y productos quimicos: ';...
               'Equipo de seguridad: ';...
               'Aceite y combustible retenido: '};
    drawCustomPieChart(X,labels,false,false);
    title('Contribuciones al peso de los elementos operacionales','interpreter','latex')
    saveFigure(ME.FiguresFolder,'OperationalPie') 
    
    %TOTAL WEIGHT
    X = [WE.Propulsion, WE.Airframe, WE.Equipment, WE.Furnishing, WE.Operational];
    labels =  {'Peso de la planta propulsora: ';...
               'Peso estructural: ';...
               'Peso de los equipos: ';...
               'Peso del mobiliario: ';...
               'Peso de los elementos operacionales: '};
    drawCustomPieChart(X,labels,false,false);
    title('Contribuciones al peso vacio de la aeronave','interpreter','latex')
    saveFigure(ME.FiguresFolder,'WeightPie') 
    
    
    
    
       
    %OEW XCG
    X = [XCG.wing1*WE.wing1, XCG.wing2*WE.wing2, XCG.VTP*WE.VTP, XCG.fuselage*WE.fuselage, ...
           XCG.undercarriageNose*WE.undercarriageNose, XCG.undercarriageMain*WE.undercarriageMain, ...
           XCG.surfacecontrols*WE.surfaceControls, XCG.engines*WE.engines, XCG.nacelles*WE.nacelles, ...
           XCG.Equipment*WE.Equipment, XCG.Furnishing*WE.Furnishing, + XCG.Operational*WE.Operational];
    labels =  {'Ala 1: ';...
               'Ala 2: ';...
               'VTP: ';...
               'Fuselaje: ';...
               'Tren de morro: ';...
               'Tren principal: ';...
               'Superficies de Control: ';...
               'Motores: ';...
               'Gondolas: ';...
               'Equipos: ';...
               'Mobiliario: ';...
               'Elementos operacionales: '};
    drawCustomPieChart(X,labels,false,false);
    title('Distribucion de centros de gravedad','interpreter','latex')
    saveFigure(ME.FiguresFolder,'XcgOEWPie') 
end

function [f] = drawCustomPieChart(X,labels,explodeFlag,partialFlag)
    
    %Create figure
    f=figure();
    
    %Check if explode graph
    if explodeFlag
        explode = ones(1,length(X));
    else
        explode = zeros(1,length(X));
    end
    
    %Check if partial graph
    if ~partialFlag
        while sum(X)<1
            X=10.*X;
        end
    end
    
    %Create pie chart
    h=pie(X,explode);
    
    %Store Precalculated Percent Values
    hText = findobj(h,'Type','text'); % text object handles
    percentValues = get(hText,'String'); % percent values
    
    %Combine Percent Values and Additional Text
        %labels = {'Item A: ';'Item B: ';'Item C: '}; % strings
    combinedtxt = strcat(labels,percentValues); % strings and percent values
    
    %Store the text Extent property values for the current labels
    oldExtents_cell = get(hText,'Extent');  % cell array
    oldExtents = cell2mat(oldExtents_cell); % numeric array
    for i=1:length(X)
        hText(i).String = combinedtxt(i);
    end
    
    %Determine Horizontal Distance to Move Each Label
    newExtents_cell = get(hText,'Extent'); % cell array
    newExtents = cell2mat(newExtents_cell); % numeric array 
    width_change = newExtents(:,3)-oldExtents(:,3);
    signValues = sign(oldExtents(:,1));
    offset = signValues.*(width_change/2);
    
    %Position New Label
    textPositions_cell = get(hText,{'Position'}); % cell array
    textPositions = cell2mat(textPositions_cell); % numeric array
    textPositions(:,1) = textPositions(:,1) + offset; % add offset 
    for i=1:length(X)
        hText(i).Position = textPositions(i,:);
    end

end

function [n_ult] = getUltimateLoad(AC, Wing, DP, CST, CF)
    % Getting the limit load
    [~, asound, P, rho] = atmosisa(DP.CruiseAltitude);
    [~, ~, ~, rho0] = atmosisa(0);
    
    %Suponiendo todo el combustible en las alas
    Wdes = (AC.Weight.MTOW-AC.Weight.MFW)*CST.GravitySI;
    
    %Load Factor for limit maneuvering
        n = 2.1+(24e3/(10e3+AC.Weight.MTOW*CF.kg2lbm));
        if n<2.5
            n = 2.5;
        elseif n>3.8
            n = 3.8;
        end
        
    %Load gust Factor 
    mu = 2*(Wdes/Wing.Sw)/(rho*Wing.CL_alpha_wf*Wing.CMA*CST.GravitySI);
    kg = 0.88*mu/(5.3+mu);
    V_EAS =  correctairspeed(DP.CruiseSpeed, asound, P, 'TAS', 'EAS');
    if DP.CruiseAltitude*CF.m2ft<20e3
        U_EAS = 50*CF.ft2m;
    elseif DP.CruiseAltitude*CF.m2ft>=20e3
        U_EAS = ((50 + (DP.CruiseAltitude*CF.m2ft-20e3)*(25-50)/(50e3-20e3)))*CF.ft2m;
    end
    ngust = 1 + kg*Wing.CL_alpha_wf*rho0*U_EAS*V_EAS*Wing.Sw/(2*Wdes);
    %Most limiting:
    n_ult = 1.5*max([n, ngust]);

end

function [chord] = getChord(y, SpanWidth, TaperRatio, Sw)
    chord = (2*Sw/((1+TaperRatio)*SpanWidth)).*(1-(2*(1-TaperRatio)/SpanWidth).*abs(y));
end

function [Error] = getWingWeight(X, AC, CST, Wing, V_D, V_stall_L, n_ult, delta_f_L )
%Wing structural weight
%Inputs:
%	Wing: wing variable to get weight
%	V_D: Diving speed in EAS --> EAS=TASÂ·sqrt(rho/rho0)
%   Wdes: maximum weight with wing fuel tanks empty or maximum zero fuel weight

EstimatedTotalWeight = X;

%Suponiendo todo el combustible en las alas
Wdes = (AC.Weight.MTOW-AC.Weight.MFW)*CST.GravitySI;

% Skin joints factor
k_no = 1 + sqrt(1.905/Wing.WingSpan);

%Taper factor
k_taper = (1+Wing.TaperRatio)^0.4;

% Engines
k_e = 1.0; % 2 engines mounted 1 if  no engines

% Undercarriage
k_uc = 0.95; %no undercarriage

% Stifnnes against flutter
k_st = 1+9.06e-4*((Wing.WingSpan*cosd(Wing.Sweep_LE))^3/Wdes)*(V_D/100/Wing.Airfoil.t_c)^2*cosd(Wing.Sweep_12); % Low subsonic

% Struts 
k_b = 1; % no hay struts

%Basic Weight:
W_w_basic = 4.58e-3*k_no*k_taper*k_e*k_uc*k_st*...
    (k_b*n_ult*(Wdes-0.8*EstimatedTotalWeight))^0.55*...
    Wing.WingSpan^1.675*Wing.Airfoil.t_c^(-0.45)*(cosd(Wing.Sweep_12))^(-1.325);

%High lift devices:
k1 = 1.30;
	%k1=1.30: Double slotted Fowler
	%k1=1.45: Triple slotted Fowler
k2 = 1.25;
	%k2=1.00: Slotted flaps with fixed vane
	%k2=1.25: Double slotted flaps with variable geometry
kf = k1*k2;
bfs = 2*(Wing.eta_outboard-Wing.eta_inboard)*(Wing.WingSpan/2)/cosd(Wing.Sweep_RE);
W_tef = Wing.Swf*2.706*kf*(Wing.Swf*bfs)^(3/16)*...
    ((1.3*V_stall_L/100)^2*sind(delta_f_L)*cosd(Wing.Sweep_RE)/Wing.Airfoil.t_c)^(3/4);

%Spoilers:
W_sp = 0.015*W_w_basic;

% Total Weight:
TotalWeight = W_w_basic + 1.2*(1.15*W_tef+W_sp);

Error = TotalWeight-EstimatedTotalWeight;
end