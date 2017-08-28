

fsolve(@(x)getFinalWeight(x,AC,ME,CST,CF,AF,delta_f_L), [20000,18000])
%% Getting the limit load
%Suponiendo todo el combustible en las alas
Wdes = (AC.Weight.MTOW-AC.Weight.FW)*CST.GravitySI;
%Load Factor for limit maneuvering
    n = 2.1+(24e3/(10e3+AC.Weight.MTOW*CF.kg2lbm));
    if n<2.5
        n = 2.5;
    elseif n>3.8
        n = 3.8;
    end
%Load gust Factor 
mu = 2*(Wdes/AC.Wing.Sw)/(ME.Cruise.Density*AC.Wing1.CL_alpha_wf*AC.Wing1.CMA*CST.GravitySI);
kg = 0.88*mu/(5.3+mu);
[~, asound, P, ~] = atmosisa(ME.Cruise.Altitude);
[~, ~, ~, rho0] = atmosisa(0);
V_EAS =  correctairspeed(ME.Cruise.Speed, asound, P, 'TAS', 'EAS');
    if ME.Cruise.Altitude*CF.m2ft<20000
    U_EAS = 50*CF.ft2m;
    elseif ME.Cruise.Altitude*CF.m2ft>=20001
    U_EAS = ((50+ (ME.Cruise.Altitude*CF.m2ft-20000)*(25-50)/(50000-20000)))*CF.ft2m;
    end
ngust = 1+ kg*rho0*AC.Wing1.CL_alpha_wf*0.5*U_EAS*V_EAS*AC.Wing.Sw/Wdes;
%Most limiting:
n_ult = 1.5*max([n, ngust]);

%% Getting the main weights:
% Fuselage
V_D = ME.Cruise.Speed/0.8;  %Design dive speed
lt = AC.Wing2.Root_LE - (AC.Wing1.Root_LE+0.25*AC.Wing1.CMA);
WE.fuselage = 0.23*sqrt(V_D*lt/(AC.Fuselage.fusWidth+AC.Fuselage.fusHeight))*AC.Fuselage.Swet^1.2;
WE.fuselage = WE.fuselage*0.8*0.7; %MTorres reduction

% Hull
WE.Hull = AC.Weight.EW * 0.12;
WE.Hull = 0.5*WE.Hull; % Composite weight reduction

% Tip floats
WE.TipFloats = 2 * 0.012 * 0.5 * AC.Weight.EW;

% Wings
WE.Wing1 = fsolve(@(x)getWing1Weight(x,AC, ME, AF,delta_f_L,n_ult, Wdes, WE),1000);
WE.Wing1 = WE.Wing1*0.7;
WE.Wing2 = fsolve(@(x)getWing2Weight(x,AC, ME, AF,delta_f_L,n_ult, Wdes, WE),1000);
WE.Wing2 = WE.Wing2*0.7;

% VTP
WE.VTP = 0.64 * (n_ult*AC.VTP.Swet^2)^0.75;

% Surface control group
cockpit = 0.046*AC.Weight.MTOW^(3/4);
ap = 9*AC.Weight.MTOW^(1/5);
system = 0.008 * AC.Weight.MTOW;    %Todos estos valores son funcion del MTOW, realmente habria que iterar
WE.Surface_controls = cockpit+ ap +system;      

% Engines
WE.Engines = 1.35 * ME.Powerplant.Number * (AC.Engine.Weight +0.109*AC.Engine.Power/CF.hp2watts);

% Nacelles
WE.Nacelles = 0.0635 * ME.Powerplant.Number*AC.Engine.Power/CF.hp2watts;

% Services 
APU = 2 * 11.7 * (0.5*ME.Passengers)^(3/5);
NAV_COM = 54.4 + 9.1*ME.Powerplant.Number + 0.006 * AC.Weight.MTOW;
hidra_elec = 0.277* AC.Weight.EW^(4/5);
airConditioning = 14 * AC.Fuselage.cabLength^1.28;
WE.Services = APU+NAV_COM+hidra_elec+airConditioning;

% Furnishing
flightDeck = 9.1*AC.Weight.EW^0.285;
passengerSeats = 22.7*ME.Passengers; %pag 96 pdf torenbeek 22.7 kg, asientos vip
galleyProvisions = 29.5+45.3+113.4; %main meal gallery, snack pantry, coffee bar
toilet = 75 + 38.5;
floorCovering = 0.94*(AC.Fuselage.cabLength*AC.Fuselage.cabWidth)^1.15; %Scf = AC.Fuselage.cabLength*AC.Fuselage.cabWidth
misc = 3.69*(ME.Cargo.Volume+AC.Fuselage.cabVolume)^1.14; % Cargo volume + Cabin volume
cargoRestraints = 1.28*ME.Cargo.Volume;
cargoHandling = 13.67*ME.Cargo.FloorArea;
WE.Furnishing = flightDeck+passengerSeats+galleyProvisions+toilet+floorCovering+...
    misc+cargoRestraints+cargoHandling;

% Operational items 
crewProvisions = 93*0+150*ME.Crew; % 2 pilots,0 azafatas
passengerCabinSupplies = (8.62+2.27)*ME.Passengers; %Long range first class
potableWaterToiletChemicals = 90.7; %Un solo baño
safetyEquipment = 3.4*ME.Passengers;
residualFuel = 0.151*(AC.Weight.FW/CST.fuelDensity)^(2/3); 
emergencyOxygen = (18.1 + 1.09*ME.Passengers);
emergencyFire = 0.003 *AC.Weight.MTOW;
emergencyEscape = 0.453*ME.Passengers;
WE.Operational = crewProvisions+passengerCabinSupplies+potableWaterToiletChemicals+...
    safetyEquipment+residualFuel+emergencyOxygen+emergencyFire+emergencyEscape;

AC.Weight.EW = WE.fuselage + WE.Hull + WE.TipFloats + WE.Wing1 + WE.Wing2+ WE.VTP + WE.Surface_controls +...
    WE.Engines + WE.Nacelles + WE.Services + WE.Furnishing + WE.Operational;
AC.Weight.MTOW = AC.Weight.EW + AC.Weight.FW + ME.Payload;

%% XCG
XCG.fuselage = 0.39 * AC.Fuselage.fusLength;
XCG.Hull = XCG.fuselage;
XCG.TipFloats = AC.Wing1.Root_LE +  AC.Wing1.RootChord*0.5;
XCG.Wing1 = AC.Wing1.Root_LE +0.4 * AC.Wing1.RootChord;
XCG.Wing2 = AC.Wing2.Root_LE +0.4 * AC.Wing2.RootChord;
XCG.VTP   = AC.VTP.Root_LE + 0.42 * AC.VTP.RootChord;
XCG.Surface_controls = (0.5*system*(AC.Wing1.Root_LE+AC.Wing1.CMA) + 0.5*system*(AC.Wing2.Root_LE+AC.Wing2.CMA)+...
    AC.Fuselage.cabinFrac*AC.Fuselage.fusLength * (ap + cockpit))/WE.Surface_controls;
XCG.Nacelles = AC.Wing1.Root_LE + AC.Wing1.RootChord - AC.Engine.Length + 0.4*AC.Engine.Length;
XCG.Engines = AC.Wing1.Root_LE + AC.Wing1.RootChord - AC.Engine.Length + 0.5*AC.Engine.Length;
% Random
XCG.Services = (APU*AC.Fuselage.fusLength + NAV_COM*AC.Fuselage.cabinFrac*AC.Fuselage.fusLength+ ...
    hidra_elec*AC.Fuselage.fusLength/2 + airConditioning*AC.Fuselage.fusLength/2)/WE.Services;
XCG.Furnishing = 0.7*AC.Fuselage.fusLength;
XCG.Operational = 0.5*AC.Fuselage.fusLength;

XCG.Total = (XCG.fuselage*WE.fuselage + XCG.Hull*WE.Hull + XCG.TipFloats*WE.TipFloats+...
    XCG.Wing1*WE.Wing1 + XCG.Wing2*WE.Wing2 + XCG.VTP*WE.VTP + XCG.Surface_controls*WE.Surface_controls+...
    XCG.Nacelles*WE.Nacelles + XCG.Engines*WE.Engines + XCG.Services*WE.Services +...
    XCG.Furnishing*WE.Furnishing + XCG.Operational*WE.Operational)/AC.Weight.EW;
%%

function [F,WE] = getFinalWeight(x,AC,ME,CST,CF,AF,delta_f_L)
AC.Weight.MTOW = x(1);
AC.Weight.EW = x(2);

%% Getting the limit load
%Suponiendo todo el combustible en las alas
Wdes = (AC.Weight.MTOW-AC.Weight.FW)*CST.GravitySI;
%Load Factor for limit maneuvering
    n = 2.1+(24e3/(10e3+AC.Weight.MTOW*CF.kg2lbm));
    if n<2.5
        n = 2.5;
    elseif n>3.8
        n = 3.8;
    end
%Load gust Factor 
mu = 2*(Wdes/AC.Wing.Sw)/(ME.Cruise.Density*AC.Wing1.CL_alpha_wf*AC.Wing1.CMA*CST.GravitySI);
kg = 0.88*mu/(5.3+mu);
[~, asound, P, ~] = atmosisa(ME.Cruise.Altitude);
[~, ~, ~, rho0] = atmosisa(0);
V_EAS =  correctairspeed(ME.Cruise.Speed, asound, P, 'TAS', 'EAS');
    if ME.Cruise.Altitude*CF.m2ft<20000
    U_EAS = 50*CF.ft2m;
    elseif ME.Cruise.Altitude*CF.m2ft>=20000
    U_EAS = ((50+ (ME.Cruise.Altitude*CF.m2ft-20000)*(25-50)/(50000-20000)))*CF.ft2m;
    end
ngust = 1+ kg*rho0* AC.Wing1.CL_alpha_wf*0.5*U_EAS*V_EAS*AC.Wing.Sw/Wdes;
%Most limiting:
n_ult = 1.5*max([n, ngust]);

%% Getting the main weights:
% Fuselage
V_D = ME.Cruise.Speed/0.8;  %Design dive speed
lt = AC.Wing2.Root_LE - (AC.Wing1.Root_LE-0.25*AC.Wing1.CMA);
WE.fuselage = 0.23*sqrt(V_D*lt/(AC.Fuselage.fusWidth+AC.Fuselage.fusHeight))*AC.Fuselage.Swet^1.2;
WE.fuselage = WE.fuselage*0.8*0.7; %MTorres reduction

% Hull
WE.Hull = AC.Weight.EW * 0.12;
WE.Hull = 0.5*WE.Hull; % Composite weight reduction

% Tip floats
WE.TipFloats = 2 * 0.012 * 0.5 * AC.Weight.EW;

% Wings
WE.Wing1 = fsolve(@(x)getWing1Weight(x,AC, ME, AF,delta_f_L,n_ult, Wdes, WE),1000);
WE.Wing1 = WE.Wing1*0.7;
WE.Wing2 = fsolve(@(x)getWing2Weight(x,AC, ME, AF,delta_f_L,n_ult, Wdes, WE),1000);
WE.Wing2 = WE.Wing2*0.7;

% VTP
WE.VTP = 0.64 * (n_ult*AC.VTP.Swet^2)^0.75;

% Surface control group
cockpit = 0.046*AC.Weight.MTOW^(3/4);
ap = 9*AC.Weight.MTOW^(1/5);
system = 0.008 * AC.Weight.MTOW;    %Todos estos valores son funcion del MTOW, realmente habria que iterar
WE.Surface_controls = cockpit+ ap +system;       clear cockpit ap system

% Engines
WE.Engines = 1.35 * ME.Powerplant.Number * (AC.Engine.Weight +0.109*AC.Engine.Power/CF.hp2watts);

% Nacelles
WE.Nacelles = 0.0635 * ME.Powerplant.Number*AC.Engine.Power/CF.hp2watts;

% Services 
APU = 2 * 11.7 * (0.5*ME.Passengers)^(3/5);
NAV_COM = 54.4 + 9.1*ME.Powerplant.Number + 0.006 * AC.Weight.MTOW;
hidra_elec = 0.277* AC.Weight.EW^(4/5);
airConditioning = 14 * AC.Fuselage.cabLength^1.28;
WE.Services = APU+NAV_COM+hidra_elec+airConditioning;

% Furnishing
flightDeck = 9.1*AC.Weight.EW^0.285;
passengerSeats = 22.7*ME.Passengers; %pag 96 pdf torenbeek 22.7 kg, asientos vip
galleyProvisions = 29.5+45.3+113.4; %main meal gallery, snack pantry, coffee bar
toilet = 75 + 38.5;
floorCovering = 0.94*(AC.Fuselage.cabLength*AC.Fuselage.cabWidth)^1.15; %Scf = AC.Fuselage.cabLength*AC.Fuselage.cabWidth
misc = 3.69*(ME.Cargo.Volume+AC.Fuselage.cabVolume)^1.14; % Cargo volume + Cabin volume
cargoRestraints = 1.28*ME.Cargo.Volume;
cargoHandling = 13.67*ME.Cargo.FloorArea;
WE.Furnishing = flightDeck+passengerSeats+galleyProvisions+toilet+floorCovering+...
    misc+cargoRestraints+cargoHandling;

% Operational items 
crewProvisions = 93*0+150*ME.Crew; % 2 pilots,0 azafatas
passengerCabinSupplies = (8.62+2.27)*ME.Passengers; %Long range first class
potableWaterToiletChemicals = 90.7; %Un solo baño
safetyEquipment = 3.4*ME.Passengers;
residualFuel = 0.151*(AC.Weight.FW/CST.fuelDensity)^(2/3); 
emergencyOxygen = (18.1 + 1.09*ME.Passengers);
emergencyFire = 0.003 *AC.Weight.MTOW;
emergencyEscape = 0.453*ME.Passengers;
WE.Operational = crewProvisions+passengerCabinSupplies+potableWaterToiletChemicals+...
    safetyEquipment+residualFuel+emergencyOxygen+emergencyFire+emergencyEscape;

AC.Weight.EW = WE.fuselage + WE.Hull + WE.TipFloats + WE.Wing1 + WE.Wing2+ WE.VTP + WE.Surface_controls +...
    WE.Engines + WE.Nacelles + WE.Services + WE.Furnishing + WE.Operational;
AC.Weight.MTOW = AC.Weight.EW + AC.Weight.FW + ME.Payload;

F(1) = x(1) - AC.Weight.MTOW;
F(2) = x(2) - AC.Weight.EW;

end
%% Auxiliar functions
%% Wing basic weight 
function F = getWing1Weight(x,AC, ME, AF,delta_f_L,n_ult, Wdes, WE)
WE.Wing1 = x;
constant = 4.58e-3;
% Skin joints factor
k_no = 1 + sqrt(1.905/AC.Wing1.WingSpan);
%Taper factor
k_taper = (1+AC.Wing1.TaperRatio)^0.4;
% Engines
k_e = 0.95; % 2 engines mounted 1 if  no engines
% Undercarriage
k_uc = 0.95; %no undercarriage
% Stifnnes against flutter
k_st = 1; % Low subsonic
% Struts 
k_b = 1; % no hay struts

%Basic Weight:
W_w_basic = constant*k_no*k_taper*k_e*k_uc*k_st*...
    (k_b*n_ult*(Wdes-0.8 *WE.Wing1))^0.55*...
    AC.Wing1.WingSpan^1.675*AF.t_c^(-0.45)*(cos(AC.Wing1.Sweep_12*pi/180))^(-1.325); %OJO!!!! El "+" de la linea anterior es un "*"

%High lift devices:
kf =1.45*1.25; % El mas tocho
bfs = 2*(AC.Wing1.eta_outboard-AC.Wing1.eta_inboard)*AC.Wing1.WingSpan/2;
W_tef = AC.Wing1.Swf*2.706*kf*(AC.Wing1.Swf*bfs)^(3/16)*...
    ((1.8*ME.Landing.Speed/100)^2*sin(delta_f_L*pi/180)/AF.t_c)^(3/4);
%Spoilers:
W_sp = 0.015*W_w_basic; 

% Total Weight:
WE.Wing1 = W_w_basic + 1.2*(W_tef+W_sp);

F = WE.Wing1-x;
end
function F = getWing2Weight(x,AC, ME, AF,delta_f_L,n_ult, Wdes, WE)
WE.Wing2 = x;
constant = 4.58e-3;
% Skin joints factor
k_no = 1 + sqrt(1.905/AC.Wing2.WingSpan);
%Taper factor
k_taper = (1+AC.Wing2.TaperRatio)^0.4;
% Engines
k_e = 1; % 0.95==2 engines mounted ;1 if  no engines
% Undercarriage
k_uc = 0.95; %no undercarriage
% Stifnnes against flutter
k_st = 1; % Low subsonic
% Struts 
k_b = 1; % no hay struts

%Basic Weight:
W_w_basic = constant*k_no*k_taper*k_e*k_uc*k_st*...
    (k_b*n_ult*(Wdes-0.8*WE.Wing2))^0.55*...
    AC.Wing2.WingSpan^1.675*AF.t_c^(-0.45)*(cos(AC.Wing2.Sweep_12*pi/180))^(-1.325);

%High lift devices:
kf =1.45*1.25; % El mas tocho
bfs = 2*(AC.Wing2.eta_outboard-AC.Wing2.eta_inboard)*AC.Wing2.WingSpan/2;
W_tef = AC.Wing2.Swf*2.706*kf*(AC.Wing2.Swf*bfs)^(3/16)*...
    ((1.8*ME.Landing.Speed/100)^2*sin(delta_f_L*pi/180)/AF.t_c)^(3/4);
%Spoilers:
W_sp = 0.015*W_w_basic;

% Total Weight:
WE.Wing2 = W_w_basic + 1.2*(W_tef+W_sp);


F = WE.Wing2-x;
end