% %% WING CONFIGURATION --> Chapter 6 & 7 Roskam; Chapter 7 Torenbeek
% 
% % Selection of wing configuration and geometric characteristics
% 
%     % Cantilever wing (without braces)
%     % Leading wing --> High[ ] / Low[ ]
%     % Rear wing --> High[ ] / Low[ ]
%     % Zero sweep[ ] / Positive sweep[X] / Negative sweep[ ]
%     % Aspect ratio
%     % Thickness ratio
%     % Airfoils
%     % Taper ratio
%     % Twist
%     % Incidence angle
%     % Dihedral angle
%     % High lift and control surface requirements
%     % Winglets
% 
 %% Airfoil characteristics
%  % Get Reynolds number at each chord location: (para qué?)0
%  [~,~,~,~,nu,~] = atmos(ME.Cruise.Altitude);
%  Reynolds = ME.Cruise.Speed .* c / nu;

 
AF.cl_alpha =1*(180/pi)*(0.991-0.15)/(8-0); %1/rad,  NACA 63A210 pag 557 theory of wing sections
AF.alpha_l0 = -1.5*pi/180;
AF.cl_max = 1.41;
AF.c_mac = -0.25;
AF.t_c   = 0.1;
AF.cli   = 0.005; %cl for minimum drag

%% Wing extra parameters
%Only for Wing1
AC.Wing1.Root_LE = 4;

%Common parameters

AC.Wing1.Torsion =-5*pi/180;  %-0.024
AC.Wing1.TaperRatio = 0.4;

% AC.Wing1.Torsion =-0.1;% -6*pi/180;  %-0.024
% AC.Wing1.TaperRatio = 0.458; % 0.458 %Menor que 0.4, argumentando que el taper ratio "efectivo" es diferente respecto a un ala de simple estrechamiento
AC.Wing1.Dihedral = 0;

DP.Stagger = 10;
DP.VerticalGap = 0;

%%
plotFlag = 0;
options = optimoptions('fsolve','FunctionTolerance',1e-12,...
                           'StepTolerance',1e-9,...
                           'Display','none');
%Get medium cruise weight:
CruiseFF = 1; %Mission Fuel Fraction
for i=1:4
    CruiseFF = CruiseFF*Parameters.fuelFraction(i).value;
end
W = CST.GravitySI*AC.Weight.MTOW*(CruiseFF + CruiseFF*Parameters.fuelFraction(5).value)/2;
AC.Wing.CLdesign = 2*W/(ME.Cruise.q* AC.Wing.Sw);
incidences = fsolve(@(incidences)getWings(incidences,AC, DP, ME, Parameters, AF, plotFlag,W), [0.1,0.1,10],options);   



%% VTP DESIGN
%Volume coefficient = Sv * lv / (Sw * bw)
coefVolume = 0.065; % From similar planes, NOT USED


lv = AC.Fuselage.fusLength - DP.x_cg; %Posicion aproximada
Ye = AC.Engine.Position(2);
Sr_S = 0.25;
AC.VTP.AspectRatio   = 1.9;
AC.VTP.TaperRatio = 0.5; %Este completamente a boleo, sugerencias?
AC.VTP.Sweep_14      = 35; %Grados
Sweep_r = 0; % ángulo de la línea de rotacion del rudder respecto a la horizontal.
kdr = 1.05; %aprox, para deflexion máxima de 30º
kv = 1; %por no ser cola en T

%Fig 24 p. 355
clear x y 
x = Ye/lv * (AC.Engine.Power/CF.hp2watts * Parameters.CL_max_TO) / (AC.Weight.MTOW - ME.Payload) ;
load('VTPsizing.mat')
y = interp1(VTPsizing(:,1),VTPsizing(:,2),x,'linear','extrap');

Sv_S = y / (kdr*kv*(Sr_S*AC.VTP.AspectRatio*cos(Sweep_r*pi/180))^(1/3));

AC.VTP.Sw = Sv_S* AC.Wing.Sw;
AC.VTP.Swet = 2*AC.VTP.Sw;
AC.VTP.WingSpan = sqrt(AC.VTP.Sw * AC.VTP.AspectRatio);
AC.VTP.CMG = AC.VTP.Sw/AC.VTP.WingSpan;
AC.VTP.RootChord = AC.VTP.CMG * 2/(1+AC.VTP.TaperRatio);
AC.VTP.TipChord = AC.VTP.RootChord * AC.VTP.TaperRatio;
AC.VTP.CMA = (2/3)*AC.VTP.RootChord*((1+AC.VTP.TaperRatio+AC.VTP.TaperRatio^2)/(1+AC.VTP.TaperRatio));
AC.VTP.Sweep_12 = atand(tand(AC.VTP.Sweep_14)+(4/AC.VTP.AspectRatio)*((1-AC.VTP.TaperRatio)/(1+AC.VTP.TaperRatio))*(0.25-0.5));

y_vtp = linspace(0, AC.VTP.WingSpan, 100);
c_vtp = AC.VTP.RootChord + (AC.VTP.TipChord-AC.VTP.RootChord)/AC.VTP.WingSpan .*y_vtp;
AC.VTP.x_ac_wf = 0.25* AC.VTP.RootChord + tan(AC.VTP.Sweep_14*pi/180)/AC.VTP.Sw *trapz(y_vtp,c_vtp.*y_vtp);

AC.VTP.Root_LE = DP.x_cg + lv - AC.VTP.x_ac_wf;
AC.VTP.t_c = 0.12; %NACA 0012

% Airfoil coordinates for fancy plotting in getWings
AC.VTP.Airfoil.designation='0012';
AC.VTP.Airfoil.wantFile = 0;
AC.VTP.Airfoil.n=30;
AC.VTP.Airfoil.HalfCosineSpacing=1;
AC.VTP.Airfoil.is_finiteTE=0;
AC.VTP.Airfoil.Points = naca4gen(AC.VTP.Airfoil);

%% Call function again for plotting
plotFlag = 1;
getWings(incidences,AC, DP, ME, Parameters, AF, plotFlag, W);
plotFlag = 0;


run('G11_polarPrediction.m')
run('H11_highLiftDevices.m')

