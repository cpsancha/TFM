%% Based on  Floating Device Optimized Testing Source Code by Alan Canamar %%

rhos = 1025; %Average Density of Salt Water [kg/m^3]
BouyR = 1.9; %Bouyancy Reserve [Percent]
HullR = 1.5; %Hull Displacement [Percent]
OutR = BouyR - HullR; %Outrigger Displacement [Percent]
 %Aircraft Total Displacement on Water [kg]
delb = HullR*AC.Weight.MTOW; %Boat Hull Total Displacement [kg]
del0 = OutR*AC.Weight.MTOW; %Outriggers Total Displacement [kg]
del1 = del0/2; %Displacement on one Outrigger [kg]
V = delb./rhos; %Displacement Volume [m^3] Volumen total de hull necesara

K =.085; %0.0675; %Satisfactory Spray
AC.Hull.Length_Beam  = mean(loadFields(SP,'Hull.Length')./(loadFields(SP,'Hull.Beam')), 'omitnan');
AC.Hull.Lf_Beam =  mean(loadFields(SP,'Hull.Lf')./(loadFields(SP,'Hull.Beam')), 'omitnan');
 %ya definido a traves de semajantes
% lf_lah = AC.Hull.Length_Beam/AC.Hull.Lf_Beam;
%lf_lah = 2.07;ac.



%Main dimensions
AC.Hull.Cv0 = K * AC.Hull.Lf_Beam^2;
AC.Hull.Beam = (delb/(rhos*AC.Hull.Cv0))^(1/3);  %Compatible with AC.Fuselage.fusWidth
AC.Hull.Length = AC.Hull.Length_Beam * AC.Hull.Beam;  %Compatible with AC.Fuselage.fusLength
AC.Hull.Lf = AC.Hull.Lf_Beam * AC.Hull.Beam;
AC.Hull.La = AC.Hull.Length-AC.Hull.Lf;

%Other dimensions
Bhh = 0.65; %Bow height increment [Percent]
Shh = 0.09; %Step Height increment [Percent]
chh = 0.08; %Chine Flare increment [Percent]
Lch = 1.7; %Flat Forebody increment [Percent]

AC.Hull.Height = AC.Hull.Beam*Bhh; %Bow Height [m]
S = AC.Hull.Beam*Shh; %Step Height [m]
ch = AC.Hull.Beam*chh; %Chine Flare [m]
Lc = AC.Hull.Beam*Lch; %Flat Forebody Length [m]
AC.Hull.Vha = (0.45)*AC.Hull.Length*AC.Hull.Beam*AC.Hull.Height; %Actual Volume [m^3]
Vha = AC.Hull.Vha;

if Vha >= V
fprintf('Vha > V = Pass');
else
fprintf('Vha < V = Fail');
end
%% 
% [GMS,GMSy,Dras] = Hydrostatics(b,Lh,h,bo,Lo,ho,y,~,b2,bstab,Lstab,dstab,l,bstabWT,LstabWT,dstabWT,VF1,VFWT1, GW, rhos)
% [GMS,GMSy,Dras] = Hydrostatics(AC.Hull.Beam,AC.Hull.Length,AC.Hull.Height,...
%     bo,Lo,ho,y,~,... %Outrigger
%     AC.Wing1.WingSpan,...
%     bstab,Lstab,dstab,l,bstabWT,LstabWT,dstabWT,...
%     VF1,VFWT1, AC.Weight.MTOW, rhos)


% [GMS,GMSy,Dras] = Hydrostatics(AC.Hull.Beam,AC.Hull.Length,AC.Hull.Height,...
%     0,0,0,0,0,... %Outrigger
%     AC.Wing1.WingSpan,...
%     0,0,0,0,0,0,0,...
%     0,0, AC.Weight.MTOW, rhos)
%% Stability 

Lh = AC.Hull.Length;
b = AC.Hull.Beam;
h = AC.Hull.Height;

delH = HullR*AC.Weight.MTOW; %Hull Displacement Weight [kg]
V = delH/rhos; %Hull Displacement Volume [m^3]
KA = 0.7; %Proportionality Coefficient
AHull = KA*Lh*b; %Area of Load Water Plane of Hull [m^2]


%%% Draft Level and Center of Bouyancy for Hull
Draft = 0.95; %Draft line [Percent]
Drah = V/AHull; %Draught [m]
KBh = Draft*Drah; %Center of Bouyancy [m]


%%%% METACENTRIC HEIGHT Transverse [x]
%%% Metacentric Height of Hull
% KGh = 0.85; %Center of Gravity from keel [m]

KGh = DP.y_cg; %Center of Gravity from keel [m]

% Ih basado en lh? depende de la constante pero a saber...
K1 = 0.036; %Proportionality Coefficient
Ih = K1*Lh*b^3; %Moment of Inertia [m^4]
BMh = (Ih/V); %Distance from CB to Metacentre [m]
BGh = KGh - KBh; %Distance from CG to CB [m]
GMh = BMh - BGh %Metacentric Height [m]


%%%% METACENTRIC HEIGHT LONGITUDINAL [y]
%%% Metacentric Height of Hull
Ihy = K1*Lh^3*b; %Moment of Inertia [m^4]
BMhy = (Ihy/V); %Distance from CB to Metacentre [m]
GMhy = BMhy - BGh; %Metacentric Height [m]

%%
fprintf(' Air Ministry Minimum Requirements for Stability')
disp(' ')
%Stability of Hull
TGMh = 4*(V)^(1/3); %Transverse Metacentric Height
LGMh = 6*(V)^(1/3); %Longitudinal Metacentric Height
if GMh >= TGMh
fprintf('GMh = Pass');
else
fprintf('GMh = Fail');
end
disp(' ')
if GMhy >= LGMh
fprintf('GMhy = Pass');
else
fprintf('GMhy = Fail');
end
disp(' ')


theta= 8;
% M_R = 0.75*AC.Weight.MTOW/CF.lbm2kg*(GMh/CF.ft2m + (AC.Weight.MTOW/CF.lbm2kg)^(1/3))*sin(theta*pi/180);
% delStabs = M_R /(AC.Wing1.WingSpan/2)

% Esta formula es una mierda.
k=0.029+0.0005*AC.Wing1.WingSpan
deltaF = AC.Weight.MTOW/(AC.Wing1.WingSpan/2)*(k*AC.Weight.MTOW^(1/3)+BGh*sin(theta*pi/180)) 

% Lo que hay que hacer es seguir las FAR:
%Con el flotador totalmente sumergido, el momento de este debe ser 1.5 el
%momento de las fuerzas desestabilizantes

% %%% Righting Moment for Hull
RMh = abs(delH*sin(theta*(pi/180))*GMh); %kg m
M_R = 1.5*RMh; %norma
delStabs = M_R /(AC.Wing1.WingSpan/2) %Desplazamiento necesario del flotador
%Sabiendo el desplazamiento, podemos dimensionar, todas las dimensiones son
%funcion de beam
C = rhos * KA* 4* Bhh * 0.5; %0.6,por meter otra reduccion de volumen respecto al cubo ideal
bstabWT = (delStabs/C)^(1/3); 
LstabWT = 4* bstabWT;
hstabWT = bstabWT * Bhh;

AC.Hull.bstabWT = bstabWT;
AC.Hull.LstabWT = LstabWT;
AC.Hull.hstabWT = hstabWT;


%% Metacentro nuevo caso:
%Ih
K1 = 0.036; %Proportionality Coefficient
Ih = K1*Lh*b^3;
%ItW
KA = 0.7;
AWT = KA*LstabWT*bstabWT;%Area of Load Water Plane WingTip Float [m^2]
IWT = K1*LstabWT*bstabWT^3;%Moment of Inertia [m^4]
ItW =2*(IWT+(AWT*(AC.Wing1.WingSpan/2)^2));%

ITW = Ih + ItW;


OutR = BouyR - HullR;
delO = OutR*AC.Weight.MTOW;
%%%%% Ver esto
Vo = delO/rhos;
% Vo =  2*KA *bStabs*LstabWT*hstabWT;
%%%%%
VT = V+Vo;
ATW = AHull+(2*AWT); %Area of Load Water Plane [m^2]
DraTW = VT/ATW; %Draught [m]
KBTW = Draft*DraTW; %Center of Bouyancy [m]

delT = BouyR*AC.Weight.MTOW;
VS = delT/rhos;
%%% Metacentric Height of Seaplane Wing Tip Float
KGSW = DP.y_cg; %Center of Gravity from keel [m]
% KGSW = CGL + (KGTW/3.5);%Center of Gravity from keel [m]
ISW = ITW; %Moment of Inertia of Trimaran[m^4]
BMSW = (ISW/VS); %Distance from CB to Metacentre [m]
BGSW = KGSW - KBTW; %Distance from CG to CB [m]
GMSW = BMSW - BGSW %Metacentric Height [m]

%Requisito EEUU
metacentroEEUU = 0.41*AC.Weight.MTOW^(1/3)
metacentroUK = 4*(V)^(1/3)

% %% Water Reaction
% Vl = sqrt(2*CST.GravitySI *  AC.Weight.MLW/( ME.Cruise.Density*AC.Wing.Sw * Parameters.CL_max_L));
% ME.Landing.Speed = Vl;
% Vso = ME.Landing.Speed/CF.kts2ms; 
% C1 = 0.012;
% W = (AC.Weight.EW/CF.lbm2kg);
% AC.Hull.Beta = 25;
% nw_lim = (C1 * Vso^2)/((tan(AC.Hull.Beta*pi/180))^2*W)^(1/3)
% 
% K2 = 0.75;
% C4 =0.078*C1;
% P = C4 * K2 * Vso^2/(tan(AC.Hull.Beta*pi/180));
% Area = AC.Hull.Beam*AC.Hull.Beam;%0.125*AHull;% AC.Hull.Beam*AC.Hull.Beam
% F=P*(Area/CF.in2m^2)
% 
% nw2 = F/W