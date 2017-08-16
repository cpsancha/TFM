%% Based on  Floating Device Optimized Testing Source Code by Alan Canamar %%

rhos = 1025; %Average Density of Salt Water [kg/m^3]
BouyR = 1.9; %Bouyancy Reserve [Percent]
HullR = 1.5; %Hull Displacement [Percent]
OutR = BouyR - HullR; %Outrigger Displacement [Percent]
delT = BouyR*AC.Weight.MTOW; %Aircraft Total Displacement on Water [kg]
delb = HullR*AC.Weight.MTOW; %Boat Hull Total Displacement [kg]
del0 = OutR*AC.Weight.MTOW; %Outriggers Total Displacement [kg]
del1 = del0/2; %Displacement on one Outrigger [kg]
V = delb./rhos; %Displacement Volume [m^3] Volumen total de hull necesara

K =.08; %0.0675; %Satisfactory Spray
AC.Hull.Length_Beam  =  mean(loadFields(SP,'Hull.Length')./(loadFields(SP,'Hull.Beam')), 'omitnan');
AC.Hull.Lf_Beam =  mean(loadFields(SP,'Hull.Lf')./(loadFields(SP,'Hull.Beam')), 'omitnan');
 %ya definido a traves de semajantes
lf_lah = AC.Hull.Length_Beam/AC.Hull.Lf_Beam;

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
Vha = (0.45)*AC.Hull.Length*AC.Hull.Beam*AC.Hull.Height; %Actual Volume [m^3]


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
%%

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
Ka = 1.3;
KGh = Ka*(h/2); %Center of Gravity from keel [m]
K1 = 0.036; %Proportionality Coefficient
Ih = K1*Lh*b^3; %Moment of Inertia [m^4]
BMh = (Ih/V); %Distance from CB to Metacentre [m]
BGh = KGh - KBh; %Distance from CG to CB [m]
GMh = BMh - BGh; %Metacentric Height [m]

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

%%% Righting Moment for Hull
RMh = delH*sin(theta*(pi/180))*GMh;



%% Meter en script de la polar
%% Flat plate drag area of the Boat Hull
%Boat Hull Geometry
%b stands for Boat Hull
rb = (b/2); %Radius of Boat Hull [m]
KA = 0.7; %Proportionality Coefficient
AHull = KA*Lh*b; %Area of Load Water Plane of Hull [m^2]
Swetb = 0.5*((pi*rb^2)+AHull+(pi*rb*Lh));%Wetted area of Boat Hull [m^2]
Qb = 1.25; %Interference Factor
Reb = (V*Lh*rho)/v; %Reynolds Number
Cfb = 0.455./(log10(Reb)).^2.58;%Friction coefficient
Amaxb = (pi*(rb/2)^2)/4; %Boat Hull Cross Area [m^2]
ldb = Lh/sqrt((4/pi)*Amaxb); %Fineness ratio
Fb = 1+(60/(ldb^3))+(ldb/400); %Boat Hull Form Factor
fb = Cfb.*Fb.*Qb.*Swetb; %Flat plate drag area [m^2]




