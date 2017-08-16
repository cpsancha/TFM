function [GMS,GMSy,Dras] = Hydrostatics(b,Lh,h,bo,Lo,ho,y,~,b2,bstab,Lstab,dstab,l,bstabWT,LstabWT,dstabWT,VF1,VFWT1, GW, rhos)


BouyR = 1.9; %Bouyancy Reserve [Percent]
HullR = 1.5; %Hull Displacement [Percent]
OutR = BouyR - HullR; %Outrigger Displacement [Percent]
delT = BouyR*GW; %Trimaran Displacement Weight [kg]
delH = HullR*GW; %Hull Displacement Weight [kg]
delO = OutR*GW; %Outrigger Displacement Weight [kg]
delO1 = delO/2; %One Outrigger Displacement Weight [kg]
delF = VF1*rhos; %Displacement Weight of Stabilizer Float [kg]
delFWT = VFWT1*rhos; %Displacement Weight of WingTip Float [kg]
V = delH/rhos; %Hull Displacement Volume [m^3]
Vo = delO/rhos; %Twin Outrigger Displacement Volume [m^3]
Vo1 = delO1/rhos; %One Outrigger Displacement Volume [m^3]
VT = V+Vo; %Trimaran Displacement Volume [m^3]
VS = delT/rhos; %Seaplane Displacement Volume [m^3]
%% Airfoil Area Calculation
KA = 0.7; %Proportionality Coefficient
AHull = KA*Lh*b; %Area of Load Water Plane of Hull [m^2]
AFloat = KA*Lo*bo; %Area of Load Water Plane Float [m^2]
ASTAB = KA*Lstab*bstab; %Area of Load Water Plane Float Stabilizer[m^2]
AWT = KA*LstabWT*bstabWT;%Area of Load Water Plane WingTip Float [m^2]
%% Stability Calculation
%%%% Draft Level and Center of Bouyancy in the Lateral Direction [x]
%%% Draft Level and Center of Bouyancy for Hull
Draft = 0.95; %Draft line [Percent]
Drah = V/AHull; %Draught [m]
KBh = Draft*Drah; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Outrigger
Drao = Vo1/AFloat; %Draught [m]
KBo = Draft*Drao; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Stabilizer Float
Drastab = VF1/ASTAB; %Draught [m]
KBstab = Draft*Drastab; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Wing Tip Float
DraWT = VFWT1/AWT; %Draught [m]
KBWT = Draft*DraWT; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Twin Outrigger
Draot = Vo/(2*AFloat); %Draught [m]
KBot = Draft*Draot; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Trimaran Outrigger
AT = AHull+(2*AFloat); %Area of Load Water Plane [m^2]
DraT = VT/AT; %Draught [m]
KBT = Draft*DraT; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Trimaran Stabilizer
ATs = AHull+(2*ASTAB); %Area of Load Water Plane [m^2]
DraTs = VT/ATs; %Draught [m]
KBTs = Draft*DraTs; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Trimaran Wing Tip Float
ATW = AHull+(2*AWT); %Area of Load Water Plane [m^2]
DraTW = VT/ATW; %Draught [m]
KBTW = Draft*DraTW; %Center of Bouyancy [m]
%%% Draft Level and Center of Bouyancy for Seaplane
As = AT; %Area of Load Water Plane [m^2]
Dras = VS/As; %Draught [m]
KBS = Draft*Dras; %Center of Bouyancy [m]
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
%%% Metacentric Height of Outrigger
% KGo = 0.37; %Center of Gravity from keel [m]
KGo = Ka*(ho/2); %Center of Gravity from keel [m]
Io = K1*Lo*bo^3; %Moment of Inertia [m^4]
BMo = (Io/Vo1); %Distance from CB to Metacentre [m]
BGo = KGo - KBo; %Distance from CG to CB [m]
GMo = BMo - BGo; %Metacentric Height [m]
%%% Metacentric Height of Stabilizer Float
% KGo = 0.37; %Center of Gravity from keel [m]
KGstab = Ka*(dstab/2); %Center of Gravity from keel [m]
Istab = K1*Lstab*bstab^3;%Moment of Inertia [m^4]
BMstab = (Istab/VF1); %Distance from CB to Metacentre [m]
BGstab = KGstab - KBstab;%Distance from CG to CB [m]
GMstab = BMstab - BGstab;%Metacentric Height [m]
%%% Metacentric Height of Wing Tip Float
% KGo = 0.37; %Center of Gravity from keel [m]
KGWT = Ka*(dstabWT/2); %Center of Gravity from keel [m]
IWT = K1*LstabWT*bstabWT^3;%Moment of Inertia [m^4]
BMWT= (IWT/VFWT1); %Distance from CB to Metacentre [m]
BGWT = KGWT - KBWT; %Distance from CG to CB [m]
GMWT = BMWT - BGWT; %Metacentric Height [m]
%%% Metacentric Height of Twin Outrigger
% KGt = 0.42; %Center of Gravity from keel [m]
KGt = KGo + 0.08; %Center of Gravity from keel [m]
It = 2*(Io+(AFloat*y^2));%Moment of Inertia [m^4]
BMt = (It/Vo); %Distance from CB to Metacentre [m]
BGt = KGt - KBot; %Distance from CG to CB [m]
GMt = BMt - BGt; %Metacentric Height [m]
%%% Metacentric Height of Twin Stabilizer
% KGt = 0.42; %Center of Gravity from keel [m]
KGts = KGstab + 0.08; %Center of Gravity from keel [m]
Its = 2*(Istab+(ASTAB*l^2));%Moment of Inertia [m^4]
BMts = (Its/(VF1*2)); %Distance from CB to Metacentre [m]
BGts = KGts - KBstab; %Distance from CG to CB [m]
GMts = BMts - BGts; %Metacentric Height [m]
%%% Metacentric Height of Twin Wing Tip Float
% KGt = 0.42; %Center of Gravity from keel [m]
KGtW = KGWT + 0.08; %Center of Gravity from keel [m]
ItW = 2*(IWT+(AWT*(b2/2)^2));%Moment of Inertia [m^4]
BMtW = (ItW/(VFWT1*2)); %Distance from CB to Metacentre [m]
BGtW = KGtW - KBWT; %Distance from CG to CB [m]
GMtW = BMtW - BGtW; %Metacentric Height [m]
%%% Metacentric Height of Trimaran Outrigger
% KGT = 0.826; %Center of Gravity from keel [m]
KGT = (KGh/2.15) + KGt; %Center of Gravity from keel [m]
IT = Ih + It; %Moment of Inertia [m^4]
BMT = (IT/VT); %Distance from CB to Metacentre [m]
BGT = KGT - KBT; %Distance from CG to CB [m]
GMT = BMT - BGT; %Metacentric Height [m]
%%% Metacentric Height of Trimaran Stabilizer
% KGT = 0.826; %Center of Gravity from keel [m]
KGTs = (KGh/2.15) + KGts;%Center of Gravity from keel [m]
ITs = Ih + Its; %Moment of Inertia [m^4]
BMTs = (ITs/VT); %Distance from CB to Metacentre [m]
BGTs = KGTs - KBTs; %Distance from CG to CB [m]
GMTs = BMTs - BGTs; %Metacentric Height [m]
%%% Metacentric Height of Trimaran Wing Tip Float
% KGT = 0.826; %Center of Gravity from keel [m]
KGTW = (KGh/2.15) + KGtW;%Center of Gravity from keel [m]
ITW = Ih + ItW; %Moment of Inertia [m^4]
BMTW = (ITW/VT); %Distance from CB to Metacentre [m]
BGTW = KGTW - KBTW; %Distance from CG to CB [m]
GMTW = BMTW - BGTW; %Metacentric Height [m]
%%% Metacentric Height of Seaplane Trimaran
KGS = 1.75; %Center of Gravity from keel [m]
% KGS = CGL + (KGT*1.75); %Center of Gravity from keel [m]
IS = IT; %Moment of Inertia of Trimaran[m^4]
BMS = (IS/VS); %Distance from CB to Metacentre [m]
BGS = KGS - KBS; %Distance from CG to CB [m]
GMS = BMS - BGS; %Metacentric Height [m]
%%% Metacentric Height of Seaplane Stabilizer
KGSs = 1.84; %Center of Gravity from keel [m]
% KGSs = CGL + (KGTs/3.5);%Center of Gravity from keel [m]
ISs = ITs; %Moment of Inertia of Trimaran[m^4]
BMSs = (ISs/VS); %Distance from CB to Metacentre [m]
BGSs = KGSs - KBTs; %Distance from CG to CB [m]
GMSs = BMSs - BGSs; %Metacentric Height [m]
%%% Metacentric Height of Seaplane Wing Tip Float
KGSW = 1.84; %Center of Gravity from keel [m]
% KGSW = CGL + (KGTW/3.5);%Center of Gravity from keel [m]
ISW = ITW; %Moment of Inertia of Trimaran[m^4]
BMSW = (ISW/VS); %Distance from CB to Metacentre [m]
BGSW = KGSW - KBTW; %Distance from CG to CB [m]
GMSW = BMSW - BGSW; %Metacentric Height [m]
%%%% METACENTRIC HEIGHT LONGITUDINAL [y]
%%% Metacentric Height of Hull
Ihy = K1*Lh^3*b; %Moment of Inertia [m^4]
BMhy = (Ihy/V); %Distance from CB to Metacentre [m]
GMhy = BMhy - BGh; %Metacentric Height [m]
%%% Metacentric Height of Outrigger
Ioy = K1*Lo^3*bo; %Moment of Inertia [m^4]
BMoy = (Ioy/Vo1); %Distance from CB to Metacentre [m]
GMoy = BMoy - BGo; %Metacentric Height [m]
%%% Metacentric Height of Twin Outrigger
Ity = 2*(Ioy+(AFloat*y^2));%Moment of Inertia [m^4]
BMty = (Ity/Vo); %Distance from CB to Metacentre [m]
GMty = BMty - BGt; %Metacentric Height [m]
%%% Metacentric Height of Trimaran
ITy = Ihy + Ity; %Moment of Inertia [m^4]
BMTy = (ITy/VT); %Distance from CB to Metacentre [m]
GMTy = BMTy - BGT; %Metacentric Height [m]
%%% Metacentric Height of Seaplane
ISy = ITy; %Moment of Inertia of Trimaran[m^4]
BMSy = (ISy/VS); %Distance from CB to Metacentre [m]
GMSy = BMSy - BGS; %Metacentric Height [m]
%% Air Ministry Minimum Requirements for Stability
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
% Stability of Outrigger
TGMo = 4*(Vo1)^(1/3); %Transverse Metacentric Height
LGMo = 6*(Vo1)^(1/3); %Longitudinal Metacentric Height
if GMo >= TGMo
fprintf('GMo = Pass');
else
fprintf('GMo = Fail');
end
disp(' ')
if GMoy >= LGMo
fprintf('GMoy = Pass');
else
fprintf('GMoy = Fail');
end
disp(' ')
%Stability of Twin Floats
TGMt = 4*(Vo)^(1/3); %Transverse Metacentric Height
LGMt = 6*(Vo)^(1/3); %Longitudinal Metacentric Height
if GMt >= TGMt
fprintf('GMt = Pass');
else
fprintf('GMt = Fail');
end
disp(' ')
if GMty >= LGMt
fprintf('GMty = Pass');
else
fprintf('GMty = Fail');
end
disp(' ')
% Stability of Trimaran
TGMT = 4*(VT)^(1/3); %Transverse Metacentric Height
LGMT = 6*(VT)^(1/3); %Longitudinal Metacentric Height
if GMT >= TGMT
fprintf('GMT = Pass');
else
fprintf('GMT = Fail');
end
disp(' ')
if GMTy >= LGMT
fprintf('GMTy = Pass');
else
fprintf('GMTy = Fail');
end
disp(' ')
% Stability of Seaplane
TGMS = 4*(VS)^(1/3); %Transverse Metacentric Height
LGMS = 6*(VS)^(1/3); %Longitudinal Metacentric Height
if GMS >= TGMS
fprintf('GMS = Pass');
else
fprintf('GMS = Fail');
end
disp(' ')
if GMSy >= LGMS
fprintf('GMSy = Pass');
else
fprintf('GMSy = Fail');
end
%% Righting Moment
theta = 0:90;
% %%% Float Track location
% RM = delT*y*cos(theta*(pi/180));
% %%% Float Track location
% RMF = delT*l*cos(theta*(pi/180));
% %%% Float Track location
% RMFWT = delT*(b2/2)*cos(theta*(pi/180));
% %%% Righting Moment for Hull
% RMh = delH*sin(theta*(pi/180))*GMh;
% %%% Righting Moment for Outrigger
% RMo = delO1*sin(theta*(pi/180))*GMo;
% %%% Righting Moment for Stabilizer Float
% RMstab = delF*sin(theta*(pi/180))*GMstab;
% %%% Righting Moment for Wing Tip Float
% RMWT = delFWT*sin(theta*(pi/180))*GMWT;
% %%% Righting Moment for Seaplane Trimaran
% RMS = delT*sin(theta*(pi/180))*GMS;
% %%% Righting Moment for Seaplane Stabilizer
% RMSs = delT*sin(theta*(pi/180))*GMSs;
% %%% Righting Moment for Seaplane Wing Tip Float
% RMSW = delT*sin(theta*(pi/180))*GMSW;
%Maximum Righting Moment and Max angle of tilting
thetamax = 90; %Max angle of Tilting [deg]
RMSmax = delT*sin(thetamax*(pi/180))*GMS; %Max Rigthing Moment [kg m]
phi = RMSmax./(delT*y);
phimax = acos(phi)*(180/pi); %Max allowed angle of Tilting
%% Stability in Wind
B = y; %Beam between the centerlines of the outer hulls [m]
CE = GMSy; %Height of the center of effort above the CG [m]
SA = 0:50; %Sail Area [m^2]
SF = 9.48*sqrt((0.5*B*GW)./(SA*CE));%Wind Speed [m/s]
%% Plots
figure
plot(theta,RMS,'k',theta,RMSs,theta,RMSW,'Linewidth',2)
xlabel('Angle of Inclination [deg]')
ylabel('Righting Moment [Kg m]')
title('Righting Moment Transverse Stability')
legend('Seaplane Trimaran','Seaplane Stabilizer','Seaplane Wing Tip Float')
end

