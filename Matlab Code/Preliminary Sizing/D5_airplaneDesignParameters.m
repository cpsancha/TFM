%********************************************************************************************************
%*                                                                                                      *
%*                                  D - AIRPLANE DESIGN PARAMETERS:                                     *
%*                                                                                                      *
%*    In this script are estimated the wing area (Sw), wing aspect ratio (A), take-off thrust (T_TO),   *
%*    and maximum lift coefficients (CL_max, CL_max_TO and CL_max_L), acording to the previously        *
%*    obtained values of the MTOW and EW.                                                               *
%*                                                                                                      *
%********************************************************************************************************

Wto_S= linspace(10,1000,100);

%% POLAR ESTIMATION
%Values needed from similar airplanes
W_S = mean(loadFields(SP,'Wing.WingLoading'),'omitnan'); %Wing loading from similar airplanes
A   = mean(loadFields(SP,'Wing.AspectRatio'),'omitnan'); %Aspect ratio from similar airplanes

%Eq 3.22
S_wet = CF.ft2m^2*10^( Parameters.Table_3_5.c + Parameters.Table_3_5.d * log10(AC.Weight.MTOW/CF.lbm2kg));
%Eq 3.21
f = CF.ft2m^2*10^(Parameters.Table_3_4.a + Parameters.Table_3_4.b * log10(S_wet/(CF.ft2m^2)));

S =  AC.Weight.MTOW/ W_S;
C_D0 = f/S; 

%Compresibility corrections: torenbeek p 172 pdf
%deltaC_D0 =0.0005 long-range cruise conditions
%deltaC_D0 =0.002 high-speed cruise conditions

k=1/(pi*A*Parameters.Table_3_6.e.clean(1));


%% 3.1 Sizing to stall speed requirements

%% 3.2 Sizing to take-off distance requirements
[~, ~, ~, rho] = atmosisa(ME.TakeOff.Altitude);
[~, ~, ~, rho0] = atmosisa(0);
sigma = rho/rho0;

S_tofl = ME.TakeOff.S_TOFL/CF.ft2m; %stofl in ft

TOP_25 = (S_tofl)/37.5; %TOP in lbs/ft^2 Eq 3.8

Cl_max_TO = 1.2;

Wto_S_to=Wto_S./(CF.lbf2N/CF.ft2m^2); %W/S to lbs/ft^2

Tto_Wto_TO = Wto_S_to./ (Cl_max_TO.*TOP_25*sigma); 

% figure()
% plot(Wto_S,Tto_Wto_TO)
% 
% hold all



%% 3.3 Sizing to landing distance requirements

%% 3.4 Sizing to climb requirementes **
Wto_S_cl = Wto_S./(CF.lbf2N/CF.ft2m^2);  %W/S to lbs/ft^2

% Wto_S_cl = [20, 30, 40, 50];

[~, ~, ~, rho] = atmosisa(ME.TakeOff.Altitude); 
[~, ~, ~, rho0] = atmosisa(0);
sigma = rho/rho0;

%%%%%%%%%%%%%%%%%
% FAR 23.65 AOE
RC = 300; %FPM
RCP = 300/33000; % hp/lbs

  % Polar corrections: Take-off flaps
e_cl = Parameters.Table_3_6.e.take_off_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end);


CL_RC_max = (3*C_D0_cl*pi*A*e_cl)^0.5;
CD_RC_max = 4*C_D0_cl;
CL32_CD = 12.1; %(CL_RC_max^(3/2)/CD_RC_max); %!!!!!!!!!!!!!!!!!!!!!!!!

% (20^0.5/(19*12.1)+0.0091)/0.8
Parameters.Cruise.n_p=0.8; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

W_P.cl.AOE = (((RCP + (Wto_S_cl).^(0.5)./(19.*CL32_CD.*sigma^0.5))./Parameters.Cruise.n_p).^-1)./1.1; % Eq 3.23 hp/lbs, Pto/Pmaxcont=1.1
P_W.cl.FAR2365 = (W_P.cl.AOE*CF.hp2watts/(CF.lbm2kg*CST.GravitySI)).^-1;
%%%%%%%%%%%%%%%%
% FAR 23.67 OEI
[~, ~, ~, rho5000] = atmosisa(5000*CF.ft2m);
sigma = rho5000/rho0;
rho5000 = rho5000/CF.slug2kq*(CF.ft2m)^3;


V_s0 = CF.fps2kts*(2.*Wto_S_cl./(rho5000*1.7)).^0.5;

RC = 0.027.*V_s0.^2;
RCP = RC./33000;

% Polar corrections: flaps most favorable, stopped propeller
e_cl = Parameters.Table_3_6.e.clean(1);
C_D0_cl = C_D0 + 0.005; %increase from stopped propeller

CL_RC_max = (3*C_D0_cl*pi*A*e_cl)^0.5;
CD_RC_max = 4*C_D0_cl;
CL32_CD = 13; %(CL_RC_max^(3/2)/CD_RC_max); !!!!!!!!!!!!!!!!!!!!!!!!!!


W_P.cl.OEI = 0.85*0.5* (((RCP + (Wto_S_cl).^(0.5)./(19.*CL32_CD.*sigma^0.5))./Parameters.Cruise.n_p).^-1);% 0.5-> 1 engine, 0.85-> sealevel
P_W.cl.FAR2367= (W_P.cl.OEI*CF.hp2watts/(CF.lbm2kg*CST.GravitySI)).^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb gradient requirements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% FAR 23.65 CGR
CGR = 1/12;

 %Polar corrections
e_cl = Parameters.Table_3_6.e.take_off_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end);

CL_cl = 1.6; %Parameters.CL_max_TO(1)-0.2;  %Take-off flaps, 0.2 margin for preventing stall p.138, p.132 !!!!!!!!!!!!!!!!

C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;

L_D_cl = 9.6; %CL_cl/C_D0_cl  !!!!!!!!!!!!!!!!!!!

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.FAR2365CGR = (1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.FAR2365CGR = (W_P.cl.FAR2365CGR*CF.hp2watts/(CF.lbm2kg*CST.GravitySI)).^-1;

%%%%%%%%%%%%%%%%
% FAR 23.77
CGR = 1/30;

% Landing flaps + Landing gear out
e_cl = Parameters.Table_3_6.e.landing_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.landing_flaps(end)+ Parameters.Table_3_6.deltaC_D0.landing_gear ;

CL_cl = 1.8; %Parameters.CL_max_L(1)-0.2;  %Take-off flaps, 0.2 margin for preventing stall p.138, p.132 !!!!!!!!!!!!!!!!
C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;
L_D_cl = 6.8; %CL_cl/C_D0_cl  !!!!!!!!!!!!!!!!!!!


CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.FAR2377 = 0.85.*(1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.FAR2377 = (W_P.cl.FAR2377*CF.hp2watts/(CF.lbm2kg*CST.GravitySI)).^-1;

%%%%%%%%%%%
% Climb requirements plot
figure(1)
 hold all

plot(Wto_S,P_W.cl.FAR2365,'DisplayName','23.65')
plot(Wto_S,P_W.cl.FAR2367,'DisplayName','23.67')
plot(Wto_S,P_W.cl.FAR2365CGR,'DisplayName','23.65 CGR')
plot(Wto_S,P_W.cl.FAR2377, 'DisplayName','23.77')
% plot(W_P.cl.OEI,Wto_S_cl)
% plot(W_P.cl.AOE,Wto_S_cl)
% plot(W_P.cl.FAR2365,Wto_S_cl)
% plot(W_P.cl.FAR2377,Wto_S_cl)

legend('show')
%% 3.5 Sizing to maneuvering requierements

%% 3.6 Sizing to cruise speed requirementes **
[T, a, P, rho] = atmosisa(ME.Cruise.Altitude);
q = 0.5*rho*ME.Cruise.Speed^2;

Wcr_WTO = 1;
for i = [1,2,3,4]
Wcr_WTO =  Wcr_WTO*Parameters.fuelFraction(i).value; % Wcrucero entre mtow
end

Pto_Pcr = 0.65; %empuje al despegue versus empuje crucero

Wcr_S = Wcr_WTO * Wto_S;

Tcr_Wcr_cr = (C_D0*q)./Wcr_S + k.*Wcr_S./q;
Tto_Wto_cr = Tcr_Wcr_cr * Wcr_WTO;
% Tto_Wto_cr = ((C_D0*q)./Wto_S + Wto_S.*(Wcr_WTO^2*k)); %?
P_W.cr = Pto_Pcr.*Tto_Wto_cr.*ME.Cruise.Speed./Parameters.Cruise.n_p

plot(Wto_S,P_W.cr)



% 
% figure(2)
% hold all
% plot(Wto_S,(C_D0*q)./Wto_S)
%  plot(Wto_S, Wto_S.*(Wcr_WTO^2*k))


% S_wet = 10^(parameters.Table_3_5.c + parameters.Table_3_5.d * ...
% log10(W_TO))  Usar si las constantes se calculan en SI

