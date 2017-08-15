%********************************************************************************************************
%*                                                                                                      *
%*                                  D - AIRPLANE DESIGN PARAMETERS:                                     *
%*                                                                                                      *
%*    In this script are estimated the wing area (Sw), wing aspect ratio (A), take-off thrust (T_TO),   *
%*    and maximum lift coefficients (CL_max, CL_max_TO and CL_max_L), acording to the previously        *
%*    obtained values of the MTOW and EW.                                                               *
%*                                                                                                      *
%********************************************************************************************************

clear Parameters.CL_max_L Parameters.CL_max_L
%% Design parameters
% Seleccionar tipo de barrido:
% 1: En alargamiento
% 2: En Cl_max_TO
% 3: EN Cl_max_L

barrido=3;

switch barrido
    case 1
A_vec = [8,10,12];
CL_max_TO_vec = 2.8;
CL_max_L_vec = 3.8;
    case 2
A_vec = 9.1;
CL_max_TO_vec = [1.6, 2.2, 2.8];
CL_max_L_vec = 3.8;
    case 3
% A_vec = 9.1;
% CL_max_TO_vec = 2.8;
% CL_max_L_vec = [2.5, 3.5, 3.8];

A_vec = 9.1;
CL_max_TO_vec = 2.15;%2.405;
CL_max_L_vec = [2.2, 2.6, 2.6];
end

%A=9.1 Clto = 2.8 clL=3.8


line={':','--','-'};
for ii=1:length(A_vec)
    for jj=1:length(CL_max_TO_vec)
        for kk=1:length(CL_max_L_vec)


    A=A_vec(ii);
    Parameters.CL_max_TO = CL_max_TO_vec(jj);
    Parameters.CL_max_L = CL_max_L_vec(kk);
    
% A = 9;
% Parameters.CL_max    = 1.5; %[1.2, 1.8]; %Doesn't affect any restriction!
% Parameters.CL_max_TO = 1.6;   %[1.6, 2.2];
% Parameters.CL_max_L  = 3.4;   %[1.8, 3.4];

%Definimos un rango de cargas alares al despegue para representar las gráficas
Wto_S= linspace(100,6000,100);

%% Polar estimation
%Values needed from similar airplanes
W_S = mean(loadFields(SP,'Wing.WingLoading'),'omitnan')./(CF.lbf2N/CF.ft2m^2);        %Wing loading from similar airplanes

%Eq 3.22
S_wet = (CF.ft2m^2)*10^(Parameters.Table_3_5.c + Parameters.Table_3_5.d*log10(AC.Weight.MTOW*CST.GravitySI*CF.N2lbf)); %[m^2]
%Eq 3.21
f = (CF.ft2m^2)*10^(Parameters.Table_3_4.a + Parameters.Table_3_4.b * log10(S_wet*(CF.m2ft^2))); %[m^2]
S =  AC.Weight.MTOW/ W_S;

C_D0 = f/S; 
k=1/(pi*A*Parameters.Table_3_6.e.clean(1));

%% Thrust model at Take-off and power to thrust relationship
Pto2Tto =2.9*CF.lbf2N/CF.hp2watts; %W to N

% T= T_to+ T_to*m*Vr, m=(0.63-1)/(100*CF.kts2ms)

%% 3.1 Sizing to stall speed requirements

%% 3.2 Sizing to take-off distance requirements


[~, ~, ~, rho] = atmosisa(ME.TakeOff.Altitude);
[~, ~, ~, rho0] = atmosisa(0);
sigma = rho/rho0;

S_tofl = 1.3*ME.TakeOff.S_TOFL/CF.ft2m; %stofl in ft
TOP_25 = (S_tofl)/37.5; %TOP in lbs/ft^2 Eq 3.8
Wto_S_to=Wto_S./(CF.lbf2N/CF.ft2m^2); %W/S to lbs/ft^2
Tto_Wto_TO = Wto_S_to./ (Parameters.CL_max_TO*TOP_25*sigma); 

P_W.take_off =Tto_Wto_TO./2.9; %HP/LBF
P_W.take_off1=P_W.take_off.*CF.hp2watts./CF.lbf2N;

Tto_Wto_TO = Tto_Wto_TO./CF.lbf2N; %lbf/N
P_W.take_off_Roskam=Tto_Wto_TO./Pto2Tto;

% figure(13)
% hold on
% plot(Wto_S, P_W.take_off_Roskam,'^')
% plot(Wto_S, P_W.take_off1,'*')
% %%%
% Take off distance method from paper 3. Design Point
%%%%
a= [12.54183, -6.77017, 0.0827, -0.90283, -0.10521, -1.42082, -3.73432, 0.28393];

%Preliminary value for the AC beam:
AC.Hull.Beam =  mean(loadFields(SP,'Hull.Beam'),'omitnan'); %Set the beam of the hull from similar planes

Cdelta0=mean(loadFields(SP,'Weight.MTOW')./(CST.WaterDensity*(loadFields(SP,'Hull.Beam')).^3), 'omitnan');
Lf_B = mean(loadFields(SP,'Hull.Lf')./(loadFields(SP,'Hull.Beam')), 'omitnan');
Beta = mean(loadFields(SP,'Hull.Beta') ,'omitnan');
C_D0_to= C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end);

%Solve a*x^2+b*x+c=0
c = -( Wto_S./(ME.TakeOff.S_TOFL)-a(7)*rho*C_D0_to-a(8) )./(rho*Parameters.CL_max_TO)...
     +a(3)*Lf_B+a(4)/cos(Beta*pi/180)+a(5)*Cdelta0+a(6);
 
T_W_to = (-a(1)+sqrt((a(1))^2-4*a(2).*c))./(2*a(2));



% plot(Wto_S,sqrt((a(1))^2-4*a(2).*c))

P_W.take_off=T_W_to./Pto2Tto;

clear a c T_W_to
% W_S = [2,6];
% a=2;
% fun =@(T_W)getStow(W_S,T_W,a);
% T_W_to = fsolve(fun,0)





%% 3.3 Sizing to landing distance requirements
[~, ~, ~, rho] = atmosisa(ME.Landing.Altitude); 

%Roskam Method
    
    Vapp_R   = sqrt(1.1*ME.Landing.S_LFL*CF.m2ft/0.3); %[kts] --> Roskam eq 3.16
    VStall_L = Vapp_R/1.3;               %[kts] --> Roskam eq 3.15
    WL_Sw    = (rho*Parameters.CL_max_L*(VStall_L*CF.kts2ms)^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingRoskam = WL_Sw/MLW_MTOW*CST.GravitySI;

    clear Vapp_R VStall_L WL_Sw

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
CL32_CD = (CL_RC_max^(3/2)/CD_RC_max); 
Parameters.Cruise.n_p=Parameters.Cruise.n_p(1); %!!!!!!!!!!!!!!!!

W_P.cl.AOE = (((RCP + (Wto_S_cl).^(0.5)./(19.*CL32_CD.*sigma^0.5))./Parameters.Cruise.n_p).^-1)./1.1; % Eq 3.23 hp/lbs, Pto/Pmaxcont=1.1
P_W.cl.FAR2365 = (W_P.cl.AOE*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;
%%%%%%%%%%%%%%%%
% FAR 23.67 OEI
[~, ~, ~, rho5000] = atmosisa(5000*CF.ft2m);
sigma = rho5000/rho0;
rho5000 = rho5000/CF.slug2kg*(CF.ft2m)^3;


V_s0 = CF.fps2kts*(2.*Wto_S_cl./(rho5000*1.7)).^0.5;

RC = 0.027.*V_s0.^2;
RCP = RC./33000;

% Polar corrections: flaps most favorable, stopped propeller
e_cl = Parameters.Table_3_6.e.clean(1);
C_D0_cl = C_D0 + 0.005; %increase from stopped propeller

CL_RC_max = (3*C_D0_cl*pi*A*e_cl)^0.5;
CD_RC_max = 4*C_D0_cl;
CL32_CD = (CL_RC_max^(3/2)/CD_RC_max);


W_P.cl.OEI = 0.85*0.5* (((RCP + (Wto_S_cl).^(0.5)./(19.*CL32_CD.*sigma^0.5))./Parameters.Cruise.n_p).^-1);% 0.5-> 1 engine, 0.85-> sealevel
P_W.cl.FAR2367= (W_P.cl.OEI*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Climb gradient requirements %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% FAR 23.65 CGR
CGR = 1/12;

 %Polar corrections
e_cl = Parameters.Table_3_6.e.take_off_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end);

CL_cl = Parameters.CL_max_TO(1)-0.2;  %Take-off flaps, 0.2 margin for preventing stall p.138, p.132 

C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;

L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.FAR2365CGR = (1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.FAR2365CGR = (W_P.cl.FAR2365CGR*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

%%%%%%%%%%%%%%%%
% FAR 23.77
CGR = 1/30;

% Landing flaps + Landing gear out
e_cl = Parameters.Table_3_6.e.landing_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.landing_flaps(end)+ Parameters.Table_3_6.deltaC_D0.landing_gear ;

CL_cl = Parameters.CL_max_L(1)-0.2;  %Take-off flaps, 0.2 margin for preventing stall p.138, p.132 
C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;
L_D_cl = CL_cl/C_D0_cl;


CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.FAR2377 = 0.85.*(1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.FAR2377 = (W_P.cl.FAR2377*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

%%%%%%%%%%%
% % Climb requirements plot
% figure(1)
%  hold all
% 
% plot(Wto_S,P_W.cl.FAR2365,'DisplayName','23.65')
% plot(Wto_S,P_W.cl.FAR2367,'DisplayName','23.67')
% plot(Wto_S,P_W.cl.FAR2365CGR,'DisplayName','23.65 CGR')
% plot(Wto_S,P_W.cl.FAR2377, 'DisplayName','23.77')
% plot(W_P.cl.OEI,Wto_S_cl)
% plot(W_P.cl.AOE,Wto_S_cl)
% plot(W_P.cl.FAR2365,Wto_S_cl)
% plot(W_P.cl.FAR2377,Wto_S_cl)

%% CS25
    % CS 25.121 (OEI) --> Transition segment climb requirement
        %Climb Gradient: Positive% for 2 engines, 0.3% for 3 engines, 0.5% for 4 engines
        %Flaps: Take-Off
        %Landing Gear: Down
        %Speed: Between Vliftoff=1.1VstallTO and V2=1.13Vstall_TO --> Vliftoff
        %Ground effect: Active
        %Weight: Maximum Take-Off Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
switch ME.Powerplant.Number
            case 2
                CGR = 0.000;
            case 3
                CGR = 0.003;
            case 4
                CGR = 0.005;
end

 %Polar corrections
e_cl = Parameters.Table_3_6.e.take_off_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end)+Parameters.Table_3_6.deltaC_D0.landing_gear;

CL_cl = Parameters.CL_max_TO(1)/(1.1^2); 

C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;

L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.CS25121tr = (1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.CS25121tr = (ME.Powerplant.Number/(ME.Powerplant.Number-1)).*(W_P.cl.CS25121tr*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

    % CS 25.111 (OEI) --> Initial climb segment requirement
        %Climb Gradient: 1.2% for 2 engines, 1.5% for 3 engines, 1.7% for 4 engines
        %Flaps: Take-Off
        %Landing Gear: Up
        %Speed: V2=1.13Vstall_TO
        %Ground effect: Active
        %Weight: Maximum Take-Off Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch ME.Powerplant.Number
            case 2
                CGR = 0.012;
            case 3
                CGR = 0.015;
            case 4
                CGR = 0.017;
        end
        
         %Polar corrections
e_cl = Parameters.Table_3_6.e.take_off_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end);

CL_cl = Parameters.CL_max_TO(1)/(1.13^2); 

C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;

L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.CS25111= (1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.CS25111 = (ME.Powerplant.Number/(ME.Powerplant.Number-1)).*(W_P.cl.CS25111*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;
       
    % CS 25.121 (OEI) --> Second segment climb requirement
        %Climb Gradient: 2.4% for 2 engines, 2.7% for 3 engines, 3.0% for 4 engines
        %Flaps: Take-Off
        %Landing Gear: Up
        %Speed: V2=1.13Vstall_TO
        %Ground effect: Not active
        %Weight: Maximum Take-Off Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch ME.Powerplant.Number
            case 2
                CGR = 0.024;
            case 3
                CGR = 0.027;
            case 4
                CGR = 0.030;
        end
        
        
e_cl = Parameters.Table_3_6.e.take_off_flaps(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.take_off_flaps(end);

CL_cl = Parameters.CL_max_TO(1)/(1.13^2); 

C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;

L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.CS25121nd = (1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.CS25121nd = (ME.Powerplant.Number/(ME.Powerplant.Number-1)).*(W_P.cl.CS25121nd*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

    % CS 25.121 (OEI) --> En-Route segment climb requirement
        %Climb Gradient: 1.2% for 2 engines, 1.5% for 3 engines, 1.7% for 4 engines
        %Flaps: Retracted
        %Landing Gear: Up
        %Speed: VFTO=1.18Vstall_TO
        %Ground effect: Not active
        %Weight: Maximum Take-Off Weight
        %Thrust: Continuous Thrust
        %Engines: N-1
        switch ME.Powerplant.Number
            case 2
                CGR = 0.012;
            case 3
                CGR = 0.015;
            case 4
                CGR = 0.017;
        end
e_cl = Parameters.Table_3_6.e.clean(1);
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.clean;

CL_cl = Parameters.CL_max_TO(1)/(1.18^2); 
C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;
L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.CS25121er = (1./(Wto_S_cl).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.CS25121er = (ME.Powerplant.Number/(ME.Powerplant.Number-1)).*(W_P.cl.CS25121er*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

%LANDING CLIMB REQUIREMENTS
    % CS 25.119 (AEO) --> Balked landing requirements
        %Climb Gradient: 3.2% for 2, 3 and 4 engines
        %Flaps: Landing
        %Landing Gear: Down
        %Speed: VL=1.3Vstall_L
        %Ground effect: Not active
        %Weight: Maximum Landing Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch ME.Powerplant.Number
            case 2
                CGR = 0.032;
            case 3
                CGR = 0.032;
            case 4
                CGR = 0.032;
        end
e_cl = Parameters.Table_3_6.e.landing_flaps;       
C_D0_cl = C_D0 + Parameters.Table_3_6.deltaC_D0.landing_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear;

CL_cl = Parameters.CL_max_L(1)/(1.3^2); 
C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;
L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.CS25119 = (1./(Wto_S_cl.*MLW_MTOW).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.CS25119 = (ME.Powerplant.Number/(ME.Powerplant.Number-1)).*(W_P.cl.CS25119*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

    % CS 25.121 (OEI) --> Balked Approach requirements
        %Climb Gradient: 2.1% for 2 engines, 2.4% for 3 engines, 2.7% for 4 engines
        %Flaps: Approach
        %Landing Gear: Up
        %Speed: VstallAprox<=1.1Vstall_L && V<=1.4VstallAprox
        %Ground effect: Not active
        %Weight: Maximum Landing Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch ME.Powerplant.Number
            case 2
                CGR = 0.021;
            case 3
                CGR = 0.024;
            case 4
                CGR = 0.027;
        end
        
e_cl = (Parameters.Table_3_6.e.landing_flaps+Parameters.Table_3_6.e.take_off_flaps)/2;       
C_D0_cl = C_D0 + (Parameters.Table_3_6.deltaC_D0.landing_flaps+Parameters.Table_3_6.deltaC_D0.landing_flaps)/2+Parameters.Table_3_6.deltaC_D0.landing_gear;

CL_cl = Parameters.CL_max_L(1)/((1.1*1.4)^2); 
C_D0_cl = C_D0_cl + 1/(pi*A*e_cl)*CL_cl^2;
L_D_cl = CL_cl/C_D0_cl ;

CGRP = (CGR + 1/L_D_cl)/CL_cl^0.5;

W_P.cl.CS25121ba = (1./(Wto_S_cl.*MLW_MTOW).^(0.5)).*18.97.*Parameters.Cruise.n_p./CGRP;
P_W.cl.CS25121ba = (ME.Powerplant.Number/(ME.Powerplant.Number-1)).*(W_P.cl.CS25121ba*(CF.hp2watts/(CF.lbm2kg*CST.GravitySI))^(-1)).^-1;

clear CGR CGRP 
%% 3.5 Sizing to maneuvering requierements

%% 3.6 Sizing to cruise speed requirementes **
[T, asound, P, rho] = atmosisa(ME.Cruise.Altitude);
sigma = rho/rho0;
q = 0.5*rho*ME.Cruise.Speed^2;

Wcr_WTO = 1;
for i = [1,2,3,4]
Wcr_WTO =  Wcr_WTO*Parameters.fuelFraction(i).value; % Wcrucero entre mtow
end

Pto_Pcr = 1/0.65; %potencia al despegue versus potencia crucero
Wcr_S = Wcr_WTO * Wto_S;
Tcr_Wcr_cr = (C_D0*q)./Wcr_S + k.*Wcr_S./q;
Tto_Wto_cr = Tcr_Wcr_cr * Wcr_WTO;
% Tto_Wto_cr = ((C_D0*q)./Wto_S + Wto_S.*(Wcr_WTO^2*k)); %?
P_W.cr = Pto_Pcr.*Tto_Wto_cr.*ME.Cruise.Speed./Parameters.Cruise.n_p;

% Roskam cruise
Ip=((loadFields(SP,'Wing.WingLoading')./CF.psf2Pa)./...
    (sigma./(loadFields(SP,'Weight.Pto_MTOW').*CF.lbf2N./CF.hp2watts))).^(1/3);
V = loadFields(SP,'Actuations.Vcruise')./CF.mph2ms;
m = mean(V./Ip, 'omitnan');

% 
% IP=polyval([Ip./V,0],10)
% Ip1=polyval([k,0],ME.Cruise.Speed)
% figure(12)
% hold on
% for i =1:length(V)
% plot(Ip(i),V(i),'x')
% end
% x=linspace(0,5,50);
% plot(x,m.*x)
% Ip1 =ME.Cruise.Speed/CF.mph2ms/m;
% Ip=Ip1;
% Wto_S_cr = Wto_S./(CF.lbf2N/CF.ft2m^2);  %W/S to lbs/ft^2
% % Ip = 1.7;
% W_P_roskam = 0.7*(1/(sigma*Ip^3)).* Wto_S_cr;
% P_W.crRoskam = (W_P_roskam.*(CF.hp2watts./CF.lbf2N).^(-1)).^-1;
% figure(1)
% plot(Wto_S,P_W.crRoskam,'DisplayName','Cruise Roskam')
% 
% legend('show')

clear m V Ip Ip1 Wto_S_cr x Wcr_S Tcr_Wcr_cr Tto_Wto_cr

%% Sizing to gust requirements
WingLoading.Gust=fsolve(@(x)gustWingLoading(x,A, CST, AC, CF, ME),1000);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design point selection

 figure(99); hold on;
    LegendStr=cell(0);

    grid on
    

    %Cruise Roskam
%             plot(Wto_S,P_W.crRoskam,'Color',Parameters.Colors(1,:));
%             LegendStr{end+1} = 'Max Speed Cruise (Roskam)';
    %Gust
            plot(WingLoading.Gust.*ones(1,100),linspace(0,100,100),'LineWidth',1.25,'Color',Parameters.Colors(12,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Gust';
    
    %Cruise
            plot(Wto_S,P_W.cr,'Color',Parameters.Colors(2,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Max Speed Cruise';          
            
    %Take-off
            plot(Wto_S,P_W.take_off,'Color',Parameters.Colors(3,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Water Take off length';
            plot(Wto_S,P_W.take_off1,'Color',Parameters.Colors(4,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Ground Take off length';
            
    %Landing
            plot(WingLoading.LandingRoskam.*ones(1,100),linspace(0,100,100),'LineWidth',1.25,'Color',Parameters.Colors(5,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Landing (Roskam)';
    %Climb
        %Take-Off
        plot(Wto_S,P_W.cl.CS25121tr,'LineWidth',1.25,'Color',Parameters.Colors(6,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Climb (Transition Segment)';

        plot(Wto_S,P_W.cl.CS25111,'LineWidth',1.25,'Color',Parameters.Colors(7,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Climb (Initial Segment)';

        plot(Wto_S,P_W.cl.CS25121nd,'LineWidth',1.25,'Color',Parameters.Colors(8,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Climb (Second Segment)';
            
        plot(Wto_S,P_W.cl.CS25121er,'LineWidth',1.25,'Color',Parameters.Colors(20,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Climb (En-Route Segment)';
            
        %Landing
        plot(Wto_S,P_W.cl.CS25119,'LineWidth',1.25,'Color',Parameters.Colors(15,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Climb (Balked Landing)';
            
        plot(Wto_S,P_W.cl.CS25121ba,'LineWidth',1.25,'Color',Parameters.Colors(18,:),'LineStyle',line{ii+jj+kk-2});
            LegendStr{end+1} = 'Climb (Balked Approach)';
            
        %Similar planes design points
    for i=1:length(SP)
        plot(SP{i}.Wing.WingLoading,SP{i}.Weight.Pto_MTOW,'x')
    end
    
     %Formating
        xlim([250,max(Wto_S)])
        ylim([0,40])
        set(gcf,'Position',[450   200   700   525])
        
        title('Design Point') 
        xlabel('Wing Loading - MTOW/Sw [N/m^2]')
        ylabel('Power/Weight_{TO} [W/N]')
        grafWidth   = 16;
        grafAR      = 0.6;
        set(gcf,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
        warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
        [~,objs]=columnlegend(2,LegendStr,'Location','northwest','FontSize',8,'boxoff');
        drawnow; % make sure everything is finished rendering 
        set(findall(objs, 'type', 'text'), 'fontsize', 8, 'interpreter', 'tex')
        warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');
        clear grafWidth grafAR showRoskamRequirements LegendStr objs
    
    if DP.selectDesignPoint == true
    %Design Point Selector
    choiceFlag = 0;
    uiwait(msgbox('Por favor, selecione de forma gráfica el punto de diseño deseado.','Selección del punto de diseño','modal'));
    while choiceFlag<1
        if choiceFlag==0
            %First, chose graphically
            [x,y] = ginput(1);
            p=plot(x,y,'o');
        end
        % Construct a questdlg with three options
        choice = questdlg(['Has selecionado una carga alar al despegue de ',num2str(x),' (kg/m^2) y un cociente de empuje/peso de ',num2str(y),...
                           '. ¿Es esto correcto, o desea modificar el valor?'],'Selección del punto de diseño','Si, es correcto',...
                           'Elegir otro punto gráficamente','Elegir otro punto numéricamente','Si, es correcto');
        % Handle response
        switch choice
            case 'Si, es correcto'
                choiceFlag = 1;
            case 'Elegir otro punto gráficamente'
                delete(p);
                choiceFlag = 0;
            case 'Elegir otro punto numéricamente'
                delete(p);
                prompt = {'Carga alar [kg/m^2]:','Cociente Empuje/Peso [-]:'};
                dlg_title = '';
                num_lines = 1;
                defaultans = {'',''};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                x = str2double(answer{1});
                y = str2double(answer{2});
                p=plot(x,y,'o');
                clear prompt dlg_title num_lines defaultans answer
                choiceFlag = -1;
        end
    end
    else
        x = 3378; %[N/m^2] WingLoad
        y =  23.9988; % 21.2[W/N] Power to Weight ratio at take-off
%         y =  21;
        plot(x,y,'ko');
    end
    
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
    saveFigure(ME.FiguresFolder,'DesignPoint')
    warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');

    clear choiceFlag choice p LegendStr

    %Store design point into AC structure
    AC.Wing.WingLoading   = x;
    AC.Weight.Pto_MTOW    = y;
    AC.Wing.AspectRatio   = A;
    AC.Wing.Sw            = AC.Weight.MTOW*CST.GravitySI/AC.Wing.WingLoading;
    AC.Engine.TotalThrust = AC.Weight.Tto_MTOW*(AC.Weight.MTOW*CST.GravitySI);
    AC.Engine.TotalPower  = AC.Weight.Pto_MTOW*AC.Weight.MTOW*CST.GravitySI;
%     Parameters.CL_max_L   = 
%     Parameters.CL_max_TO  = 
    
    clear A x y rho rho5000 rho0 sigma
    
        end
    end
end
%% Auxiliary Functions
    
function [F] = gustWingLoading(x,A, CST, AC, CF, ME)

switch ME.MissionType
    case 5
    case 11
[~, asound, P, rho] = atmosisa(ME.Cruise.Altitude);
[~, ~, ~, rho0] = atmosisa(0);

% M = ME.Cruise.Speed/asound
% beta = sqrt(1-M^2)
% k = 4*pi/beta/(2*pi)
% flecha=0; %RAD
% a = 2*pi*A/(2+sqrt((A^2*beta^2/k^2)*(1+(tan(flecha))^2/beta^2)+4)) %DA MUY ALTO


a = 5; %Raymer p.311 High AR unswept
c =sqrt((AC.Weight.MTOW*CST.GravitySI/x)/A);


mu = 2*x/(rho*a*c*CST.GravitySI);
kg = 0.88*mu/(5.3+mu);

%Load Factor for limit maneuvering
    n = 2.1+(24e3/(10e3+AC.Weight.MTOW*CF.kg2lbm));
    if n<2.5
        n = 2.5;
    elseif n>3.8
        n = 3.8;
    end

%Reference Gust Speed
    if ME.Cruise.Altitude*CF.m2ft<20000
    U_EAS = 50*CF.ft2m;
    elseif ME.Cruise.Altitude*CF.m2ft>=20000
    U_EAS = ((50+ (ME.Cruise.Altitude*CF.m2ft-20000)*(25-50)/(50000-20000)))*CF.ft2m;
    end

%Equivalent Cruise Airspeed    
V_EAS =  correctairspeed(ME.Cruise.Speed, asound, P, 'TAS', 'EAS');

F= kg*rho0*U_EAS*V_EAS*a/(2*x)-(n-1);
% W_S_Gust = rho0*a*U*V_EAS/(2*(n-1));  
end
end