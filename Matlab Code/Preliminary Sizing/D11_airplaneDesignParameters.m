%********************************************************************************************************
%*                                                                                                      *
%*                                  D - AIRPLANE DESIGN PARAMETERS:                                     *
%*                                                                                                      *
%*    In this script are estimated the wing area (Sw), wing aspect ratio (A), take-off thrust (T_TO),   *
%*    and maximum lift coefficients (CL_max, CL_max_TO and CL_max_L), acording to the previously        *
%*    obtained values of the MTOW and EW.                                                               *
%*                                                                                                      *
%********************************************************************************************************


%% Design parameters
A = 9;
Parameters.CL_max    = 1.5; %[1.2, 1.8];
Parameters.CL_max_TO = 2;   %[1.6, 2.2];
Parameters.CL_max_L  = 3;   %[1.8, 3.4];

%Definimos un rango de cargas alares al despegue para representar las gráficas
Wto_S= linspace(100,10000,100);

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


%% 3.1 Sizing to stall speed requirements

%% 3.2 Sizing to take-off distance requirements
[~, ~, ~, rho] = atmosisa(ME.TakeOff.Altitude);
[~, ~, ~, rho0] = atmosisa(0);
sigma = rho/rho0;

S_tofl = ME.TakeOff.S_TOFL/CF.ft2m; %stofl in ft
TOP_25 = (S_tofl)/37.5; %TOP in lbs/ft^2 Eq 3.8
Wto_S_to=Wto_S./(CF.lbf2N/CF.ft2m^2); %W/S to lbs/ft^2
Tto_Wto_TO = Wto_S_to./ (Parameters.CL_max_TO*TOP_25*sigma); 

Tto_Wto_TO = Tto_Wto_TO/CF.lbf2N;
Pto2Tto =2.9*CF.lbf2N/CF.hp2watts;

P_W.take_off_Roskam=Tto_Wto_TO./Pto2Tto;
%%%
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

% W_S = [2,6];
% a=2;
% fun =@(T_W)getStow(W_S,T_W,a);
% T_W_to = fsolve(fun,0)





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
MLW_MTOW = mean(loadFields(SP,'Weight.MLW_MTOW'),'omitnan');
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

%% 3.5 Sizing to maneuvering requierements

%% 3.6 Sizing to cruise speed requirementes **
[T, asound, P, rho] = atmosisa(ME.Cruise.Altitude);
q = 0.5*rho*ME.Cruise.Speed^2;

Wcr_WTO = 1;
for i = [1,2,3,4]
Wcr_WTO =  Wcr_WTO*Parameters.fuelFraction(i).value; % Wcrucero entre mtow
end

Pto_Pcr = 0.65; %potencia al despegue versus potencia crucero

Wcr_S = Wcr_WTO * Wto_S;

Tcr_Wcr_cr = (C_D0*q)./Wcr_S + k.*Wcr_S./q;
Tto_Wto_cr = Tcr_Wcr_cr * Wcr_WTO;
% Tto_Wto_cr = ((C_D0*q)./Wto_S + Wto_S.*(Wcr_WTO^2*k)); %?
P_W.cr = Pto_Pcr.*Tto_Wto_cr.*ME.Cruise.Speed./Parameters.Cruise.n_p;

% 
% figure(1)
% hold all
% plot(Wto_S,P_W.cr,'DisplayName','cr')


% Roskam cruise

Wto_S_cr = Wto_S./(CF.lbf2N/CF.ft2m^2);  %W/S to lbs/ft^2
[~, ~, ~, rho] = atmosisa(ME.Cruise.Altitude);
sigma = rho/rho0;
Ip = 1.7;

W_P_roskam = 0.7*(1/(sigma*Ip^3)).* Wto_S_cr;
P_W.crRoskam = (W_P_roskam.*(CF.hp2watts./(CF.lbm2kg.*CST.GravitySI))^(-1)).^-1;
% figure(1)
% plot(Wto_S,P_W.crRoskam,'DisplayName','Cruise Roskam')
% 
% legend('show')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Design point selection

   figure(); hold on;
    LegendStr=cell(0);
    

    
    %Cruise Roskam
            plot(Wto_S,P_W.crRoskam,'Color',Parameters.Colors(1,:));
            LegendStr{end+1} = 'Max Speed Cruise (Roskam)';
            
    %Cruise
            plot(Wto_S,P_W.cr,'Color',Parameters.Colors(2,:));
            LegendStr{end+1} = 'Max Speed Cruise';          
            
    %Take-off
            plot(Wto_S,P_W.take_off,'Color',Parameters.Colors(3,:));
            LegendStr{end+1} = 'Take off length';
            plot(Wto_S,P_W.take_off_Roskam,'Color',Parameters.Colors(4,:));
            LegendStr{end+1} = 'Take off length (Roskam)';
    %Climb
        %Take-Off
        plot(Wto_S,P_W.cl.CS25121tr,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
            LegendStr{end+1} = 'Climb (Transition Segment)';

        plot(Wto_S,P_W.cl.CS25111,'LineWidth',1.25,'Color',Parameters.Colors(7,:));
            LegendStr{end+1} = 'Climb (Initial Segment)';

        plot(Wto_S,P_W.cl.CS25121nd,'LineWidth',1.25,'Color',Parameters.Colors(8,:));
            LegendStr{end+1} = 'Climb (Second Segment)';
            
        plot(Wto_S,P_W.cl.CS25121er,'LineWidth',1.25,'Color',Parameters.Colors(9,:));
            LegendStr{end+1} = 'Climb (En-Route Segment)';
            
        %Landing
        plot(Wto_S,P_W.cl.CS25119,'LineWidth',1.25,'Color',Parameters.Colors(10,:));
            LegendStr{end+1} = 'Climb (Balked Landing)';
            
        plot(Wto_S,P_W.cl.CS25121ba,'LineWidth',1.25,'Color',Parameters.Colors(11,:));
            LegendStr{end+1} = 'Climb (Balked Approach)';
     %Formating
        xlim([250,max(Wto_S)])
        ylim([0,60])
        set(gcf,'Position',[450   200   700   525])
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
        x = 2000; %[N/m^2] WingLoad
        y =  20; %[W/N] Power to Weight ratio at take-off
        plot(x,y,'o');
    end
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
    saveFigure(ME.FiguresFolder,'DesignPoint')
    warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');

    clear choiceFlag choice p LegendStr

    %Store design point into AC structure
    AC.Wing.WingLoading   = x;
    AC.Weight.Tto_MTOW    = y;
    AC.Wing.AspectRatio   = A;
    AC.Wing.Sw            = AC.Weight.MTOW/AC.Wing.WingLoading;
    AC.Engine.TotalThrust = AC.Weight.Tto_MTOW*(AC.Weight.MTOW*CST.GravitySI);