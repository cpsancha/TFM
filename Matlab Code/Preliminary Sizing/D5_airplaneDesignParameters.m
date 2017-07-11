%********************************************************************************************************
%*                                                                                                      *
%*                                  D - AIRPLANE DESIGN PARAMETERS:                                     *
%*                                                                                                      *
%*    In this script are estimated the wing area (Sw), wing aspect ratio (A), take-off thrust (T_TO),   *
%*    and maximum lift coefficients (CL_max, CL_max_TO and CL_max_L), acording to the previously        *
%*    obtained values of the MTOW and EW.                                                               *
%*                                                                                                      *
%********************************************************************************************************

%% Definimos un rango de ratios Thrust/Weight al despegue para representar las gráficas --> eje y
    ThrustWeight_TO = linspace(0,1,100);


    
%% DENSITY
    % Calculamos las distintas densidades para despegue, crucero, y nivel de mar
    [~, ~, ~, rho0]   = atmosisa(0);
    [~, ~, ~, rho]    = atmosisa(ME.Cruise.Altitude);
    [~, ~, ~, rho_TO] = atmosisa(ME.TakeOff.Altitude);
    [~, ~, ~, rho_L]  = atmosisa(ME.TakeOff.Altitude); %<--Supongo que es la misma que en despegue


    
%% POLAR ESTIMATION
    %Values needed from similar airplanes
    W_S_temp = mean(loadFields(SP,'Wing.WingLoading'),'omitnan'); %Wing loading from similar airplanes

    %Roskam --> Eq 3.22
    S_wet = (CF.ft2m^2)*10^(Parameters.Table_3_5.c + Parameters.Table_3_5.d*log10(AC.Weight.MTOW*CST.GravitySI*CF.N2lbf)); %[m^2]
    %Roskam --> Eq 3.21
    f = (CF.ft2m^2)*10^(Parameters.Table_3_4.a + Parameters.Table_3_4.b * log10(S_wet*(CF.m2ft^2))); %[m^2]

    Sw =  AC.Weight.MTOW/ W_S_temp;
    C_D0 = f/Sw; 

    %Compresibility corrections: torenbeek p 172 pdf
    %deltaC_D0 =0.0005 long-range cruise conditions
    %deltaC_D0 =0.002 high-speed cruise conditions
    
    %Polar coefficients as function of CL to use with polyval --> CD = CD0 + k*CL^2
    Polar.LowSpeedClean   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0];
    Polar.HighSpeedClean  = [];
    Polar.TakeOffGearUp   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.take_off_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.take_off_flaps];
    Polar.TakeOffGearDown = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.take_off_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.take_off_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear];
    Polar.LandingGearUp   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.landing_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.landing_flaps];
    Polar.LandingGearDown = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.landing_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.landing_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear];
    
    

    
%% 3.1 SIZING TO STALL SPEED REQUIREMENTS
    %Cruise stall speed
    if ~isnan(DP.StallSpeed)
        WingLoading.Stall = ((DP.StallSpeed^2)*rho*DP.CLmax)/(2*CST.GravitySI); %[kg/m^2]
    end
    
    %Landing Stall Speed
    if ~isnan(DP.StallSpeed_L)
        WingLoading.Stall_L = ((DP.StallSpeed_L^2)*rho_L*DP.CLmax_L)/(2*CST.GravitySI); %[kg/m^2]
    end
    
    %Take-Off Stall Speed
    if ~isnan(DP.StallSpeed_TO)
        WingLoading.Stall_TO = ((DP.StallSpeed_TO^2)*rho_TO*DP.CLmax_TO)/(2*CST.GravitySI); %[kg/m^2]
    end



%% 3.2 SIZING TO TAKE-OFF DISTANCE REQUIREMENTS
    sigma = rho_TO/rho0; %Gravity coefficient
%Roskam method
    TOP25R = ME.TakeOff.S_TOFL*CF.m2ft/37.5; %[lbf/ft^2]
    WingLoading.TakeOffRoskam = TOP25R*sigma*DP.CLmax_TO*ThrustWeight_TO*(CF.lbf2N/(CF.ft2m^2))/CST.GravitySI; %[kg/m^2]
    
%Similar planes method
    TOP25_SP = loadFields(SP,'Wing.WingLoading')./(sigma.*loadFields(SP,'Wing.CLmax_TO').*loadFields(SP,'Weight.Tto_MTOW')./CST.GravitySI); %[kg/m^2]
    TOP25=polyval([loadFields(SP,'Actuations.Sto')'\TOP25_SP',0],ME.TakeOff.S_TOFL);
    WingLoading.TakeOffSP = TOP25*sigma*DP.CLmax_TO*ThrustWeight_TO; %[kg/m^2]
    
    
%Obtain static constant from similar planes
    figure()
    plot(loadFields(SP,'Actuations.Sto'),TOP25_SP,'*'); hold on
    plot([min(loadFields(SP,'Actuations.Sto')),max(loadFields(SP,'Actuations.Sto'))],(loadFields(SP,'Actuations.Sto')'\TOP25_SP').*[min(loadFields(SP,'Actuations.Sto')),max(loadFields(SP,'Actuations.Sto'))],'--')
    plot(ME.TakeOff.S_TOFL,TOP25R*CF.lbf2N/CST.GravitySI/(CF.ft2m^2),'+')
    plot(ME.TakeOff.S_TOFL,TOP25,'o')
    xlabel('Take-Off Field Length (TOFL) [m]');
    ylabel('Take-Off Parameter [kg/m^2]');
    h = findobj(gca,'Type','line');
    legend([h(4),h(2),h(1)],{'Similar Planes','Roskam Method','Similar Planes Method'},'Location','southeast')
    legend('boxoff')
    
clear sigma h TOP25 TOP25_SP TOP25R



%% 3.3 SIZING TO LANDING DISTANCE REQUIREMENTS
%Roskam Method
    WL_WTO_R = 0.85;
    Vapp_R   = sqrt(ME.Landing.S_LFL*CF.m2ft/0.3); %[kts] --> Roskam eq 3.16
    VStall_L = Vapp_R/1.3;                 %[kts] --> Roskam eq 3.15
    WL_Sw    = (rho_L*DP.CLmax_L*(VStall_L*CF.kts2ms)^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingRoskam = WL_Sw/WL_WTO_R;

%Similar Planes Method
    WL_WTO_SP = mean(loadFields(SP,'Weight.MLW')./loadFields(SP,'Weight.MTOW')); %From SP: min-->0.7900, max-->0.9267, mean-->0.8753
    Vapp_SP   = sqrt((loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*ME.Landing.S_LFL);
    VStall_L  = Vapp_SP/1.3; %[m/s]
    WL_Sw     = (rho_L*DP.CLmax_L*VStall_L^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingSP = WL_Sw/WL_WTO_SP;



%Obtain static constant from similar planes
    figure()
    plot(loadFields(SP,'Actuations.Sl'),loadFields(SP,'Actuations.Vapprox').^2,'*'); hold on
    plot([450,max(loadFields(SP,'Actuations.Sl'))],(loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*[450,max(loadFields(SP,'Actuations.Sl'))],'--')
    plot(ME.Landing.S_LFL,(Vapp_R*CF.kts2ms).^2,'+')
    plot(ME.Landing.S_LFL,Vapp_SP.^2,'o')
    xlabel('Landing Field Length (LFL) [m]','Interpreter','Tex');
    ylabel('Square of Approach Speed   V_{approx}^2   [m^2/s^2]');
    h = findobj(gca,'Type','line');
    legend([h(4),h(2),h(1)],{'Similar Planes','Roskam Method','Similar Planes Method'},'Location','southeast')
    legend('boxoff')

clear WL_WTO_R WL_WTO_SP Vapp_R Vapp_SP VStall_L WL_Sw h



%% 3.4 SIZING TO CLIMB REQUIREMENTS
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



%% DESIGN POINT
    figure()
    %Take-Off Field Length
        plot(WingLoading.TakeOffRoskam,ThrustWeight_TO); hold on;
        plot(WingLoading.TakeOffSP,ThrustWeight_TO);
    %Take-Off Stall Speed
        if isfield(WingLoading,'Stall_TO')
            plot(WingLoading.Stall_TO.*ones(1,length(ThrustWeight_TO)),ThrustWeight_TO);
        end
    %Landing Field Length
        plot(WingLoading.LandingRoskam.*ones(1,length(ThrustWeight_TO)),ThrustWeight_TO);
        plot(WingLoading.LandingSP.*ones(1,length(ThrustWeight_TO)),ThrustWeight_TO);
    %Formating
        xlabel('Wing Loading - MTOW/Sw [kg/m^2]')
        ylabel('Thrust/Weight_{TO} [-]')
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
    
    clear choiceFlag choice p

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

