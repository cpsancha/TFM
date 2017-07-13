%********************************************************************************************************
%*                                                                                                      *
%*                                  D - AIRPLANE DESIGN PARAMETERS:                                     *
%*                                                                                                      *
%*    In this script are estimated the wing area (Sw), wing aspect ratio (A), take-off thrust (T_TO),   *
%*    and maximum lift coefficients (CL_max, CL_max_TO and CL_max_L), acording to the previously        *
%*    obtained values of the MTOW and EW.                                                               *
%*                                                                                                      *
%********************************************************************************************************

%% PLOT RANGE DEFINITION
    %Definimos un rango de cargas alares Weight/Swing al despegue para representar las gráficas --> eje x
    W_S_TO = linspace(0,450,100);
    %Definimos un rango de ratios Thrust/Weight al despegue para representar las gráficas --> eje y
    T_W_TO = linspace(0,1,100);
    

    
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
    Parameters.Polar.LowSpeedCruise  = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0];
    Parameters.Polar.LongRangeCruise = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0+0.0005];
    Parameters.Polar.HighSpeedCruise = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0+0.0020];
    Parameters.Polar.TakeOffGearUp   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.take_off_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.take_off_flaps];
    Parameters.Polar.TakeOffGearDown = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.take_off_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.take_off_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear];
    Parameters.Polar.LandingGearUp   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.landing_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.landing_flaps];
    Parameters.Polar.LandingGearDown = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.landing_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.landing_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear];
    
    clear f C_D0 W_S_temp
    
    
    
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
    WingLoading.TakeOffRoskam = TOP25R*sigma*DP.CLmax_TO*T_W_TO*(CF.lbf2N/(CF.ft2m^2))/CST.GravitySI; %[kg/m^2]
    
%Similar planes method
    TOP25_SP = loadFields(SP,'Wing.WingLoading')./(sigma.*loadFields(SP,'Wing.CLmax_TO').*loadFields(SP,'Weight.Tto_MTOW')./CST.GravitySI); %[kg/m^2]
    TOP25=polyval([loadFields(SP,'Actuations.Sto')'\TOP25_SP',0],ME.TakeOff.S_TOFL);
    WingLoading.TakeOffSP = TOP25*sigma*DP.CLmax_TO*T_W_TO; %[kg/m^2]
    
    
%Obtain static constant from similar planes
    figure()
    plot(loadFields(SP,'Actuations.Sto'),TOP25_SP,'*','Color',Parameters.Colors(1,:)); hold on
    plot([min(loadFields(SP,'Actuations.Sto')),max(loadFields(SP,'Actuations.Sto'))],(loadFields(SP,'Actuations.Sto')'\TOP25_SP').*...
         [min(loadFields(SP,'Actuations.Sto')),max(loadFields(SP,'Actuations.Sto'))],'--','Color',Parameters.Colors(2,:))
    plot(ME.TakeOff.S_TOFL,TOP25R*CF.lbf2N/CST.GravitySI/(CF.ft2m^2),'+','Color',Parameters.Colors(3,:))
    plot(ME.TakeOff.S_TOFL,TOP25,'o','Color',Parameters.Colors(4,:))
    xlabel('Take-Off Field Length (TOFL) [m]');
    ylabel('Take-Off Parameter [kg/m^2]');
    h = findobj(gca,'Type','line');
    legend([h(4),h(2),h(1)],{'Similar Planes','Roskam Method','Similar Planes Method'},'Location','southeast')
    legend('boxoff')
    saveFigure(ME.FiguresFolder,'TOP25')
    
clear sigma h TOP25 TOP25_SP TOP25R



%% 3.3 SIZING TO LANDING DISTANCE REQUIREMENTS
%Roskam Method
    WL_WTO_R = DP.MLW_MTOW;
    Vapp_R   = sqrt(ME.Landing.S_LFL*CF.m2ft/0.3); %[kts] --> Roskam eq 3.16
    VStall_L = Vapp_R/1.3;                 %[kts] --> Roskam eq 3.15
    WL_Sw    = (rho_L*DP.CLmax_L*(VStall_L*CF.kts2ms)^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingRoskam = WL_Sw/WL_WTO_R;

%Similar Planes Method
    WL_WTO_SP = DP.MLW_MTOW;
    Vapp_SP   = sqrt((loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*ME.Landing.S_LFL);
    VStall_L  = Vapp_SP/1.3; %[m/s]
    WL_Sw     = (rho_L*DP.CLmax_L*VStall_L^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingSP = WL_Sw/WL_WTO_SP;

%Obtain static constant from similar planes
    figure()
    plot(loadFields(SP,'Actuations.Sl'),loadFields(SP,'Actuations.Vapprox').^2,'*','Color',Parameters.Colors(1,:)); hold on
    plot([450,max(loadFields(SP,'Actuations.Sl'))],(loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*...
         [450,max(loadFields(SP,'Actuations.Sl'))],'--','Color',Parameters.Colors(2,:))
    plot(ME.Landing.S_LFL,(Vapp_R*CF.kts2ms).^2,'+','Color',Parameters.Colors(3,:))
    plot(ME.Landing.S_LFL,Vapp_SP.^2,'o','Color',Parameters.Colors(4,:))
    xlabel('Landing Field Length (LFL) [m]','Interpreter','Tex');
    ylabel('Square of Approach Speed   V_{approx}^2   [m^2/s^2]');
    h = findobj(gca,'Type','line');
    legend([h(4),h(2),h(1)],{'Similar Planes','Roskam Method','Similar Planes Method'},'Location','southeast')
    legend('boxoff')
    saveFigure(ME.FiguresFolder,'Vapprox')

clear WL_WTO_R WL_WTO_SP Vapp_R Vapp_SP VStall_L WL_Sw h



%% 3.4 SIZING TO CLIMB REQUIREMENTS --> pag 143
%TAKE-OFF REQUIREMENTS
    % CS 25.121 (OEI) --> Transition segment climb requirement
        %Climb Gradient: Positive% for 2 engines, 0.3% for 3 engines, 0.5% for 4 engines
        %Flaps: Take-Off
        %Landing Gear: Down
        %Speed: Between Vliftoff=1.1VstallTO and V2=1.13Vstall_TO --> Vliftoff
        %Ground effect: Active
        %Weight: Maximum Take-Off Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch DP.NumberEngines
            case 2
                CGR = 0.000;
            case 3
                CGR = 0.003;
            case 4
                CGR = 0.005;
        end
        CL = DP.CLmax_TO/(1.1^2);
        CD = polyval(Parameters.Polar.TakeOffGearDown,CL);
        Efficiency = CL/CD;
        ThrustWeight_TO.TransitionSegment = (DP.NumberEngines/(DP.NumberEngines-1))*(CGR + 1/Efficiency);
        ThrustWeight_TO.TransitionSegment = ThrustWeight_TO.TransitionSegment/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
       
    % CS 25.111 (OEI) --> Initial climb segment requirement
        %Climb Gradient: 1.2% for 2 engines, 1.5% for 3 engines, 1.7% for 4 engines
        %Flaps: Take-Off
        %Landing Gear: Up
        %Speed: V2=1.13Vstall_TO
        %Ground effect: Active
        %Weight: Maximum Take-Off Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch DP.NumberEngines
            case 2
                CGR = 0.012;
            case 3
                CGR = 0.015;
            case 4
                CGR = 0.017;
        end
        CL = DP.CLmax_TO/(1.13^2);
        CD = polyval(Parameters.Polar.TakeOffGearUp,CL);
        Efficiency = CL/CD;
        ThrustWeight_TO.InitialSegment = (DP.NumberEngines/(DP.NumberEngines-1))*(CGR + 1/Efficiency);
        ThrustWeight_TO.InitialSegment = ThrustWeight_TO.InitialSegment/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
        
    % CS 25.121 (OEI) --> Second segment climb requirement
        %Climb Gradient: 2.4% for 2 engines, 2.7% for 3 engines, 3.0% for 4 engines
        %Flaps: Take-Off
        %Landing Gear: Up
        %Speed: V2=1.13Vstall_TO
        %Ground effect: Not active
        %Weight: Maximum Take-Off Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch DP.NumberEngines
            case 2
                CGR = 0.024;
            case 3
                CGR = 0.027;
            case 4
                CGR = 0.030;
        end
        CL = DP.CLmax_TO/(1.13^2);
        CD = polyval(Parameters.Polar.TakeOffGearUp,CL);
        Efficiency = CL/CD;
        ThrustWeight_TO.SecondSegment = (DP.NumberEngines/(DP.NumberEngines-1))*(CGR + 1/Efficiency);
        ThrustWeight_TO.SecondSegment = ThrustWeight_TO.SecondSegment/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
        
    % CS 25.121 (OEI) --> En-Route segment climb requirement
        %Climb Gradient: 1.2% for 2 engines, 1.5% for 3 engines, 1.7% for 4 engines
        %Flaps: Retracted
        %Landing Gear: Up
        %Speed: VFTO=1.18Vstall_TO
        %Ground effect: Not active
        %Weight: Maximum Take-Off Weight
        %Thrust: Continuous Thrust
        %Engines: N-1
        switch DP.NumberEngines
            case 2
                CGR = 0.012;
            case 3
                CGR = 0.015;
            case 4
                CGR = 0.017;
        end
        CL = DP.CLmax/(1.18^2);
        CD = polyval(Parameters.Polar.LowSpeedCruise,CL);
        Efficiency = CL/CD;
        ThrustWeight_TO.EnRouteSegment = (DP.NumberEngines/(DP.NumberEngines-1))*(CGR + 1/Efficiency);
        ThrustWeight_TO.EnRouteSegment = ThrustWeight_TO.EnRouteSegment/0.94; %Ratio between continuous thrust and take-off thrust --> 0.94 (Roskam)
        ThrustWeight_TO.EnRouteSegment = ThrustWeight_TO.EnRouteSegment/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
        
        
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
        switch DP.NumberEngines
            case 2
                CGR = 0.032;
            case 3
                CGR = 0.032;
            case 4
                CGR = 0.032;
        end
        CL = DP.CLmax_L/(1.3^2);
        CD = polyval(Parameters.Polar.LandingGearDown,CL);
        Efficiency = CL/CD;
        T_W_L = (CGR + 1/Efficiency);
        ThrustWeight_TO.BalkedLanding = T_W_L/DP.MLW_MTOW;
        ThrustWeight_TO.BalkedLanding = ThrustWeight_TO.BalkedLanding/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
        
    % CS 25.121 (OEI) --> Balked Approach requirements
        %Climb Gradient: 2.1% for 2 engines, 2.4% for 3 engines, 2.7% for 4 engines
        %Flaps: Approach
        %Landing Gear: Up
        %Speed: VstallAprox<=1.1Vstall_L && V<=1.4VstallAprox
        %Ground effect: Not active
        %Weight: Maximum Landing Weight
        %Thrust: Take-Off Thrust
        %Engines: N-1
        switch DP.NumberEngines
            case 2
                CGR = 0.021;
            case 3
                CGR = 0.024;
            case 4
                CGR = 0.027;
        end
        CL = DP.CLmax_L/((1.1*1.4)^2);
        CD = polyval([mean([Parameters.Polar.TakeOffGearUp(1),Parameters.Polar.LandingGearUp(1)]),0,...
                      mean([Parameters.Polar.TakeOffGearUp(3),Parameters.Polar.LandingGearUp(3)])],CL);
        Efficiency = CL/CD;
        T_W_L = (DP.NumberEngines/(DP.NumberEngines-1))*(CGR + 1/Efficiency);
        ThrustWeight_TO.BalkedApproach = T_W_L/DP.MLW_MTOW;
        ThrustWeight_TO.BalkedApproach = ThrustWeight_TO.BalkedApproach/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
        
 clear CGR CL CD Efficiency T_W_L              



%% DESIGN POINT
    figure()
    LegendSTR=cell(0);
    
    %Take-Off Field Length
        plot(WingLoading.TakeOffRoskam,T_W_TO,'Color',Parameters.Colors(1,:)); hold on;
            LegendSTR{end+1} = 'Take-Off (Roskam)';
            
        plot(WingLoading.TakeOffSP,T_W_TO,'Color',Parameters.Colors(2,:));
            LegendSTR{end+1} = 'Take-Off (Similar Planes)';
            
    %Take-Off Stall Speed
        if isfield(WingLoading,'Stall_TO')
            plot(WingLoading.Stall_TO.*ones(1,length(T_W_TO)),T_W_TO,'Color',Parameters.Colors(3,:));
                LegendSTR{end+1} = 'Take-Off Stall Speed';
        end
        
    %Landing Field Length
        plot(WingLoading.LandingRoskam.*ones(1,length(T_W_TO)),T_W_TO,'Color',Parameters.Colors(4,:));
            LegendSTR{end+1} = 'Landing (Roskam)';
            
        plot(WingLoading.LandingSP.*ones(1,length(T_W_TO)),T_W_TO,'Color',Parameters.Colors(5,:));
            LegendSTR{end+1} = 'Landing (Similar Planes)';
            
    %Climb Requirements
        %Take-Off
        plot(W_S_TO,ThrustWeight_TO.TransitionSegment.*ones(1,length(W_S_TO)),'Color',Parameters.Colors(6,:));
            LegendSTR{end+1} = 'Climb (Transition Segment)';

        plot(W_S_TO,ThrustWeight_TO.InitialSegment.*ones(1,length(W_S_TO)),'Color',Parameters.Colors(7,:));
            LegendSTR{end+1} = 'Climb (Initial Segment)';

        plot(W_S_TO,ThrustWeight_TO.SecondSegment.*ones(1,length(W_S_TO)),'Color',Parameters.Colors(8,:));
            LegendSTR{end+1} = 'Climb (Second Segment)';
            
        plot(W_S_TO,ThrustWeight_TO.EnRouteSegment.*ones(1,length(W_S_TO)),'Color',Parameters.Colors(9,:));
            LegendSTR{end+1} = 'Climb (En-Route Segment)';
            
        %Landing
        plot(W_S_TO,ThrustWeight_TO.BalkedLanding.*ones(1,length(W_S_TO)),'Color',Parameters.Colors(10,:));
            LegendSTR{end+1} = 'Climb (Balked Landing)';
            
        plot(W_S_TO,ThrustWeight_TO.BalkedApproach.*ones(1,length(W_S_TO)),'Color',Parameters.Colors(11,:));
            LegendSTR{end+1} = 'Climb (Balked Approach)';
            
    %Formating
        xlim([0,max(W_S_TO)])
        xlabel('Wing Loading - MTOW/Sw [kg/m^2]')
        ylabel('Thrust/Weight_{TO} [-]')
        grafWidth   = 16;
        grafAR      = 0.6;
        set(gcf,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
        warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
        [h,objs]=columnlegend(2,LegendSTR,'Location','northwest','FontSize',8,'boxoff');
        drawnow; % make sure everything is finished rendering 
        set(findall(objs, 'type', 'text'), 'fontsize', 8, 'interpreter', 'tex')
        warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');
        clear grafWidth grafAR
        
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
    warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
    saveFigure(ME.FiguresFolder,'DesignPoint')
    warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');

    clear choiceFlag choice p LegendStr

    
%% 3.5 Sizing to maneuvering requierements

%% 3.6 Sizing to cruise speed requirementes **
% [T, a, P, rho] = atmosisa(ME.Cruise.Altitude);
% q = 0.5*rho*ME.Cruise.Speed^2;
% 
% Wcr_WTO = 1;
% for i = [1,2,3,4]
% Wcr_WTO =  Wcr_WTO*Parameters.fuelFraction(i).value; % Wcrucero entre mtow
% end
% 
% Pto_Pcr = 0.65; %empuje al despegue versus empuje crucero
% 
% Wcr_S = Wcr_WTO * Wto_S;
% 
% Tcr_Wcr_cr = (C_D0*q)./Wcr_S + k.*Wcr_S./q;
% Tto_Wto_cr = Tcr_Wcr_cr * Wcr_WTO;
% % Tto_Wto_cr = ((C_D0*q)./Wto_S + Wto_S.*(Wcr_WTO^2*k)); %?
% P_W.cr = Pto_Pcr.*Tto_Wto_cr.*ME.Cruise.Speed./Parameters.Cruise.n_p
% 
% plot(Wto_S,P_W.cr)
% 
% 

% 
% figure(2)
% hold all
% plot(Wto_S,(C_D0*q)./Wto_S)
%  plot(Wto_S, Wto_S.*(Wcr_WTO^2*k))


% S_wet = 10^(parameters.Table_3_5.c + parameters.Table_3_5.d * ...
% log10(W_TO))  Usar si las constantes se calculan en SI






clear rho rho0 rho_TO rho_L
