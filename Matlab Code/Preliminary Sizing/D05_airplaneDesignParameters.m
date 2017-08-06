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
    if isfield(DP,'WingLoading')
        W_S_temp = DP.WingLoading;
    else
        W_S_temp = mean(loadFields(SP,'Wing.WingLoading'),'omitnan'); %Wing loading from similar airplanes
    end

    %Roskam --> Eq 3.22
    S_wet = (CF.ft2m^2)*10^(Parameters.Table_3_5.c + Parameters.Table_3_5.d*log10(AC.Weight.MTOW*CST.GravitySI*CF.N2lbf)); %[m^2]
    %Roskam --> Eq 3.21
    f = (CF.ft2m^2)*10^(Parameters.Table_3_4.a + Parameters.Table_3_4.b * log10(S_wet*(CF.m2ft^2))); %[m^2]

    Sw =  AC.Weight.MTOW/ W_S_temp;
    C_D0 = f/Sw; 

    %Compresibility corrections: Torenbeek pag.172 pdf
        %deltaC_D0 =0.0005 long-range cruise conditions
        %deltaC_D0 =0.0020 high-speed cruise conditions
    %Compresibility corrections for a Boeing 727: Roskam Part I pag.166
        %deltaC_D0 =0.0000 Mach = 0.70
        %deltaC_D0 =0.0001 Mach = 0.75
        %deltaC_D0 =0.0003 Mach = 0.80
        %deltaC_D0 =0.0011 Mach = 0.85
        %deltaC_D0 =0.0030 Mach = 0.90
    
    %Polar coefficients as function of CL to use with polyval --> CD = CD0 + k*CL^2
    Parameters.Polar.LowSpeedCruise  = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0];
    Parameters.Polar.LongRangeCruise = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0+0.0005]; %Compresibility effects from Torenbeek
    Parameters.Polar.HighSpeedCruise = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.clean),0,C_D0+0.0020]; %Compresibility effects from Torenbeek
    Parameters.Polar.TakeOffGearUp   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.take_off_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.take_off_flaps];
    Parameters.Polar.TakeOffGearDown = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.take_off_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.take_off_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear];
    Parameters.Polar.LandingGearUp   = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.landing_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.landing_flaps];
    Parameters.Polar.LandingGearDown = [inv(pi*DP.AspectRatio*Parameters.Table_3_6.e.landing_flaps),0,C_D0+Parameters.Table_3_6.deltaC_D0.landing_flaps+Parameters.Table_3_6.deltaC_D0.landing_gear];
    
    clear f C_D0 W_S_temp S_wet Sw
    
    
    
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
    TOP25R = DP.TOFL*CF.m2ft/37.5; %[lbf/ft^2]
    WingLoading.TakeOffRoskam = TOP25R*sigma*DP.CLmax_TO*T_W_TO*(CF.lbf2N/(CF.ft2m^2))/CST.GravitySI; %[kg/m^2]
    
%Similar planes method
    TOP25_SP = loadFields(SP,'Wing.WingLoading')./(sigma.*loadFields(SP,'Wing.CLmax_TO').*loadFields(SP,'Weight.Tto_MTOW')./CST.GravitySI); %[kg/m^2]
    TOP25=polyval([loadFields(SP,'Actuations.Sto')'\TOP25_SP',0],DP.TOFL);
    WingLoading.TakeOffSP = TOP25*sigma*DP.CLmax_TO*T_W_TO; %[kg/m^2]
    
    
%Obtain static constant from similar planes
if DP.ShowReportFigures
    figure()
    plot(loadFields(SP,'Actuations.Sto'),TOP25_SP,'*','Color',Parameters.Colors(1,:)); hold on
    m = loadFields(SP,'Actuations.Sto')'\TOP25_SP';
    plot([min(loadFields(SP,'Actuations.Sto')),max(loadFields(SP,'Actuations.Sto'))],...
    m.*[min(loadFields(SP,'Actuations.Sto')),max(loadFields(SP,'Actuations.Sto'))],'--','LineWidth',1.5,'Color',Parameters.Colors(2,:))
    plot(DP.TOFL,TOP25R*CF.lbf2N/CST.GravitySI/(CF.ft2m^2),'+','LineWidth',2,'Color',Parameters.Colors(3,:))
    plot(DP.TOFL,TOP25,'o','LineWidth',2,'Color',Parameters.Colors(4,:))
    xlabel('Take-Off Field Length (TOFL) [m]','Interpreter','latex');
    ylabel('Take-Off Parameter (TOP) [kg/m^2]','Interpreter','latex');
    h = findobj(gca,'Type','line');
    legend([h(4),h(2),h(1)],{'Similar Planes','Roskam Method','Similar Planes Method'},'Location','southeast')
    legend('boxoff')
    x = min(loadFields(SP,'Actuations.Sto')) + 0.35*(max(loadFields(SP,'Actuations.Sto'))-min(loadFields(SP,'Actuations.Sto')));
    y = m*x;
    txt0 = '$$\ \ \leftarrow$$ $$TOP=A \cdot TOFL$$';
    txt1 = strcat('$$\ \ \ \ \ \ A=',num2str(m),'$$');
    txt2 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(getR2(loadFields(SP,'Actuations.Sto'),TOP25_SP,[m,0])),'$$');
    text(x,y,txt0,'Interpreter','latex','FontSize',11)
    text(x,y-25,txt1,'Interpreter','latex','FontSize',11)
    text(x,y-50,txt2,'Interpreter','latex','FontSize',11)
    clear x y m
    saveFigure(ME.FiguresFolder,'TOP25')
end
    
clear sigma h TOP25 TOP25_SP TOP25R



%% 3.3 SIZING TO LANDING DISTANCE REQUIREMENTS
%Roskam Method
    WL_WTO_R = DP.MLW_MTOW;
    Vapp_R   = sqrt(DP.LFL*CF.m2ft/0.3); %[kts] --> Roskam eq 3.16
    VStall_L = Vapp_R/1.3;               %[kts] --> Roskam eq 3.15
    WL_Sw    = (rho_L*DP.CLmax_L*(VStall_L*CF.kts2ms)^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingRoskam = WL_Sw/WL_WTO_R;

%Similar Planes Method
    WL_WTO_SP = DP.MLW_MTOW;
    Vapp_SP   = sqrt((loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*DP.LFL);
    VStall_L  = Vapp_SP/1.3; %[m/s]
    WL_Sw     = (rho_L*DP.CLmax_L*VStall_L^2)/(2*CST.GravitySI); %[kg/m^2]
    WingLoading.LandingSP = WL_Sw/WL_WTO_SP;

%Obtain static constant from similar planes
if DP.ShowReportFigures
    figure()
    plot(loadFields(SP,'Actuations.Sl'),loadFields(SP,'Actuations.Vapprox').^2,'*','Color',Parameters.Colors(1,:)); hold on
    m = loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)';
    plot([450,max(loadFields(SP,'Actuations.Sl'))],m.*[450,max(loadFields(SP,'Actuations.Sl'))],'--','LineWidth',1.5,'Color',Parameters.Colors(2,:))
    plot(DP.LFL,(Vapp_R*CF.kts2ms).^2,'+','LineWidth',2,'Color',Parameters.Colors(3,:))
    plot(DP.LFL,Vapp_SP.^2,'o','LineWidth',2,'Color',Parameters.Colors(4,:))
    xlabel('Landing Field Length (LFL) [m]','Interpreter','Tex');
    ylabel('Square of Approach Speed   V_{approx}^2   [m^2/s^2]');
    h = findobj(gca,'Type','line');
    legend([h(4),h(2),h(1)],{'Similar Planes','Roskam Method','Similar Planes Method'},'Location','southeast')
    legend('boxoff')
    x = 525;
    y = m*x;
    txt0 = '$$\ \ \leftarrow$$ $$V_{approx}^2=A \cdot LFL$$';
    txt1 = strcat('$$\ \ \ \ \ \ A=',num2str(m),'$$');
    txt2 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(getR2(loadFields(SP,'Actuations.Sl'),loadFields(SP,'Actuations.Vapprox').^2,[m,0])),'$$');
    text(x,y,txt0,'Interpreter','latex','FontSize',11)
    text(x,y-150,txt1,'Interpreter','latex','FontSize',11)
    text(x,y-300,txt2,'Interpreter','latex','FontSize',11)
    clear x y m txt0 txt1 txt2 txt3
    saveFigure(ME.FiguresFolder,'Vapprox')
end

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
        switch DP.EngineNumber
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
        ThrustWeight_TO.TransitionSegment = (DP.EngineNumber/(DP.EngineNumber-1))*(CGR + 1/Efficiency);
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
        switch DP.EngineNumber
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
        ThrustWeight_TO.InitialSegment = (DP.EngineNumber/(DP.EngineNumber-1))*(CGR + 1/Efficiency);
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
        switch DP.EngineNumber
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
        ThrustWeight_TO.SecondSegment = (DP.EngineNumber/(DP.EngineNumber-1))*(CGR + 1/Efficiency);
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
        switch DP.EngineNumber
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
        ThrustWeight_TO.EnRouteSegment = (DP.EngineNumber/(DP.EngineNumber-1))*(CGR + 1/Efficiency);
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
        switch DP.EngineNumber
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
        switch DP.EngineNumber
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
        T_W_L = (DP.EngineNumber/(DP.EngineNumber-1))*(CGR + 1/Efficiency);
        ThrustWeight_TO.BalkedApproach = T_W_L/DP.MLW_MTOW;
        ThrustWeight_TO.BalkedApproach = ThrustWeight_TO.BalkedApproach/0.8;  %28ºC (50ºF) Correction to take into account the worst case possible (CS 25.101)
        
 clear CGR CL CD Efficiency T_W_L              




%% 3.5 SIZING TO MANEUVERING REQUIREMENTS
%Para aviones de transporte comercial y privado, no es restrictiva, se usa mucho para aviones militares



%% 3.6 SIZING TO CRUISE SPEED REQUIREMENTS
%Estimate compresibility effects from Roskam depending on max speed
    %deltaC_D0 =0.0000 Mach = 0.70
    %deltaC_D0 =0.0001 Mach = 0.75
    %deltaC_D0 =0.0003 Mach = 0.80
    %deltaC_D0 =0.0011 Mach = 0.85
    %deltaC_D0 =0.0030 Mach = 0.90
    
    Parameters.deltaCD0_Compresibility = [0.70, 0.0000;
                                          0.75, 0.0001;
                                          0.80, 0.0003;
                                          0.85, 0.0011;
                                          0.90, 0.0030];
    [~,a,~,~] = atmosisa(DP.CruiseAltitude);
    if DP.MaxSpeed/a < Parameters.deltaCD0_Compresibility(1,1)
        deltaCD0 = 0;
    elseif DP.MaxSpeed/a > Parameters.deltaCD0_Compresibility(end,1)
        deltaCD0 = inf;
    else
        deltaCD0  = interp1(Parameters.deltaCD0_Compresibility(:,1),Parameters.deltaCD0_Compresibility(:,2),DP.MaxSpeed/a);
    end

    Parameters.Polar.MaxSpeedCruise = [Parameters.Polar.LowSpeedCruise(1),0,Parameters.Polar.LowSpeedCruise(3)+deltaCD0];
    
    W_WTO = prod([Parameters.fuelFraction(1:4).value]);
    sigma = rho / rho0;
    T_T0 = (sigma.^0.7).*(1-exp((DP.CruiseAltitude-17000)./2000)); %<-- Modelo de variación del empuje con la altura, Flight Dynamics - Robert F. Stengel
    ThrustWeight_TO.MaxSpeedCruise = (1./T_T0).*((0.5.*rho.*DP.MaxSpeed^2.*Parameters.Polar.MaxSpeedCruise(1)./(W_S_TO.*CST.GravitySI))+((W_WTO.^2).*...
                                     (2.*W_S_TO.*CST.GravitySI.*Parameters.Polar.MaxSpeedCruise(3)./(rho.*DP.MaxSpeed.^2))));
    
    
    clear a deltaCD0 W_WTO sigma T_T0



%% SIZE TO CEILING REQUIREMENTS


%% SIZE TO GUST REQUIREMENTS

options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
[WingLoading.Gust,~,exitflag,~] = fsolve(@(x)gustWingLoading(x, DP.AspectRatio, CST, AC, CF, ME), 1000, options); %WingLoading in N/m^2
WingLoading.Gust = WingLoading.Gust ./ CST.GravitySI; %WingLoading in kg/m^2
if ~isequal(exitflag,1)
    disp('El solver de la limitacion por rafagas no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear exitflag options
end


%% DESIGN POINT
if DP.ShowReportFigures
    figure(); hold on;
    LegendStr=cell(0);
    
    %Take-Off Field Length
    if DP.showRoskamRequirements == true
        plot(WingLoading.TakeOffRoskam,T_W_TO,'--','LineWidth',1.25,'Color',Parameters.Colors(1,:));
            LegendStr{end+1} = 'Take-Off (Roskam)';
    end
            
        plot(WingLoading.TakeOffSP,T_W_TO,'LineWidth',1.25,'Color',Parameters.Colors(2,:));
            LegendStr{end+1} = 'Take-Off (Similar Planes)';
            
    %Take-Off Stall Speed
        if isfield(WingLoading,'Stall_TO')
            plot(WingLoading.Stall_TO.*ones(1,length(T_W_TO)),T_W_TO,'LineWidth',1.25,'Color',Parameters.Colors(3,:));
                LegendStr{end+1} = 'Take-Off Stall Speed';
        end
        
    %Landing Field Length
    if DP.showRoskamRequirements == true
        plot(WingLoading.LandingRoskam.*ones(1,length(T_W_TO)),T_W_TO,'--','LineWidth',1.25,'Color',Parameters.Colors(4,:));
            LegendStr{end+1} = 'Landing (Roskam)';
    end
            
        plot(WingLoading.LandingSP.*ones(1,length(T_W_TO)),T_W_TO,'LineWidth',1.25,'Color',Parameters.Colors(5,:));
            LegendStr{end+1} = 'Landing (Similar Planes)';
            
    %Climb Requirements
        %Take-Off
        plot(W_S_TO,ThrustWeight_TO.TransitionSegment.*ones(1,length(W_S_TO)),'LineWidth',1.25,'Color',Parameters.Colors(6,:));
            LegendStr{end+1} = 'Climb (Transition Segment)';

        plot(W_S_TO,ThrustWeight_TO.InitialSegment.*ones(1,length(W_S_TO)),'LineWidth',1.25,'Color',Parameters.Colors(7,:));
            LegendStr{end+1} = 'Climb (Initial Segment)';

        plot(W_S_TO,ThrustWeight_TO.SecondSegment.*ones(1,length(W_S_TO)),'LineWidth',1.25,'Color',Parameters.Colors(8,:));
            LegendStr{end+1} = 'Climb (Second Segment)';
            
        plot(W_S_TO,ThrustWeight_TO.EnRouteSegment.*ones(1,length(W_S_TO)),'LineWidth',1.25,'Color',Parameters.Colors(9,:));
            LegendStr{end+1} = 'Climb (En-Route Segment)';
            
        %Landing
        plot(W_S_TO,ThrustWeight_TO.BalkedLanding.*ones(1,length(W_S_TO)),'LineWidth',1.25,'Color',Parameters.Colors(10,:));
            LegendStr{end+1} = 'Climb (Balked Landing)';
            
        plot(W_S_TO,ThrustWeight_TO.BalkedApproach.*ones(1,length(W_S_TO)),'LineWidth',1.25,'Color',Parameters.Colors(11,:));
            LegendStr{end+1} = 'Climb (Balked Approach)';
            
    %Max Speed Cruise
        plot(W_S_TO,ThrustWeight_TO.MaxSpeedCruise,'LineWidth',1.25,'Color',Parameters.Colors(12,:));
            LegendStr{end+1} = 'Max Speed Cruise';
            
    %Max gust factor
    plot(WingLoading.Gust.*ones(1,length(T_W_TO)),T_W_TO,'LineWidth',1.25,'Color',Parameters.Colors(13,:));
            LegendStr{end+1} = 'Gust Loading';
        
    %Available Engines
        for i=1:length(Parameters.EngineOptions)
            plot(W_S_TO,(DP.EngineNumber.*Parameters.EngineOptions(i).Thrust.*1e3./(AC.Weight.MTOW.*CST.GravitySI)).*ones(1,length(W_S_TO)),...
                 'LineWidth',1,'LineStyle','--','Color',Parameters.Colors(13+i,:));
            LegendStr{end+1} = ['Engine: ',char(Parameters.EngineOptions(i).Model)]; %#ok<SAGROW>
        end
        
    %Similar Planes
        plot(loadFields(SP,'Wing.WingLoading'),loadFields(SP,'Engine.TotalThrust')./(loadFields(SP,'Weight.MTOW').*CST.GravitySI),...
             '*','Color',Parameters.Colors(14+i,:))
        clear i
        
    %Formating
        xlim([200,max(W_S_TO)])
        ylim([0.20,0.75])
        set(gcf,'Position',[450   200   700   525])
        xlabel('Wing Loading - MTOW/S_{w} [kg/m^2]')
        ylabel('Thrust/Weight_{TO} [-]')
        grafWidth   = 16;
        grafAR      = 0.6;
        set(gcf,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
        warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
        [~,objs]=columnlegend(2,LegendStr,'Location','northwest','FontSize',8,'boxoff');
        drawnow; % make sure everything is finished rendering 
        set(findall(objs, 'type', 'text'),'FontName','Times new Roman','FontSize', 8, 'interpreter', 'tex')
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
        x =  382.5; %[kg/m^2] WingLoad
        for i=1:length(Parameters.EngineOptions)
            usedEngine(i) = strcmp(Parameters.EngineOptions(i).Model,DP.EngineModel); %#ok<SAGROW>
        end
        y =  DP.EngineNumber*Parameters.EngineOptions(find(usedEngine,1)).Thrust*1e3/(AC.Weight.MTOW*CST.GravitySI); %[-] Thrust to Weight ratio at take-off --> Snecma Silvercrest 2D
        plot(x,y,'o');
    end
    saveFigure(ME.FiguresFolder,'DesignPoint')
    clear choiceFlag choice p LegendStr i
else
    x =  382.5; %[kg/m^2] WingLoad
    for i=1:length(Parameters.EngineOptions)
        usedEngine(i) = strcmp(Parameters.EngineOptions(i).Model,DP.EngineModel);  %#ok<SAGROW>
    end
    y =  DP.EngineNumber*Parameters.EngineOptions(find(usedEngine,1)).Thrust*1e3/(AC.Weight.MTOW*CST.GravitySI); %[-] Thrust to Weight ratio at take-off --> Snecma Silvercrest 2D
end


    %Store design point into AC structure
    AC.Wing.WingLoading   = x;
    AC.Weight.Tto_MTOW    = y;
    AC.Wing.AspectRatio   = DP.AspectRatio;
    AC.Wing.Sw            = AC.Weight.MTOW/AC.Wing.WingLoading;
    AC.Wing.WingSpan      = sqrt(AC.Wing.AspectRatio*AC.Wing.Sw);
    AC.Wing.TaperRatio    = DP.TaperRatio;
    AC.Wing.RootChord     = (2/(1+AC.Wing.TaperRatio))*sqrt(AC.Wing.Sw/AC.Wing.AspectRatio);
    AC.Wing.TipChord      = AC.Wing.TaperRatio*AC.Wing.RootChord;
    AC.Wing.CMG           = AC.Wing.RootChord*((1+AC.Wing.TaperRatio)/2);
    AC.Wing.CMA           = (2/3)*AC.Wing.RootChord*((1+AC.Wing.TaperRatio+AC.Wing.TaperRatio^2)/(1+AC.Wing.TaperRatio));
    AC.Wing.CLdesign      = 2*AC.Weight.MTOW*prod([Parameters.fuelFraction(1:4).value])*CST.GravitySI/(rho*AC.Wing.Sw*DP.CruiseSpeed^2);
    AC.Wing.CLmax         = DP.CLmax;
    AC.Wing.CLmax_TO      = DP.CLmax_TO;
    AC.Wing.CLmax_L       = DP.CLmax_L;
    AC.Engine             = Parameters.EngineOptions(find(usedEngine,1));
    AC.Engine.Number      = DP.EngineNumber;
    AC.Engine.TotalThrust = AC.Weight.Tto_MTOW*(AC.Weight.MTOW*CST.GravitySI);
    AC.Engine.TotalWeight = AC.Engine.Number*AC.Engine.Weight;

    

clear rho rho0 rho_TO rho_L T_W_TO W_S_TO x y usedEngine ThrustWeight_TO WingLoading i


%% OTHER USEFUL FUNCTIONS
function [F] = gustWingLoading(x, A, CST, AC, CF, ME) %-->x=carga alar en N/m^2

    [~, asound, P, rho] = atmosisa(ME.Cruise.Altitude);
    [~, ~, ~, rho0] = atmosisa(0);

    % M = ME.Cruise.Speed/asound
    % beta = sqrt(1-M^2)
    % k = 4*pi/beta/(2*pi)
    % flecha=0; %RAD
    % a = 2*pi*A/(2+sqrt((A^2*beta^2/k^2)*(1+(tan(flecha))^2/beta^2)+4)) %DA MUY ALTO


    a = 3.8; %[Cl_alpha] Raymer p.311 High AR swept wing
    CMG =sqrt((AC.Weight.MTOW*CST.GravitySI)/(x*A));


    mu = 2*x/(rho*a*CMG*CST.GravitySI);
    kg = 0.88*mu/(5.3+mu);

    %Load Factor for limit maneuvering
    n = 2.1+(24e3/(10e3+AC.Weight.MTOW*CF.kg2lbm));
    if n<2.5
        n = 2.5;
    elseif n>3.8
        n = 3.8;
    end

    %Reference Gust Speed
    if ME.Cruise.Altitude*CF.m2ft<20e3
        U_EAS = 50*CF.ft2m;
    elseif ME.Cruise.Altitude*CF.m2ft>=20e3
        U_EAS = ((50 + (ME.Cruise.Altitude*CF.m2ft-20e3)*(25-50)/(50e3-20e3)))*CF.ft2m;
    end

    %Equivalent Cruise Airspeed
    V_EAS =  correctairspeed(ME.Cruise.Speed, asound, P, 'TAS', 'EAS');
    
    F = kg*rho0*U_EAS*V_EAS*a/(2*x)-(n-1);
    % W_S_Gust = rho0*a*U*V_EAS/(2*(n-1));

end
