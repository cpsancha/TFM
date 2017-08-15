%% WING CONFIGURATION --> Chapter 6 & 7 Roskam; Chapter 7 Torenbeek
% Selection of wing configuration and geometric characteristics
    % Cantilever wing (without braces)
    % Leading wing --> High[ ] / Low[ ]
    % Rear wing --> High[ ] / Low[ ]
    % Zero sweep[ ] / Positive sweep[X] / Negative sweep[ ]
    % Aspect ratio
    % Thickness ratio
    % Airfoils
    % Taper ratio
    % Twist
    % Incidence angle
    % Dihedral angle
    % High lift and control surface requirements
    % Winglets
    

    
%% AIRFOIL SELECTION (Normal to the wing)
    %Airfoil Model
        AC.Wing1.Airfoil.Name = 'NASA SC(3)-0712(B)'; %<--NASA Technical Memorandum 86371 (and 86370)
        AC.Wing2.Airfoil.Name = 'NASA SC(3)-0712(B)'; %<--NASA Technical Memorandum 86371 (and 86370)
    %Airfoil design Cl
        AC.Wing1.Airfoil.Cldesign = 0.7;
        AC.Wing2.Airfoil.Cldesign = 0.7;
    %Airfoil thickness ratio (t/c)
        AC.Wing1.Airfoil.t_c = 0.12;
        AC.Wing2.Airfoil.t_c = 0.12;
    %Mach of divergence
%         AC.Wing1.Airfoil.MachDiv = 0.726;
%         AC.Wing2.Airfoil.MachDiv = 0.726;
    %Cl_max
        AC.Wing1.Airfoil.Cl_max = 1.7;
        AC.Wing2.Airfoil.Cl_max = 1.7;
    %Transition point/chord
        AC.Wing1.Airfoil.transition_c = 0.05;
        AC.Wing2.Airfoil.transition_c = 0.05;
        

    
        
        
    
%% DEFINE GEOMETRICAL ASPECTS OF EACH TANDEM WING FROM THE GENERAL WING VALUES
%Surface
    AC.Wing1.Sw          = AC.Wing.Sw/2;
    AC.Wing2.Sw          = AC.Wing.Sw/2;
    
%Aspect Ratio
    AC.Wing1.AspectRatio = AC.Wing.AspectRatio;
    AC.Wing2.AspectRatio = AC.Wing.AspectRatio;
    
%Span
    AC.Wing1.WingSpan = sqrt(AC.Wing1.AspectRatio*AC.Wing1.Sw);
    AC.Wing2.WingSpan = sqrt(AC.Wing2.AspectRatio*AC.Wing2.Sw);
    
%Taper Ratio
    AC.Wing1.TaperRatio = DP.TaperRatio;
    AC.Wing2.TaperRatio = DP.TaperRatio;
    
%Twist at the tip
    AC.Wing1.TipTwist = DP.TipTwist;
    AC.Wing2.TipTwist = DP.TipTwist;
    
%Root Chord
    AC.Wing1.RootChord = (2/(1+AC.Wing1.TaperRatio))*sqrt(AC.Wing1.Sw/AC.Wing1.AspectRatio);
    AC.Wing2.RootChord = (2/(1+AC.Wing2.TaperRatio))*sqrt(AC.Wing2.Sw/AC.Wing2.AspectRatio);
    
%Tip Chord
    AC.Wing1.TipChord = AC.Wing1.TaperRatio*AC.Wing1.RootChord;
    AC.Wing2.TipChord = AC.Wing2.TaperRatio*AC.Wing2.RootChord;
    
%CMG
    AC.Wing1.CMG = AC.Wing1.RootChord*((1+AC.Wing1.TaperRatio)/2);
    AC.Wing2.CMG = AC.Wing2.RootChord*((1+AC.Wing2.TaperRatio)/2);
    
%CMA
    AC.Wing1.CMA = (2/3)*AC.Wing1.RootChord*((1+AC.Wing1.TaperRatio+AC.Wing1.TaperRatio^2)/(1+AC.Wing1.TaperRatio));
    AC.Wing2.CMA = (2/3)*AC.Wing2.RootChord*((1+AC.Wing2.TaperRatio+AC.Wing2.TaperRatio^2)/(1+AC.Wing2.TaperRatio));

%Sweep_14
    AC.Wing1.Sweep_14 = DP.Sweep_14;
    AC.Wing2.Sweep_14 = DP.Sweep_14;
    
%Sweep_12
    AC.Wing1.Sweep_12 = atand(tand(AC.Wing1.Sweep_14)+(4/AC.Wing1.AspectRatio)*((1-AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*(0.25-0.5));
    AC.Wing2.Sweep_12 = atand(tand(AC.Wing2.Sweep_14)+(4/AC.Wing2.AspectRatio)*((1-AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*(0.25-0.5));

%Sweep_LE
    AC.Wing1.Sweep_LE = atand(tand(AC.Wing1.Sweep_14)+(4/AC.Wing1.AspectRatio)*((1-AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*(0.25-0));
    AC.Wing2.Sweep_LE = atand(tand(AC.Wing2.Sweep_14)+(4/AC.Wing2.AspectRatio)*((1-AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*(0.25-0));
    
%Sweep_RE
    AC.Wing1.Sweep_RE = atand(tand(AC.Wing1.Sweep_14)+(4/AC.Wing1.AspectRatio)*((1-AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*(0.25-1));
    AC.Wing2.Sweep_RE = atand(tand(AC.Wing2.Sweep_14)+(4/AC.Wing2.AspectRatio)*((1-AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*(0.25-1));

%Dihedral
    AC.Wing1.Dihedral = DP.Dihedral;
    AC.Wing2.Dihedral = DP.Dihedral;
 
%Longitudinal Position of the leading edge at root
    AC.Wing1.Root_LE = DP.Wing1LongPos;
    AC.Wing2.Root_LE = AC.Wing1.Root_LE+AC.Wing1.RootChord+DP.Stagger;
    
%Y_CMA
    AC.Wing1.CMA_b = (AC.Wing1.WingSpan/6)*((1+2*AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio));
    AC.Wing2.CMA_b = (AC.Wing2.WingSpan/6)*((1+2*AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio));

%X_CMA_LE
    AC.Wing1.CMA_LE = AC.Wing1.Root_LE + (AC.Wing1.WingSpan/6)*((1+2*AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*tand(AC.Wing1.Sweep_LE);
    AC.Wing2.CMA_LE = AC.Wing2.Root_LE + (AC.Wing2.WingSpan/6)*((1+2*AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*tand(AC.Wing2.Sweep_LE);
    
%X_CMA_1/4
    AC.Wing1.CMA_14 = AC.Wing1.CMA_LE + AC.Wing1.CMA/4;
    AC.Wing2.CMA_14 = AC.Wing2.CMA_LE + AC.Wing2.CMA/4;

%X_ac    
    AC.Wing1.x_ac_w = AC.Wing1.CMA_14;
    AC.Wing2.x_ac_w = AC.Wing2.CMA_14;
    
%Fuel Volume
    AC.Wing1.fuelVolume = 0.54*AC.Wing1.Sw^2/AC.Wing1.WingSpan*AC.Wing1.Airfoil.t_c*((1+AC.Wing1.TaperRatio+AC.Wing1.TaperRatio^2)/(1+AC.Wing1.TaperRatio)^2);
    AC.Wing1.fuelVolume = 0.54*AC.Wing2.Sw^2/AC.Wing2.WingSpan*AC.Wing2.Airfoil.t_c*((1+AC.Wing2.TaperRatio+AC.Wing2.TaperRatio^2)/(1+AC.Wing2.TaperRatio)^2);

%Distance from Root_LE to Tip_LE
    AC.Wing1.TipSweep = (AC.Wing1.WingSpan/2)*tand(AC.Wing1.Sweep_LE);
    AC.Wing2.TipSweep = (AC.Wing2.WingSpan/2)*tand(AC.Wing2.Sweep_LE);

%Reynolds
        [rho,a,T,P,nu,~] = atmos(DP.CruiseAltitude); %#ok<ASGLU>
        CruiseMach1 = sqrt(DP.CruiseSpeed^2*Parameters.q1_qinf)/a;
        CruiseMach2 = sqrt(DP.CruiseSpeed^2*Parameters.q2_qinf)/a;
        Beta1       = sqrt(1-(CruiseMach1)^2); %Prandtl compresibility correction
        Beta2       = sqrt(1-(CruiseMach2)^2); %Prandtl compresibility correction
        AC.Wing1.Reynolds  = sqrt(DP.CruiseSpeed^2*Parameters.q1_qinf) * AC.Wing1.CMA / nu;
        AC.Wing2.Reynolds  = sqrt(DP.CruiseSpeed^2*Parameters.q2_qinf) * AC.Wing2.CMA / nu;
    
%Wing design CL
    AC.Wing1.CLdesign = AC.Wing1.Airfoil.Cldesign*(cosd(AC.Wing1.Sweep_14))^2;
    AC.Wing2.CLdesign = AC.Wing2.Airfoil.Cldesign*(cosd(AC.Wing2.Sweep_14))^2;

%Incidence of the wing respect the fuselage
    AC.Wing1.Incidence = DP.Incidence_1;
    AC.Wing2.Incidence = DP.Incidence_2;
    
%Wing AoA at root section
%     AC.Wing1.Root_AoA  = AC.Wing1.Incidence + AC.Fuselage.fuselage_AoA;
%     AC.Wing2.Root_AoA  = AC.Wing2.Incidence + AC.Fuselage.fuselage_AoA;

%Wing t_c
    AC.Wing1.t_c = AC.Wing1.Airfoil.t_c*cosd(AC.Wing1.Sweep_14);
    AC.Wing2.t_c = AC.Wing2.Airfoil.t_c*cosd(AC.Wing2.Sweep_14);
    
%Mach of divergence acording to maximum thickness, CL, and sweep angle
    options = optimoptions('fsolve',...
                           'StepTolerance',1e-9,...
                           'Display','none');
    [AC.Wing1.Airfoil.MachDiv,~,exitflag1,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, AC.Wing1.Airfoil.t_c, 0,...
                                               AC.Wing1.CLdesign,'Supercritical'),0.8,options);
    [AC.Wing1.Airfoil.MachDiv,~,exitflag2,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, AC.Wing2.Airfoil.t_c, 0,...
                                               AC.Wing2.CLdesign,'Supercritical'),0.8,options);
    [AC.Wing1.MachDiv,~,exitflag3,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, AC.Wing1.Airfoil.t_c, AC.Wing1.Sweep_14,...
                                               AC.Wing1.CLdesign,'Supercritical'),0.8,options);
    [AC.Wing2.MachDiv,~,exitflag4,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, AC.Wing2.Airfoil.t_c, AC.Wing2.Sweep_14,...
                                               AC.Wing2.CLdesign,'Supercritical'),0.8,options);
    if ~isequal(exitflag1,1)
        disp('El solver del mach de divergencia del perfil del ala 1 no ha logrado converger correctamente. Se debería revisar el resultado.')
        pause
    elseif ~isequal(exitflag2,1)
        disp('El solver del mach de divergencia del perfil del ala 2 no ha logrado converger correctamente. Se debería revisar el resultado.')
        pause
    elseif ~isequal(exitflag3,1)
        disp('El solver del mach de divergencia del ala 1 no ha logrado converger correctamente. Se debería revisar el resultado.')
        pause
    elseif ~isequal(exitflag4,1)
        disp('El solver del mach de divergencia del ala 2 no ha logrado converger correctamente. Se debería revisar el resultado.')
        pause    
    else
        clear exitflag1 exitflag2 exitflag3 exitflag4 options
    end

    %Show Divergence Mach depending on sweep
    if DP.ShowReportFigures
        sweepArray = [15,20,25,30,32.5,35,37.5];
        t_c = linspace(10,18,5);
        LegendStr = cell(0);
        options = optimoptions('fsolve',...
                               'StepTolerance',1e-9,...
                               'Display','none');
        figure()
        hold on
        for i=1:length(sweepArray)
            for j=1:length(t_c)
                [MachDiv(j),~,exitFlag,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, t_c(j)./100, sweepArray(i),AC.Wing1.CLdesign,'Supercritical'),0.8,options); %#ok<SAGROW>
                if ~isequal(exitFlag,1)
                    disp('El solver del mach de divergencia al generar la figura no ha logrado converger correctamente. Se debería revisar el resultado.')
                end
            end
            plot(t_c,MachDiv,'LineWidth',1.25,'Color',Parameters.Colors(i,:))
            LegendStr{end+1} = ['$\Lambda_{1/4}=',num2str(sweepArray(i)),'^o$']; %#ok<SAGROW>
        end
        plot(AC.Wing1.Airfoil.t_c*100,AC.Wing1.MachDiv,'o','LineWidth',1.25,'Color',Parameters.Colors(i+1,:))
        LegendStr{end+1}='Design Point';
        title('$M_{dd}\ en\ funcion\ de\ la\ flecha\ y\ el\ espesor\ relativo\ del\ perfil$','interpreter','latex')
        xlabel('$t/c\ [-]$','interpreter','latex')
        ylabel('$M_{dd}\ [-]$','interpreter','latex')
        legend(LegendStr,'Location','northeast','interpreter','latex')
        legend('boxoff')
        saveFigure(ME.FiguresFolder,'SweepDecision')
        clear sweepArray t_c i j exitFlag LegendStr MachDiv options
    end
    
    
    
%% LIFTING PROPERTIES OF AIRFOIL SECTIONS
%Wing1
    [AC.Wing1.Airfoil.Cl_alpha,...
     AC.Wing1.Airfoil.alpha_zeroLift,...
     AC.Wing1.Airfoil.Cm_ac,...
     AC.Wing1.Airfoil.Polar,~,~,~] = getAirfoilData(CruiseMach1*cosd(AC.Wing1.Sweep_12),AC.Wing1.Reynolds,Parameters.Colors,ME.FiguresFolder,false);
%Wing2
    [AC.Wing2.Airfoil.Cl_alpha,...
     AC.Wing2.Airfoil.alpha_zeroLift,...
     AC.Wing2.Airfoil.Cm_ac,...
     AC.Wing2.Airfoil.Polar,~,~,~] = getAirfoilData(CruiseMach2*cosd(AC.Wing2.Sweep_12),AC.Wing2.Reynolds,Parameters.Colors,ME.FiguresFolder,false);



%% WING LIFT AND LIFT DISTRIBUTIONS
%LIFT-CURVE SLOPE
%Torenbeek Method (Torenbeek E-4.1.b pag 473)
    EffectiveSweep1 = atand(tand(AC.Wing1.Sweep_12)/Beta1);
    EffectiveSweep2 = atand(tand(AC.Wing2.Sweep_12)/Beta2);
    k1 = (AC.Wing1.Airfoil.Cl_alpha)/(2*pi);
    k2 = (AC.Wing2.Airfoil.Cl_alpha)/(2*pi);
    CL_alpha_Torenbeek = (2*pi/Beta1)/((2/(Beta1*AC.Wing1.AspectRatio))+sqrt((1/((k1*cosd(EffectiveSweep1))^2))+(2/(Beta1*AC.Wing1.AspectRatio))^2));
%DATCOM Method - Straight-Tapered Wings - Method 1 - pag 4.1.3.2-3
    CL_alpha_DATCOM    = (2*pi*AC.Wing1.AspectRatio)/(2+sqrt((AC.Wing1.AspectRatio*Beta1/k1)^2*(1+(tand(AC.Wing1.Sweep_12)/Beta1)^2)+4));
%Polhamus Method - NACA Technical Note 3911
    CL_alpha_Polhamus  = (AC.Wing1.Airfoil.Cl_alpha*AC.Wing1.AspectRatio) / ((AC.Wing1.Airfoil.Cl_alpha/pi) + ...
                         sqrt((AC.Wing1.AspectRatio/cosd(AC.Wing1.Sweep_12))^2+(AC.Wing1.Airfoil.Cl_alpha/pi)^2 - ...
                         (AC.Wing1.AspectRatio*CruiseMach1)^2));
    if abs(CL_alpha_Torenbeek-CL_alpha_DATCOM) > 1e-6
        if ~ismember('AC:notDefinedCoG',ME.errorList)
            ME.errorList{end+1} = 'AC:CL_alpha_Torenbeek_DATCOM';
            warning('AC:CL_alpha_Torenbeek_DATCOM','La pendiente de la curva de sustentacion por el metodo de Torenbeek no es igual que por el metodo de las DATCOM. Elige cual se usa.')
        end
    elseif abs(CL_alpha_Torenbeek-CL_alpha_Polhamus) > 1e-6
        if ~ismember('AC:notDefinedCoG',ME.errorList)
            ME.errorList{end+1} = 'AC:CL_alpha_Torenbeek_DATCOM';
            warning('AC:CL_alpha_Torenbeek_Polhamus','La pendiente de la curva de sustentacion por el metodo de Torenbeek no es igual que por el metodo de Polhamus. Elige cual se usa.')
        end
    end
%Uses Torenbeek Method, because yes... to be decided when loaded airfoil data
    AC.Wing1.CL_alpha_w  = (2*pi/Beta1)/((2/(Beta1*AC.Wing1.AspectRatio))+sqrt((1/((k1*cosd(EffectiveSweep1))^2))+(2/(Beta1*AC.Wing1.AspectRatio))^2));
    AC.Wing2.CL_alpha_w  = (2*pi/Beta2)/((2/(Beta2*AC.Wing2.AspectRatio))+sqrt((1/((k2*cosd(EffectiveSweep2))^2))+(2/(Beta2*AC.Wing2.AspectRatio))^2));
    clear k1 k2 CL_alpha_Torenbeek CL_alpha_DATCOM CL_alpha_Polhamus
    
%SPANWISE LIFT DISTRIBUTION <-- Diederich Semiempirical Method from NACA TechNote 2751
%Create spanwise coordinates
    y1   = linspace(0,AC.Wing1.WingSpan/2,100);
    y2   = linspace(0,AC.Wing2.WingSpan/2,100);
%Adimensional spanwise coordinate
    AC.Wing1.eta = y1./(AC.Wing1.WingSpan/2);
    AC.Wing2.eta = y2./(AC.Wing2.WingSpan/2);
%Load graph from NACA Technical Note 2751
    run('.\Digitalized Data\Lift_Distribution_Function_f.m')
    EffectiveSweep1 = atand(tand(AC.Wing1.Sweep_14)/Beta1);
    EffectiveSweep2 = atand(tand(AC.Wing2.Sweep_14)/Beta2);
    if EffectiveSweep1 < 30
        f_low   = interp1(f_00deg(:,1),f_00deg(:,2),AC.Wing1.eta,'linear','extrap');
        f_hight = interp1(f_30deg(:,1),f_30deg(:,2),AC.Wing1.eta,'linear','extrap');
        f1 = f_low+(f_hight-f_low)./30.*EffectiveSweep1;
    elseif EffectiveSweep1 < 45
        f_low   = interp1(f_30deg(:,1),f_30deg(:,2),AC.Wing1.eta,'linear','extrap');
        f_hight = interp1(f_45deg(:,1),f_45deg(:,2),AC.Wing1.eta,'linear','extrap');
        f1 = f_low+(f_hight-f_low)./15.*(EffectiveSweep1-30);
    elseif EffectiveSweep1 < 60
        f_low   = interp1(f_45deg(:,1),f_45deg(:,2),AC.Wing1.eta,'linear','extrap');
        f_hight = interp1(f_60deg(:,1),f_60deg(:,2),AC.Wing1.eta,'linear','extrap');
        f1 = f_low+(f_hight-f_low)./15.*(EffectiveSweep1-45);
    else
        warning('Incorrect value for Effective Sweep of wing 1.')
    end

    if EffectiveSweep2 < 30
        f_low   = interp1(f_00deg(:,1),f_00deg(:,2),AC.Wing2.eta,'linear','extrap');
        f_hight = interp1(f_30deg(:,1),f_30deg(:,2),AC.Wing2.eta,'linear','extrap');
        f2 = f_low+(f_hight-f_low)./30.*EffectiveSweep2;
    elseif EffectiveSweep2 < 45
        f_low   = interp1(f_30deg(:,1),f_30deg(:,2),AC.Wing2.eta,'linear','extrap');
        f_hight = interp1(f_45deg(:,1),f_45deg(:,2),AC.Wing2.eta,'linear','extrap');
        f2 = f_low+(f_hight-f_low)./15.*(EffectiveSweep2-30);
    elseif EffectiveSweep2 < 60
        f_low   = interp1(f_45deg(:,1),f_45deg(:,2),AC.Wing2.eta,'linear','extrap');
        f_hight = interp1(f_60deg(:,1),f_60deg(:,2),AC.Wing2.eta,'linear','extrap');
        f2 = f_low+(f_hight-f_low)./15.*(EffectiveSweep2-45);
    else
        warning('Incorrect value for Effective Sweep of wing 2.')
    end
    f1(end) = 0;
    f2(end) = 0;
    
%Plan-Form parameter for Diederich Method in Torenbeek Fig E-5:
    DiederichCoordinate1 = (2*pi*AC.Wing1.AspectRatio)/(AC.Wing1.Airfoil.Cl_alpha*cosd(AC.Wing1.Sweep_14)); %Plan-Form parameter (F)
    DiederichCoordinate2 = (2*pi*AC.Wing2.AspectRatio)/(AC.Wing2.Airfoil.Cl_alpha*cosd(AC.Wing2.Sweep_14)); %Plan-Form parameter (F)

%Load digitized points, make a 3rd order polynomic regression (better
%extrapolation than higher order(situational but who cares)) and evaluate in the point of interest:
    run ('C1C2C3C4f0f30.m')
    C1_1 = polyval( polyfit(C1x,C1y,3) , DiederichCoordinate1);
    C2_1 = polyval( polyfit(C2x,C2y,3) , DiederichCoordinate1);
    C3_1 = polyval( polyfit(C3x,C3y,3) , DiederichCoordinate1);
    C4_1 = polyval( polyfit(C4x,C4y,3) , DiederichCoordinate1);
    C1_2 = polyval( polyfit(C1x,C1y,3) , DiederichCoordinate2);
    C2_2 = polyval( polyfit(C2x,C2y,3) , DiederichCoordinate2);
    C3_2 = polyval( polyfit(C3x,C3y,3) , DiederichCoordinate2);
    C4_2 = polyval( polyfit(C4x,C4y,3) , DiederichCoordinate2);  

%Calculate the chord
    AC.Wing1.c = getChord(y1, AC.Wing1.WingSpan, AC.Wing1.TaperRatio, AC.Wing1.Sw);
    AC.Wing2.c = getChord(y2, AC.Wing2.WingSpan, AC.Wing2.TaperRatio, AC.Wing2.Sw);
    
%Aditional lift distribution calculation (La)
    La1 = C1_1.*AC.Wing1.c./AC.Wing1.CMG + C2_1.*4./pi.*sqrt(1-AC.Wing1.eta.^2) + C3_1.*f1;
    La2 = C1_2.*AC.Wing2.c./AC.Wing2.CMG + C2_2.*4./pi.*sqrt(1-AC.Wing2.eta.^2) + C3_2.*f2;

 %Jone's edge velocity correction E=semiperimeter/semispan
    AC.Wing1.E = 1+((2*AC.Wing1.TaperRatio)/(AC.Wing1.AspectRatio*(1+AC.Wing1.TaperRatio)));
    AC.Wing2.E = 1+((2*AC.Wing2.TaperRatio)/(AC.Wing2.AspectRatio*(1+AC.Wing2.TaperRatio)));

%Linear Lofted (Geometrical) Twist    
%     Twist1 = AC.Wing1.TipTwist .* (AC.Wing1.TaperRatio.*AC.Wing1.eta./(1-(1-AC.Wing1.TaperRatio).*AC.Wing1.eta));
%     Twist2 = AC.Wing2.TipTwist .* (AC.Wing2.TaperRatio.*AC.Wing2.eta./(1-(1-AC.Wing2.TaperRatio).*AC.Wing2.eta));

%Linear Twist
    Twist1 = AC.Wing1.TipTwist .* AC.Wing1.eta;
    Twist2 = AC.Wing2.TipTwist .* AC.Wing2.eta;
    
%Local aerodynamic twist at the station for which Clb=0 if TipTwist=1º    
    twist_0lift_1 = -trapz(AC.Wing1.eta,(Twist1./AC.Wing1.TipTwist).*La1); %Torenbeek Eq. E-16 [adimensional, must be multiplied by TipTwist]
    twist_0lift_2 = -trapz(AC.Wing2.eta,(Twist2./AC.Wing2.TipTwist).*La2); %Torenbeek Eq. E-16 [adimensional, must be multiplied by TipTwist]
    
%Basic lift distribution (Lb)    
    Lb1 = La1.*C4_1.*cos(EffectiveSweep1).*((Twist1./AC.Wing1.TipTwist)+twist_0lift_1).*Beta1.*AC.Wing1.E;
    Lb2 = La2.*C4_2.*cos(EffectiveSweep2).*((Twist2./AC.Wing2.TipTwist)+twist_0lift_2).*Beta2.*AC.Wing2.E;

%Assumption of CL of the wing to calculate the real one
    CL_wing1 = 1;
    CL_wing2 = 1;

%Coefficient of Aditional lift
    Cla1 = CL_wing1.*AC.Wing1.CMG.*La1./AC.Wing1.c;
    Cla2 = CL_wing2.*AC.Wing2.CMG.*La2./AC.Wing2.c;
    
%Coefficient of Basic Lift
    Clb1 = Lb1.*deg2rad(AC.Wing1.TipTwist).*AC.Wing1.Airfoil.Cl_alpha.*AC.Wing1.CMG./(AC.Wing1.E.*AC.Wing1.c);
    Clb2 = Lb2.*deg2rad(AC.Wing2.TipTwist).*AC.Wing2.Airfoil.Cl_alpha.*AC.Wing2.CMG./(AC.Wing2.E.*AC.Wing2.c);

%Zero Lift Angle
    AC.Wing1.alpha_zeroLift = AC.Wing1.Airfoil.alpha_zeroLift + twist_0lift_1*AC.Wing1.TipTwist; %[degrees]
    AC.Wing2.alpha_zeroLift = AC.Wing2.Airfoil.alpha_zeroLift + twist_0lift_2*AC.Wing2.TipTwist; %[degrees]
    
%Maximum Wing Lift
    AC.Wing1.CLmax = min((AC.Wing1.Airfoil.Cl_max*cosd(AC.Wing1.Sweep_14)-Clb1)./Cla1);
    AC.Wing2.CLmax = min((AC.Wing2.Airfoil.Cl_max*cosd(AC.Wing2.Sweep_14)-Clb2)./Cla2);

%Total Lift Distribution
    Cl1  = AC.Wing1.CLmax.*Cla1 + Clb1;
    Cl2  = AC.Wing2.CLmax.*Cla2 + Clb2;

%Display Wing Lift Distribution    
if DP.ShowReportFigures
    %Wing 1
        figure()
        hold on
        [~,index] = max(Cl1);
        plot(AC.Wing1.eta(index),Cl1(index),'o','LineWidth',1.25,'Color',Parameters.Colors(1,:))
        plot(AC.Wing1.eta,ones(1,length(AC.Wing1.eta)).*AC.Wing1.Airfoil.Cl_max*cosd(AC.Wing1.Sweep_14),'--','LineWidth',1.25,'Color',Parameters.Colors(2,:));
        plot(AC.Wing1.eta,Cl1,'LineWidth',1.25,'Color',Parameters.Colors(3,:));
        plot(AC.Wing1.eta,Clb1,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
        plot(AC.Wing1.eta,AC.Wing1.CLmax.*Cla1,'LineWidth',1.25,'Color',Parameters.Colors(8,:));
        plot(AC.Wing1.eta,Cla1,'LineWidth',1.25,'Color',Parameters.Colors(4,:));
        legend('First point of stall','Maximum lift of the airfoil','Total lift distribution','Basic lift distribution','Aditional lift distribution','Aditional lift distribution for C_{L_{max}}=1','Location','southwest')
        legend('boxoff')
        xlabel('$\frac{y}{b/2}$','interpreter','latex')
        ylabel('$C_l$','interpreter','latex')
        title(['Spanwise Lift Distribution of wing 1 for $C_{Lmax}=',num2str(AC.Wing1.CLmax),'$ and $\varepsilon_t=',num2str(AC.Wing1.TipTwist),'^o$'],'interpreter','latex')
        saveFigure(ME.FiguresFolder,'SpanwiseLiftDistribution_1')
    %Wing 2
        figure()
        hold on
        [~,index] = max(Cl2);
        plot(AC.Wing2.eta(index),Cl2(index),'o','LineWidth',1.25,'Color',Parameters.Colors(1,:))
        plot(AC.Wing2.eta,ones(1,length(AC.Wing2.eta)).*AC.Wing2.Airfoil.Cl_max*cosd(AC.Wing2.Sweep_14),'--','LineWidth',1.25,'Color',Parameters.Colors(2,:));
        plot(AC.Wing2.eta,Cl2,'LineWidth',1.25,'Color',Parameters.Colors(3,:));
        plot(AC.Wing2.eta,Clb2,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
        plot(AC.Wing2.eta,AC.Wing2.CLmax.*Cla2,'LineWidth',1.25,'Color',Parameters.Colors(8,:));
        plot(AC.Wing2.eta,Cla2,'LineWidth',1.25,'Color',Parameters.Colors(4,:));
        legend('First point of stall','Maximum lift of the airfoil','Total lift distribution','Basic lift distribution','Aditional lift distribution','Aditional lift distribution for C_{L_{max}}=1','Location','southwest')
        legend('boxoff')
        xlabel('$\frac{y}{b/2}$','interpreter','latex')
        ylabel('$C_l$','interpreter','latex')
        title(['Spanwise Lift Distribution of wing 2 for $C_{Lmax}=',num2str(AC.Wing2.CLmax),'$ and $\varepsilon_t=',num2str(AC.Wing2.TipTwist),'^o$'],'interpreter','latex')
        saveFigure(ME.FiguresFolder,'SpanwiseLiftDistribution_2')
end
    

%Deflexión de estela para la segunda ala
%Adimensional parameters
    r = DP.Stagger / (AC.Wing1.WingSpan/2);
    m = DP.VerticalGap / (AC.Wing1.WingSpan/2);
%Downwash gradient
    deltaE_deltaAlpha = 1.75 * (AC.Wing1.CL_alpha_w/(pi*AC.Wing1.AspectRatio*(AC.Wing1.TaperRatio*r)^0.25*(1+abs(m))));
    %1-deltaE_deltaAlpa = 0.9*(1-deltaE_deltaAlpha) <-- 10% reduction acording to Torenbeek E-10.1

%Lift Coefficient in the linear range
    AC.Wing1.CL_w = AC.Wing1.CL_alpha_w * deg2rad(AC.Wing1.Incidence+AC.Fuselage.fuselage_AoA-AC.Wing1.Airfoil.alpha_zeroLift-twist_0lift_1*AC.Wing1.TipTwist);
    AC.Wing2.CL_w = AC.Wing2.CL_alpha_w * deg2rad(AC.Wing2.Incidence+AC.Fuselage.fuselage_AoA*0.9*(1-deltaE_deltaAlpha)-...
                    AC.Wing2.Airfoil.alpha_zeroLift-twist_0lift_2*AC.Wing2.TipTwist)+AC.Wing2.deltaCLdeltaE;

    

    
%% PITCHING MOMENT OF THE WING
%Contribution of the spanwise airfoil camber (Cm_ac_basic)
    Cm_ac_basic1 = 2/(AC.Wing1.Sw*AC.Wing1.CMA)*trapz(y1,AC.Wing1.Airfoil.Cm_ac.*AC.Wing1.c.^2);
    Cm_ac_basic2 = 2/(AC.Wing2.Sw*AC.Wing2.CMA)*trapz(y2,AC.Wing2.Airfoil.Cm_ac.*AC.Wing2.c.^2);
    
%Contribution of the basic lift distribution due to twist (deltaEpsilon*Cm_ac)    
    deltaEpsilonCm_ac1 = ((AC.Wing1.AspectRatio*AC.Wing1.CMG*tand(AC.Wing1.Sweep_14))/(2*AC.Wing1.CMA)).*trapz(AC.Wing1.eta,Clb1.*AC.Wing1.c.*AC.Wing1.eta./AC.Wing1.CMG);
    deltaEpsilonCm_ac2 = ((AC.Wing2.AspectRatio*AC.Wing2.CMG*tand(AC.Wing2.Sweep_14))/(2*AC.Wing2.CMA)).*trapz(AC.Wing2.eta,Clb2.*AC.Wing2.c.*AC.Wing2.eta./AC.Wing2.CMG);

%Pitching moment on the aerodynamic center (Cm_ac_w)    
    AC.Wing1.Cm_ac_w = Cm_ac_basic1 + deltaEpsilonCm_ac1;
    AC.Wing2.Cm_ac_w = Cm_ac_basic2 + deltaEpsilonCm_ac2;
    
%Pitching moment coefficient of the wings
    AC.Weight.x_cg = DP.x_cg;  %WARNING!!!  <-- TO BE REMOVED
    AC.Weight.y_cg = DP.y_cg;
    AC.Weight.z_cg = DP.z_cg;
    
    if ~ismember('AC:notDefinedCoG',ME.errorList)
        ME.errorList{end+1} = 'AC:notDefinedCoG';
        warning('AC:notDefinedCoG','Hay que poner bien la posición del centro de gravedad para los coeficientes de momentos')
    end
    
    AC.Wing1.Cm_w = AC.Wing1.Cm_ac_w + AC.Wing1.CL_w * ((AC.Weight.x_cg - AC.Wing1.CMA_14)/AC.Wing1.CMA);
    AC.Wing2.Cm_w = AC.Wing2.Cm_ac_w + AC.Wing2.CL_w * ((AC.Weight.x_cg - AC.Wing2.CMA_14)/AC.Wing2.CMA);
    
    

%% WING/FUSELAGE INTERFERENCE EFFECTS ON LIFT
%Net wing surface (Without fuselage)
    AC.Wing1.Snet = AC.Wing1.Sw - ((AC.Wing1.RootChord+getChord(AC.Fuselage.fusWidth/2,AC.Wing1.WingSpan,AC.Wing1.TaperRatio,AC.Wing1.Sw))/2)*AC.Fuselage.fusWidth;
    AC.Wing2.Snet = AC.Wing2.Sw - ((AC.Wing2.RootChord+getChord(AC.Fuselage.fusWidth/2,AC.Wing2.WingSpan,AC.Wing2.TaperRatio,AC.Wing2.Sw))/2)*AC.Fuselage.fusWidth;

%Interference parameters
    KI_1 = (1 + 2.15 * AC.Fuselage.fusWidth / AC.Wing1.WingSpan) * AC.Wing1.Snet/AC.Wing1.Sw + ...
            pi * AC.Fuselage.fusWidth^2 / (2 * AC.Wing1.CL_alpha_w * AC.Wing1.Sw);
    KI_2 = (1 + 2.15 * AC.Fuselage.fusWidth / AC.Wing2.WingSpan) * AC.Wing2.Snet/AC.Wing2.Sw + ...
            pi * AC.Fuselage.fusWidth^2 / (2 * AC.Wing2.CL_alpha_w * AC.Wing2.Sw);
    KII_1 = (1 + 0.7 * AC.Fuselage.fusWidth / AC.Wing1.WingSpan ) * AC.Wing1.Snet/AC.Wing1.Sw;
    KII_2 = (1 + 0.7 * AC.Fuselage.fusWidth / AC.Wing2.WingSpan ) * AC.Wing2.Snet/AC.Wing2.Sw;
    deltazCL_1 = -0.1 * AC.Wing1.RootChord * AC.Fuselage.fusWidth / AC.Wing1.Sw;
    deltazCL_2 = -0.1 * AC.Wing2.RootChord * AC.Fuselage.fusWidth / AC.Wing2.Sw;

%Slope of the lift curve
    AC.Wing1.CL_alpha_wf = KI_1 * AC.Wing1.CL_alpha_w;
    AC.Wing2.CL_alpha_wf = KI_2 * AC.Wing2.CL_alpha_w;
    
%Downwash gradient
    deltaE_deltaAlpha = 1.75 * (AC.Wing1.CL_alpha_wf/(pi*AC.Wing1.AspectRatio*(AC.Wing1.TaperRatio*r)^0.25*(1+abs(m))));
    %1-deltaE_deltaAlpa = 0.9*(1-deltaE_deltaAlpha) <-- 10% reduction acording to Torenbeek E-10.1
    
%Lift coefficient taking into acount the fuselage interference
    AC.Wing1.CL_wf = AC.Wing1.CL_alpha_wf * deg2rad(AC.Fuselage.fuselage_AoA - twist_0lift_1*AC.Wing1.TipTwist + (KII_1/KI_1) * ...
                    (AC.Wing1.Incidence - AC.Wing1.Airfoil.alpha_zeroLift)) + deltazCL_1;
    AC.Wing2.CL_wf = AC.Wing2.CL_alpha_wf * deg2rad(AC.Fuselage.fuselage_AoA*0.9*(1-deltaE_deltaAlpha) - twist_0lift_2*AC.Wing2.TipTwist + ...
                    (KII_2/KI_2) * (AC.Wing2.Incidence - AC.Wing2.Airfoil.alpha_zeroLift)) + deltazCL_2 + AC.Wing2.deltaCLdeltaE;
    


                
                
%% NACELLES CONTRIBUTION
%Efect on the aerodynamic center for rear fuselage located nacelles
    deltaNac_1 = -2.5*AC.Engine.Number*AC.Engine.Diameter^2*AC.Engine.Length/(AC.Wing1.Sw*AC.Wing1.CMA*AC.Wing1.CL_alpha_wf);
    deltaNac_2 = -2.5*AC.Engine.Number*AC.Engine.Diameter^2*AC.Engine.Length/(AC.Wing2.Sw*AC.Wing2.CMA*AC.Wing2.CL_alpha_wf);
    
    AC.Wing1.x_ac_wf = AC.Wing1.x_ac_wf + deltaNac_1;
    AC.Wing2.x_ac_wf = AC.Wing2.x_ac_wf + deltaNac_2;
    
    
    
%% WING/FUSELAGE INTERFERENCE CORRECTION TO PITCHING MOMENT
%Forward shift of the aerodynamic center due to fuselage sections forward and after the wing
    lengthFusNose1 = AC.Wing1.Root_LE + AC.Fuselage.fusWidth/2*tand(AC.Wing1.Sweep_LE);
    lengthFusNose2 = AC.Wing2.Root_LE + AC.Fuselage.fusWidth/2*tand(AC.Wing2.Sweep_LE);
    deltaF1ac_1 = - (1.8/AC.Wing1.CL_alpha_wf) * (AC.Fuselage.fusWidth*AC.Fuselage.fusHeight*lengthFusNose1/AC.Wing1.Sw);
    deltaF1ac_2 = - (1.8/AC.Wing2.CL_alpha_wf) * (AC.Fuselage.fusWidth*AC.Fuselage.fusHeight*lengthFusNose2/AC.Wing1.Sw);

%Aerodynamic center correction due to lift loss in the region where the wing/fuselage lift carry-over is concentrated
    deltaF2ac_1 = 0.273/(1+AC.Wing1.TaperRatio) * AC.Fuselage.fusWidth * AC.Wing1.CMG/AC.Wing1.CMA *...
                 (AC.Wing1.WingSpan-AC.Fuselage.fusWidth)/(AC.Wing1.WingSpan+2.15*AC.Fuselage.fusWidth) * tand(AC.Wing1.Sweep_14);
    deltaF2ac_2 = 0.273/(1+AC.Wing2.TaperRatio) * AC.Fuselage.fusWidth * AC.Wing2.CMG/AC.Wing2.CMA *...
                 (AC.Wing2.WingSpan-AC.Fuselage.fusWidth)/(AC.Wing2.WingSpan+2.15*AC.Fuselage.fusWidth) * tand(AC.Wing2.Sweep_14);

%Aerodynamic center    
    AC.Wing1.x_ac_wf = AC.Wing1.x_ac_w + deltaF1ac_1 + deltaF2ac_1; 
    AC.Wing2.x_ac_wf = AC.Wing2.x_ac_w + deltaF1ac_2 + deltaF2ac_2; 

%Fuselage contribution to the pitching moment (from Munk's Theory, Ref: E-33)
    CL0_1 = AC.Wing1.CL_alpha_wf * deg2rad(0 - twist_0lift_1*AC.Wing1.TipTwist + (KII_1/KI_1) * ...
           (AC.Wing1.Incidence - AC.Wing1.Airfoil.alpha_zeroLift)) + deltazCL_1;
    CL0_2 = AC.Wing2.CL_alpha_wf * deg2rad(0 - twist_0lift_2*AC.Wing2.TipTwist + (KII_2/KI_2) * ...
           (AC.Wing2.Incidence - AC.Wing2.Airfoil.alpha_zeroLift)) + deltazCL_2;

    deltaFCmac_1 = -1.8 * (1-2.5*AC.Fuselage.fusWidth/AC.Fuselage.fusLength) * (pi*AC.Fuselage.fusWidth*AC.Fuselage.fusHeight*...
                    AC.Fuselage.fusLength/(4*AC.Wing1.Sw*AC.Wing1.CMA)) * (CL0_1/AC.Wing1.CL_alpha_wf);
    deltaFCmac_2 = -1.8 * (1-2.5*AC.Fuselage.fusWidth/AC.Fuselage.fusLength) * (pi*AC.Fuselage.fusWidth*AC.Fuselage.fusHeight*...
                    AC.Fuselage.fusLength/(4*AC.Wing2.Sw*AC.Wing2.CMA)) * (CL0_2/AC.Wing2.CL_alpha_wf);
                
%Pitching moment coefficient in the aerodynamic center
    AC.Wing1.Cm_ac_wf = AC.Wing1.Cm_ac_w + deltaFCmac_1;
    AC.Wing2.Cm_ac_wf = AC.Wing2.Cm_ac_w + deltaFCmac_2;
    
%Pitching moment coefficient with fuselage interference
    AC.Wing1.Cm_wf = AC.Wing1.Cm_ac_wf + AC.Wing1.CL_wf*((AC.Weight.x_cg-AC.Wing1.x_ac_wf)/AC.Wing1.CMA);
    AC.Wing2.Cm_wf = AC.Wing2.Cm_ac_wf + AC.Wing2.CL_wf*((AC.Weight.x_cg-AC.Wing2.x_ac_wf)/AC.Wing2.CMA);
    
    
clear lengthFusNose1 lengthFusNose2 deltaF1ac_1 deltaF1ac_2 deltaF2ac_1 deltaF2ac_2 designWeight CL0 CL0_1 CL0_2 deltaFCmac_1 deltaFCmac_2

    
    
    
%% LIFT OF THE COMPLETE AIRCRAFT
    %Coeficiente de sustentación y de momentos del avión completo
    AC.Wing.CL_wf = Parameters.q1_qinf * AC.Wing1.Sw / AC.Wing.Sw * AC.Wing1.CL_wf + Parameters.q2_qinf * AC.Wing2.Sw / AC.Wing.Sw * AC.Wing2.CL_wf;
    AC.Wing.Cm_wf = Parameters.q1_qinf * AC.Wing1.Sw / AC.Wing.Sw * AC.Wing1.CMA / AC.Wing.CMA * AC.Wing1.Cm_ac_wf + ...
                    Parameters.q2_qinf * AC.Wing2.Sw / AC.Wing.Sw * AC.Wing2.CMA / AC.Wing.CMA * AC.Wing2.Cm_ac_wf + ...
                    Parameters.q1_qinf * AC.Wing1.Sw / AC.Wing.Sw * AC.Wing1.CL_wf * (AC.Weight.x_cg - AC.Wing1.x_ac_wf) + ...
                    Parameters.q2_qinf * AC.Wing2.Sw / AC.Wing.Sw * AC.Wing2.CL_wf * (AC.Weight.x_cg - AC.Wing2.x_ac_wf);



%% STORE TOTAL LIFT DISTRIBUTION
%Coefficient of Aditional lift
    AC.Wing1.cla = AC.Wing1.CL_wf.*AC.Wing1.CMG.*La1./AC.Wing1.c;
    AC.Wing2.cla = AC.Wing2.CL_wf.*AC.Wing2.CMG.*La2./AC.Wing2.c;
    
%Coefficient of Basic Lift
    AC.Wing1.clb = Lb1.*deg2rad(AC.Wing1.TipTwist).*AC.Wing1.Airfoil.Cl_alpha.*AC.Wing1.CMG./(AC.Wing1.E.*AC.Wing1.c);
    AC.Wing2.clb = Lb2.*deg2rad(AC.Wing2.TipTwist).*AC.Wing2.Airfoil.Cl_alpha.*AC.Wing2.CMG./(AC.Wing2.E.*AC.Wing2.c);

%Total Lift Distribution
    AC.Wing1.cl = AC.Wing1.cla + AC.Wing1.clb;
    AC.Wing2.cl = AC.Wing2.cla + AC.Wing2.clb;






%% PLOT WING LAYOUT
if DP.ShowAircraftLayout
%Import fuselage layout and tail cone layout from data file
sr=which(mfilename);
i=max(strfind(lower(sr),lower('MTORRES')))+6;
if i>0
  sr=sr(1:i);
else
  error('Cannot locate MTORRES directory. You must add the path to the MTorres directory.')
end
FuselageFile = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'fuselage.dat');
AfterbodyFile = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'tailCoordinates.dat');
[Xfus,Yfus]  = importFuselage(FuselageFile);
[Xtail,Ytail]  = importFuselage(AfterbodyFile);
clear sr i FuselageFile AfterbodyFile


%Create figure and plotting  
figure()
    hold on
    axis equal
    %COCKPIT
        plot(Xfus, Yfus,'k')
        plot(Xfus,-Yfus,'k')
    %CABIN
        plot([AC.Fuselage.ln,AC.Fuselage.ln+AC.Fuselage.cabLength],[ AC.Fuselage.fusWidth/2, AC.Fuselage.fusWidth/2],'k')
        plot([AC.Fuselage.ln,AC.Fuselage.ln+AC.Fuselage.cabLength],[-AC.Fuselage.fusWidth/2,-AC.Fuselage.fusWidth/2],'k')
    %TAIL CONE
        plot(AC.Fuselage.fusLength-AC.Fuselage.la+Xtail, Ytail,'k')
        plot(AC.Fuselage.fusLength-AC.Fuselage.la+Xtail,-Ytail,'k')
    %WING1
        %root chord
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.RootChord],[ 0, 0],'r')
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.RootChord],[-0,-0],'r')
        %tip chord
            plot([AC.Wing1.Root_LE+AC.Wing1.TipSweep,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ AC.Wing1.WingSpan/2, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.Root_LE+AC.Wing1.TipSweep,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-AC.Wing1.WingSpan/2,-AC.Wing1.WingSpan/2],'r')
        %leading edge
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.TipSweep],[ 0, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.TipSweep],[-0,-AC.Wing1.WingSpan/2],'r')
        %trailing edge
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ 0, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-0,-AC.Wing1.WingSpan/2],'r')
        %position c/4
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord/4,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[ 0, AC.Wing1.WingSpan/2],'g')
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord/4,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[-0,-AC.Wing1.WingSpan/2],'g')
        %MAC
            plot([AC.Wing1.CMA_LE,AC.Wing1.CMA_LE+AC.Wing1.CMA],[ AC.Wing1.CMA_b, AC.Wing1.CMA_b],'m')
            plot([AC.Wing1.CMA_LE,AC.Wing1.CMA_LE+AC.Wing1.CMA],[-AC.Wing1.CMA_b,-AC.Wing1.CMA_b],'m')
    %WING2
        %root chord
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.RootChord],[ 0, 0],'b')
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.RootChord],[-0,-0],'b')
        %tip chord
            plot([AC.Wing2.Root_LE+AC.Wing2.TipSweep,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[ AC.Wing2.WingSpan/2, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.Root_LE+AC.Wing2.TipSweep,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[-AC.Wing2.WingSpan/2,-AC.Wing2.WingSpan/2],'b')
        %leading edge
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.TipSweep],[ 0, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.TipSweep],[-0,-AC.Wing2.WingSpan/2],'b')
        %trailing edge
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[ 0, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[-0,-AC.Wing2.WingSpan/2],'b')
        %position c/4
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord/4,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[ 0, AC.Wing2.WingSpan/2],'g')
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord/4,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[-0,-AC.Wing2.WingSpan/2],'g')
        %MAC
            plot([AC.Wing2.CMA_LE,AC.Wing2.CMA_LE+AC.Wing2.CMA],[ AC.Wing2.CMA_b, AC.Wing2.CMA_b],'m')
            plot([AC.Wing2.CMA_LE,AC.Wing2.CMA_LE+AC.Wing2.CMA],[-AC.Wing2.CMA_b,-AC.Wing2.CMA_b],'m')
    %WEIGHT
        %Center of gravity
            plot(AC.Weight.x_cg,AC.Weight.y_cg,'*')
    clear Xfus Yfus
end     
    
    
clear a T P rho nu EffectiveSweep1 EffectiveSweep2 CruiseMach1 CruiseMach2 Xtail Ytail
clear f_00deg f_30deg f_45deg f_60deg index DiederichCoordinate1 DiederichCoordinate2
clear C1x C2x C3x C4x C1y C2y C3y C4y f0x f0y f30x f30y f_low f_hight Twist1 Twist2
clear C1_1 C1_2 C2_1 C2_2 C3_1 C3_2 f1 f2 C4_1 C4_2 CL_wing1 CL_wing2 Clb1 Clb2
clear La1 La2 Lb1 Lb2 twist_0lift_1 twist_0lift_2 Cl1 Cl2 y1 y2 Cla1 Cla2 Beta1 Beta2
clear Cm_ac_basic1 Cm_ac_basic2 deltaEpsilonCm_ac1 deltaEpsilonCm_ac2 Cm_ac_w1 Cm_ac_w2
clear KI_1 KI_2 KII_1 KII_2 deltazCL_1 deltazCL_2 deltaNac_1 deltaNac_2 deltaE_deltaAlpha m r
    
    
%% USEFUL FUNCTIONS
function Error = getDivergenceMach(DivergenceMach,t_c_airfoil,Sweep,CLwing,AirfoilType)
    switch lower(AirfoilType)
        case 'conventional'
            Mstar = 1.00 - 0.25*CLwing*(cosd(Sweep))^-2;
        case 'high-speed'
            Mstar = 1.05 - 0.25*CLwing*(cosd(Sweep))^-2;
        case 'supercritical'
            Mstar = 1.15 - 0.25*CLwing*(cosd(Sweep))^-2;
    end
    maximum_t_c_wing = (0.3/DivergenceMach) * ...
                       ((1/(DivergenceMach*cosd(Sweep))-DivergenceMach*cosd(Sweep))^(1/3)) * ...
                       ((1-((5+(DivergenceMach*cosd(Sweep))^2)/(5+Mstar^2))^3.5)^(2/3));
    Error = t_c_airfoil - maximum_t_c_wing / cosd(Sweep);
end

function chord = getChord(y, SpanWidth, TaperRatio, Sw)
    chord = (2*Sw/((1+TaperRatio)*SpanWidth)).*(1-(2*(1-TaperRatio)/SpanWidth).*abs(y));
end


