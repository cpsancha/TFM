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

%Dihedral
    AC.Wing1.Dihedral = DP.Dihedral;
    AC.Wing2.Dihedral = DP.Dihedral;
 
%Longitudinal Position of the leading edge at root
    AC.Wing1.LongPos = DP.Wing1LongPos;
    AC.Wing2.LongPos = AC.Wing1.LongPos+AC.Wing1.RootChord+DP.Stagger;
    
%Y_CMA
    AC.Wing1.CMA_b = (AC.Wing1.WingSpan/6)*((1+2*AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio));
    AC.Wing2.CMA_b = (AC.Wing2.WingSpan/6)*((1+2*AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio));

%X_CMA_LE
    AC.Wing1.CMA_LE = (AC.Wing1.WingSpan/6)*((1+2*AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*tand(AC.Wing1.Sweep_LE);
    AC.Wing2.CMA_LE = (AC.Wing2.WingSpan/6)*((1+2*AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*tand(AC.Wing2.Sweep_LE);

%X_Tip_LE
    AC.Wing1.TipSweep = (AC.Wing1.WingSpan/2)*tand(AC.Wing1.Sweep_LE);
    AC.Wing2.TipSweep = (AC.Wing2.WingSpan/2)*tand(AC.Wing2.Sweep_LE);

%Reynolds
        [rho,a,T,P,nu,~] = atmos(DP.CruiseAltitude);
        CruiseMach = DP.CruiseSpeed/a;
        Beta1       = sqrt(1-(CruiseMach)^2); %Prandtl compresibility correction
        Beta2       = sqrt(1-(CruiseMach)^2); %Prandtl compresibility correction
        Reynolds1  = DP.CruiseSpeed * AC.Wing1.CMA / nu;
        Reynolds2  = DP.CruiseSpeed * AC.Wing2.CMA / nu;
    
%Wing CL
    AC.Wing1.CLdesign = AC.Wing1.Airfoil.Cldesign*(cosd(AC.Wing1.Sweep_14))^2;
    AC.Wing2.CLdesign = AC.Wing2.Airfoil.Cldesign*(cosd(AC.Wing2.Sweep_14))^2;

%Wing t_c
    AC.Wing1.t_c = AC.Wing1.Airfoil.t_c*cosd(AC.Wing1.Sweep_14);
    AC.Wing2.t_c = AC.Wing2.Airfoil.t_c*cosd(AC.Wing2.Sweep_14);
    
%Mach of divergence acording to maximum thickness, CL, and sweep angle
    options = optimoptions('fsolve',...
                           'StepTolerance',1e-9,...
                           'Display','none');
    [AC.Wing1.MachDiv,~,exitflag1,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, AC.Wing1.Airfoil.t_c, AC.Wing1.Sweep_14,...
                                               AC.Wing1.CLdesign,'Supercritical'),0.8,options);
    [AC.Wing2.MachDiv,~,exitflag2,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, AC.Wing2.Airfoil.t_c, AC.Wing2.Sweep_14,...
                                               AC.Wing2.CLdesign,'Supercritical'),0.8,options);
    if ~isequal(exitflag1,1)
        disp('El solver del mach de divergencia del ala 1 no ha logrado converger correctamente. Se debería revisar el resultado.')
        pause
    elseif ~isequal(exitflag2,1)
        disp('El solver del mach de divergencia del ala 2 no ha logrado converger correctamente. Se debería revisar el resultado.')
        pause    
    else
        clear exitflag1 exitflag2 options
    end

    
    
    
%% LIFTING PROPERTIES OF AIRFOIL SECTIONS
%TBD... interpolación en Mach y Re
load('.\Temporary Stuff\M05_Re04.mat')
[fit,R2] = polyfitR2(deg2rad(alpha_Cn_M05_Re04(:,1)),alpha_Cn_M05_Re04(:,2),1); %#ok<ASGLU>
AC.Wing1.Airfoil.Cl_alpha = fit(1); %<--Cl_alpha a Mach = Mach_infty * cos(Sweep_14)
AC.Wing1.Airfoil.alpha_zeroLift = rad2deg(-fit(2)/fit(1));
AC.Wing2.Airfoil.Cl_alpha = fit(1);
AC.Wing2.Airfoil.alpha_zeroLift = rad2deg(-fit(2)/fit(1));
clear fit R2 alpha_Cn_M05_Re04 Cd_Cn_M05_Re04 Cm_Cn_M05_Re04



%% WING LIFT AND LIFT DISTRIBUTIONS
%LIFT-CURVE SLOPE
%Torenbeek Method (Torenbeek E-4.1.b pag 473)
    EffectiveSweep1 = atand(tand(AC.Wing1.Sweep_14)/Beta1);
    EffectiveSweep2 = atand(tand(AC.Wing2.Sweep_14)/Beta2);
    k1 = (AC.Wing1.Airfoil.Cl_alpha)/(2*pi);
    k2 = (AC.Wing2.Airfoil.Cl_alpha)/(2*pi);
    CL_alpha_Torenbeek = (2*pi/Beta1)/((2/(Beta1*AC.Wing1.AspectRatio))+sqrt((1/((k1*cosd(EffectiveSweep1))^2))+(2/(Beta1*AC.Wing1.AspectRatio))^2));
%DATCOM Method - Straight-Tapered Wings - Method 1 - pag 4.1.3.2-3
    CL_alpha_DATCOM    = (2*pi*AC.Wing1.AspectRatio)/(2+sqrt((AC.Wing1.AspectRatio*Beta1/k1)^2*(1+(tand(AC.Wing1.Sweep_12)/Beta1)^2)+4));
%Polhamus Method - NACA Technical Note 3911
    CL_alpha_Polhamus  = (AC.Wing1.Airfoil.Cl_alpha*AC.Wing1.AspectRatio) / ((AC.Wing1.Airfoil.Cl_alpha/pi) + ...
                         sqrt((AC.Wing1.AspectRatio/cosd(AC.Wing1.Sweep_12))^2+(AC.Wing1.Airfoil.Cl_alpha/pi)^2 - ...
                         (AC.Wing1.AspectRatio*CruiseMach)^2));
    if abs(CL_alpha_Torenbeek-CL_alpha_DATCOM) > 1e-6
%         warning('La pendiente de la curva de sustentacion por el metodo de Torenbeek no es igual que por el metodo de las DATCOM. Elige cual se usa.')
    elseif abs(CL_alpha_Torenbeek-CL_alpha_Polhamus) > 1e-6
        warning('La pendiente de la curva de sustentacion por el metodo de Torenbeek no es igual que por el metodo de Polhamus. Elige cual se usa.')
    end
%Uses Torenbeek Method, because yes... to be decided when loaded airfoil data
    AC.Wing1.CL_alpha  = (2*pi*AC.Wing1.AspectRatio)/(2+sqrt((AC.Wing1.AspectRatio*Beta1/k1)^2*(1+(tand(AC.Wing1.Sweep_12)/Beta1)^2)+4));
    AC.Wing2.CL_alpha  = (2*pi*AC.Wing2.AspectRatio)/(2+sqrt((AC.Wing2.AspectRatio*Beta2/k2)^2*(1+(tand(AC.Wing2.Sweep_12)/Beta2)^2)+4));
    clear k1 k2 CL_alpha_Torenbeek CL_alpha_DATCOM CL_alpha_Polhamus
    
%SPANWISE LIFT DISTRIBUTION
%Create spanwise coordinates
    y1   = linspace(0,AC.Wing1.WingSpan/2,100);
    y2   = linspace(0,AC.Wing2.WingSpan/2,100);
%Adimensional spanwise coordinate
    eta1 = y1./(AC.Wing1.WingSpan/2);
    eta2 = y2./(AC.Wing2.WingSpan/2);
%Load graph from NACA Technical Note 2751
    load('.\Temporary Stuff\Lift_Distribution.mat')
    if EffectiveSweep1 < 30
        f_low   = interp1(f_00deg(:,1),f_00deg(:,2),eta1);
        f_hight = interp1(f_30deg(:,1),f_30deg(:,2),eta1);
        f1 = f_low+(f_hight-f_low)./30.*EffectiveSweep1;
    elseif EffectiveSweep1 < 45
        f_low   = interp1(f_30deg(:,1),f_30deg(:,2),eta1);
        f_hight = interp1(f_45deg(:,1),f_45deg(:,2),eta1);
        f1 = f_low+(f_hight-f_low)./15.*(EffectiveSweep1-30);
    elseif EffectiveSweep1 < 60
        f_low   = interp1(f_45deg(:,1),f_45deg(:,2),eta1);
        f_hight = interp1(f_60deg(:,1),f_60deg(:,2),eta1);
        f1 = f_low+(f_hight-f_low)./15.*(EffectiveSweep1-45);
    else
        warning('Incorrect value for Effective Sweep of wing 1.')
    end

    if EffectiveSweep2 < 30
        f_low   = interp1(f_00deg(:,1),f_00deg(:,2),eta1);
        f_hight = interp1(f_30deg(:,1),f_30deg(:,2),eta1);
        f2 = f_low+(f_hight-f_low)./30.*EffectiveSweep2;
    elseif EffectiveSweep2 < 45
        f_low   = interp1(f_30deg(:,1),f_30deg(:,2),eta1);
        f_hight = interp1(f_45deg(:,1),f_45deg(:,2),eta1);
        f2 = f_low+(f_hight-f_low)./15.*(EffectiveSweep2-30);
    elseif EffectiveSweep2 < 60
        f_low   = interp1(f_45deg(:,1),f_45deg(:,2),eta1);
        f_hight = interp1(f_60deg(:,1),f_60deg(:,2),eta1);
        f2 = f_low+(f_hight-f_low)./15.*(EffectiveSweep2-45);
    else
        warning('Incorrect value for Effective Sweep of wing 2.')
    end

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
    c1 = getChord(y1, AC.Wing1.WingSpan, AC.Wing1.TaperRatio, AC.Wing1.Sw);
    c2 = getChord(y2, AC.Wing2.WingSpan, AC.Wing2.TaperRatio, AC.Wing2.Sw);
    


    
%Aditional lift distribution calculation La
% La1 = C1_1.*c./AC.Wing1.CMG+(C2+C3).*4./pi.*sqrt(1-eta.^2)+C3.*f; %<--Está mal, no? sobra el último término, creo que has mezclado las fórmulas E13 y E14
% 
% %Basic lift distribution Lb. Si no tienes torsion no hace falta calcularla
% %porque luego se multiplica por la torsion en la punta pa k kieres saber
% %eso jaja salu2
% epsilon_epsilont = eta; %Torsión lineal desde 0 hasta la torsion en punta epsilon=epsilont*eta;
% alpha_0_1 = -trapz(eta,epsilon_epsilont.*La1); %Eq E-16
% lambda_beta = atan((tan(AC.Wing1.Sweep_12)/ME.Cruise.beta)); %Valor nulo, flecha 1/2 nula
% Lb = La1.*C4.*cos(lambda_beta).*(epsilon_epsilont+alpha_0_1).*ME.Cruise.beta*E;
% 
% % figure()
% % hold on
% % plot(eta, La)
% % plot(eta,Lb)
% 
% %cl distribution
% % CL = 2*AC.Weight.MTOW*CST.GravitySI/(ME.Cruise.Density*ME.Cruise.Speed^2*AC.Wing.Sw); %Coeficiente de sustentacion orientativo en crucero
% CL = 1;
% cl = (La1.*CL+AC.Wing1.Torsion.*cl_alpha.*Lb./E)./c.*AC.Wing1.CMG;
% cla = CL.*AC.Wing1.CMG.*La1./c;
% clb = AC.Wing1.Torsion .* cl_alpha .* AC.Wing1.CMG .* Lb ./ (E.*c);
    
    

clear f_00deg f_30deg f_45deg f_60deg DiederichCoordinate1 DiederichCoordinate2
clear C1x C2x C3x C4x C1y C2y C3y C4y f0x f0y f30x f30y f_low f_hight



%% PLOT WING LAYOUT
%Import fuselage layout from data file
sr=which(mfilename);
i=max(strfind(lower(sr),lower('MTORRES')))+6;
if i>0
  sr=sr(1:i);
else
  error('Cannot locate MTORRES directory. You must add the path to the MTorres directory.')
end
FuselageFile = fullfile(sr,'Matlab Code',filesep,'Temporary Stuff',filesep,'fuselage.dat');
[Xfus,Yfus]  = importFuselage(FuselageFile);
clear sr i FuselageFile
   

%Create figure and plotting  
figure()
    hold on
    axis equal
    %COCKPIT
        plot(Xfus, Yfus,'k')
        plot(Xfus,-Yfus,'k')
    %CABIN
        plot([Xfus(end),DP.fusLength],[ Yfus(end), Yfus(end)],'k')
        plot([Xfus(end),DP.fusLength],[-Yfus(end),-Yfus(end)],'k')
    %WING1
        %root chord
            plot([AC.Wing1.LongPos,AC.Wing1.LongPos+AC.Wing1.RootChord],[ 0, 0],'r')
            plot([AC.Wing1.LongPos,AC.Wing1.LongPos+AC.Wing1.RootChord],[-0,-0],'r')
        %tip chord
            plot([AC.Wing1.LongPos+AC.Wing1.TipSweep,AC.Wing1.LongPos+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ AC.Wing1.WingSpan/2, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.LongPos+AC.Wing1.TipSweep,AC.Wing1.LongPos+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-AC.Wing1.WingSpan/2,-AC.Wing1.WingSpan/2],'r')
        %leading edge
            plot([AC.Wing1.LongPos,AC.Wing1.LongPos+AC.Wing1.TipSweep],[ 0, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.LongPos,AC.Wing1.LongPos+AC.Wing1.TipSweep],[-0,-AC.Wing1.WingSpan/2],'r')
        %trailing edge
            plot([AC.Wing1.LongPos+AC.Wing1.RootChord,AC.Wing1.LongPos+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ 0, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.LongPos+AC.Wing1.RootChord,AC.Wing1.LongPos+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-0,-AC.Wing1.WingSpan/2],'r')
        %position c/4
            plot([AC.Wing1.LongPos+AC.Wing1.RootChord/4,AC.Wing1.LongPos+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[ 0, AC.Wing1.WingSpan/2],'g')
            plot([AC.Wing1.LongPos+AC.Wing1.RootChord/4,AC.Wing1.LongPos+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[-0,-AC.Wing1.WingSpan/2],'g')
    %WING2
        %root chord
            plot([AC.Wing2.LongPos,AC.Wing2.LongPos+AC.Wing2.RootChord],[ 0, 0],'b')
            plot([AC.Wing2.LongPos,AC.Wing2.LongPos+AC.Wing2.RootChord],[-0,-0],'b')
        %tip chord
            plot([AC.Wing2.LongPos+AC.Wing2.TipSweep,AC.Wing2.LongPos+AC.Wing2.TipSweep+AC.Wing2.TipChord],[ AC.Wing2.WingSpan/2, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.LongPos+AC.Wing2.TipSweep,AC.Wing2.LongPos+AC.Wing2.TipSweep+AC.Wing2.TipChord],[-AC.Wing2.WingSpan/2,-AC.Wing2.WingSpan/2],'b')
        %leading edge
            plot([AC.Wing2.LongPos,AC.Wing2.LongPos+AC.Wing2.TipSweep],[ 0, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.LongPos,AC.Wing2.LongPos+AC.Wing2.TipSweep],[-0,-AC.Wing2.WingSpan/2],'b')
        %trailing edge
            plot([AC.Wing2.LongPos+AC.Wing2.RootChord,AC.Wing2.LongPos+AC.Wing2.TipSweep+AC.Wing2.TipChord],[ 0, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.LongPos+AC.Wing2.RootChord,AC.Wing2.LongPos+AC.Wing2.TipSweep+AC.Wing2.TipChord],[-0,-AC.Wing2.WingSpan/2],'b')
        %position c/4
            plot([AC.Wing2.LongPos+AC.Wing2.RootChord/4,AC.Wing2.LongPos+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[ 0, AC.Wing2.WingSpan/2],'g')
            plot([AC.Wing2.LongPos+AC.Wing2.RootChord/4,AC.Wing2.LongPos+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[-0,-AC.Wing2.WingSpan/2],'g')
    clear Xfus Yfus
     
    
    
clear T P rho nu Reynolds1 Reynolds2 CruiseMach EffectiveSweep1 EffectiveSweep2 Beta1 Beta2 y1 y2 eta1 eta2
    
    
    
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