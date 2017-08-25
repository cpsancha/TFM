%% AIRCRAFT POLAR --> Torenbeek Apendice F

%Define range of weights to calculate polar
AircraftWeight = linspace(AC.Weight.EW, AC.Weight.MTOW, 5);

options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
for i=1:length(AircraftWeight)
    %Trim aircraft for that weight
    [~,~,exitflag,~] = fsolve(@(X)trimAircraft(X, AircraftWeight(i), AC, ME, DP, Parameters, CST, CF),[DP.fuselage_AoA, 0],options);
    if ~isequal(exitflag,1)
        error('El solver del trimado no ha logrado converger correctamente. Se debería revisar el resultado.')
    else
        clear exitflag
    end
    
    %Store input polar data
    Polar.Weight(i)       = AircraftWeight(i);
    Polar.Fuselage_AoA(i) = AC.Fuselage.fuselage_AoA;
    Polar.Wing1.CL(i)     = AC.Wing1.CL_wf;
    Polar.Wing1.Cma(i)    = AC.Wing1.Cm_ac_wf;
    Polar.Wing2.CL(i)     = AC.Wing2.CL_wf;
    Polar.Wing2.Cma(i)    = AC.Wing2.Cm_ac_wf;
    Polar.Wing2.deltaCLdeltaE(i) = AC.Wing2.deltaCLdeltaE;
    Polar.CL(i) = Parameters.q1_qinf * AC.Wing1.Sw/AC.Wing.Sw * AC.Wing1.CL_wf + Parameters.q2_qinf * AC.Wing2.Sw/AC.Wing.Sw * AC.Wing2.CL_wf;
    %Get polar and store outputs
    [Polar.CD(i), Polar.D{i}, Polar.DS{i}, ME] = getPolar(AC, ME, DP, Parameters, CST, CF, false);
end
%Create polar regression
[Polar.PolarFit,Polar.Rsquared] = polyfitR2(Polar.CL,Polar.CD,2);



if DP.ShowReportFigures
%Show Aircraft Polar
figure()
    hold on
    plot(Polar.CL,Polar.CD)
    txt=['$\ \leftarrow\ \ C_{D}=',num2str(Polar.PolarFit(3)),'',num2str(Polar.PolarFit(2)),...
         'C_{L}+',num2str(Polar.PolarFit(1)),'C_{L}^2','$'];
    text(Polar.CL(round(0.1*length(Polar.CL))),Polar.CD(round(0.1*length(Polar.CL))),txt,'HorizontalAlignment','left','Interpreter','Latex')
    txt=['$\ \ \ \ \ \ \ R^{2}=',num2str(Polar.Rsquared),'$'];
    text(Polar.CL(round(0.1*length(Polar.CL))),Polar.CD(round(0.1*length(Polar.CL)))-0.002,txt,'HorizontalAlignment','left','Interpreter','Latex')
    plot(linspace(min(Polar.CL),max(Polar.CL),50), polyval(Parameters.Polar.LongRangeCruise,linspace(min(Polar.CL),max(Polar.CL),50)))
    txt=['$C_{D}=',num2str(Parameters.Polar.LongRangeCruise(3)),'+',num2str(Parameters.Polar.LongRangeCruise(1)),'C_{L}^2','\ \ \rightarrow\ \ \ \ $'];
    text(min(Polar.CL)+0.6*(max(Polar.CL)-min(Polar.CL)),polyval(Parameters.Polar.LongRangeCruise,min(Polar.CL)+0.6*(max(Polar.CL)-min(Polar.CL))),txt,'HorizontalAlignment','right','Interpreter','Latex')
    legend('Calculated Polar','Design Polar','Location','southeast')
    legend('boxoff')
    title('Aircraft polar for the valid range of weights','interpreter','latex')
    xlabel('C_L')
    ylabel('C_D')
    xlim([min(Polar.CL)-0.005,max(Polar.CL)+0.005])
    saveFigure(ME.FiguresFolder,'AircraftPolar')
    clear txt
end



[Polar.maxEfficiency, Polar.maxEfficiencyIndex] = max(Polar.CL./Polar.CD);
[~,~,exitflag,~] = fsolve(@(X)trimAircraft(X, AircraftWeight(Polar.maxEfficiencyIndex), AC, ME, DP, Parameters, CST, CF),[DP.fuselage_AoA, 0],options);
if ~isequal(exitflag,1)
    error('El solver del trimado no ha logrado converger correctamente. Se debería revisar el resultado.')
else
    clear exitflag
end
getPolar(AC, ME, DP, Parameters, CST, CF, DP.ShowReportFigures);


% AC.Fuselage.fuselage_AoA = 0;
% AC.Wing2.deltaCLdeltaE = 0;


clear i X options AircraftWeight 

%% USEFUL FUNCTIONS
function [Error, ME] = trimAircraft(X, AircraftWeight, AC, ME, DP, Parameters, CST, CF)
  %Parse inputs
    AC.Fuselage.fuselage_AoA = X(1);
    AC.Wing2.deltaCLdeltaE   = X(2);
    
  %Run wing's script
    ME = wingsDesign(AC, ME, DP, Parameters, CST, CF);
    
  %Necessary calculation
	[rho,~,~,~,~,~] = atmos(DP.CruiseAltitude);
 	CL0             = AircraftWeight*CST.GravitySI / (0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw); 
    
  %Parse outputs
    Error(1) =  AC.Wing.CL_wf - CL0;
    Error(2) =  AC.Wing.Cm_wf;   
end

function [CD, D, DS, ME] = getPolar(AC, ME, DP, Parameters, CST, CF, dragPlotFlag) %#ok<INUSL>
%METHOD A --> A. Garner [Ref. F-39]
%Spanwise center of pressure    
    eta_cp_1 = trapz(AC.Wing1.eta,AC.Wing1.eta.*AC.Wing1.cla.*AC.Wing1.c./AC.Wing1.CL_wf./AC.Wing1.CMG);
    eta_cp_2 = trapz(AC.Wing2.eta,AC.Wing2.eta.*AC.Wing2.cla.*AC.Wing2.c./AC.Wing2.CL_wf./AC.Wing2.CMG);
    
%Increment of plane wing vortex-induced drag coefficient due to additional lift
	delta_1 = 46.264.*(eta_cp_1 - 4/3/pi).^2;
    delta_2 = 46.264.*(eta_cp_2 - 4/3/pi).^2;
    
%Vortex-Induced Drag of untwisted plane wings
    D.CDv_1 = (1+delta_1)*AC.Wing1.CL_wf^2/(pi*AC.Wing1.AspectRatio);
    D.CDv_2 = (1+delta_2)*AC.Wing2.CL_wf^2/(pi*AC.Wing2.AspectRatio);

%INDUCED DRAG DUE TO TWIST CORRECTION
    %Valores de v y w de las figuras 11 y 12 de la Ref. F-48, para A=8, \lambda=0.6
    v = 0.0015;
    w = 0.0039;
    if ~isequal(DP.AspectRatio,8) || ~isequal(DP.TaperRatio,0.6)
        if ~ismember('AC:v_w_parameters_modification',ME.errorList)
            ME.errorList{end+1} = 'AC:v_w_parameters_modification';
            warning('AC:v_w_parameters_modification',['Se ha modificado el alargamiento o el estrechamiento, es necesario actualizar los'...
                    ' valores de los parámetros v y w que determinan la resistencia inducida debida a la torsion.'])
        end
    end
    D.deltaEpsilonCDv_1 = AC.Wing1.CL_wf*(deg2rad(AC.Wing1.TipTwist)*AC.Wing1.CL_alpha_wf/AC.Wing1.E)*v+(deg2rad(AC.Wing1.TipTwist)*AC.Wing1.CL_alpha_wf/AC.Wing1.E)^2*w;
    D.deltaEpsilonCDv_2 = AC.Wing2.CL_wf*(deg2rad(AC.Wing2.TipTwist)*AC.Wing2.CL_alpha_wf/AC.Wing2.E)*v+(deg2rad(AC.Wing2.TipTwist)*AC.Wing2.CL_alpha_wf/AC.Wing2.E)^2*w;
    
%INDUCED DRAG DUE TO WING TIP CORRECTION
    %No se corrije porq no se sabe como
    
%INDUCED DRAG DUE TO FUSELAGE LIFT CORRECTION
    DS.fuselageInduced = 0.15*AC.Fuselage.fuselage_AoA^2*AC.Fuselage.Volume^(2/3);

%NACELLES CONTRIBUTION
    %No sense to estimate, interferences will be much larger
 
%PLOT EACH CONTRIBUTION
if dragPlotFlag
    X = [Parameters.q1_qinf*AC.Wing1.Sw*D.CDv_1,Parameters.q1_qinf*AC.Wing1.Sw*D.deltaEpsilonCDv_1,...
         Parameters.q2_qinf*AC.Wing2.Sw*D.CDv_2,Parameters.q2_qinf*AC.Wing2.Sw*D.deltaEpsilonCDv_2,...
         DS.fuselageInduced];
    labels =  {'Sustentación ala delantera: ';'Torsión ala delantera: ';...
               'Sustentación ala trasera: ';'Torsión ala trasera: ';...
               'Sustentación del fuselaje: '};
    fInduced = drawCustomPieChart(X,labels,false,false);
    title('Contribuciones a la resistencia inducida','interpreter','latex')
    saveFigure(ME.FiguresFolder,'InducedDragPie') 
end
        
    
%%PROFILE DRAF OF SMOOTH, ISOLATED MAJOR COMPONENTS --> FLAT PLATE ANALOGY
%WING SECTIONS
    [rho,~,~,~,nu,~] = atmos(DP.CruiseAltitude);
    
%Minimum Cd profile
    Phi_w_1 = 2.7*AC.Wing1.Airfoil.t_c + 100* AC.Wing1.Airfoil.t_c^4;
    Phi_w_2 = 2.7*AC.Wing2.Airfoil.t_c + 100* AC.Wing2.Airfoil.t_c^4;
    Cf_airfoil_1 = getCf ( AC.Wing1.Reynolds, AC.Wing1.Airfoil.transition_c);
    Cf_airfoil_2 = getCf ( AC.Wing2.Reynolds, AC.Wing2.Airfoil.transition_c);
    D.Cdp_min_1 = 2*Cf_airfoil_1*(1+Phi_w_1);
    D.Cdp_min_2 = 2*Cf_airfoil_2*(1+Phi_w_2);

    
%Cl for minimum Cd profile
    D.Cli_1     = AC.Wing1.Airfoil.Cldesign;
    D.Cli_2     = AC.Wing2.Airfoil.Cldesign;
 
    
    
%Extrapolated profile drag increment at stalling angle of attack
    W_WTO = prod([Parameters.fuelFraction(1:4).value]);
    HighCLSpeed = sqrt((2*DP.Wing1_Wing2*AC.Weight.MTOW*W_WTO*CST.GravitySI)/(rho*AC.Wing1.Sw*AC.Wing1.Airfoil.Cl_max*cosd(AC.Wing1.Sweep_14)*0.75));
    AirfoilReynolds1  = sqrt(HighCLSpeed^2*Parameters.q1_qinf) * AC.Wing1.CMA / nu;
    AirfoilReynolds2  = sqrt(HighCLSpeed^2*Parameters.q2_qinf) * AC.Wing2.CMA / nu;
    if AirfoilReynolds1 < 1e7
        deltaLCdp_ref_1 = (67*AC.Wing1.Airfoil.Cl_max)/(log10(AirfoilReynolds1))^4.5 - 0.0046*(1+2.75*AC.Wing1.Airfoil.t_c);
    else
        deltaLCdp_ref_1 = 0.01*AC.Wing1.Airfoil.Cl_max - 0.0046*(1+2.75*AC.Wing1.Airfoil.t_c);
    end
    if AirfoilReynolds2 < 1e7
        deltaLCdp_ref_2 = (67*AC.Wing2.Airfoil.Cl_max)/(log10(AirfoilReynolds2))^4.5 - 0.0046*(1+2.75*AC.Wing2.Airfoil.t_c);
    else
        deltaLCdp_ref_2 = 0.01*AC.Wing2.Airfoil.Cl_max - 0.0046*(1+2.75*AC.Wing2.Airfoil.t_c);
    end
    
    run generalized_Profile_Drag
    cl_1 = linspace(-0.25,AC.Wing1.Airfoil.Cl_max,50);
    cl_2 = linspace(-0.25,AC.Wing2.Airfoil.Cl_max,50);
    D.deltaLCdp_1 = deltaLCdp_ref_1 .* interp1(generalizedProfileDrag(:,1),generalizedProfileDrag(:,2),...
                    ((cl_1-D.Cli_1)./(AC.Wing1.Airfoil.Cl_max-D.Cli_1)).^2); %#ok<NODEF>
    D.deltaLCdp_2 = deltaLCdp_ref_2 .* interp1(generalizedProfileDrag(:,1),generalizedProfileDrag(:,2),...
                    ((cl_2-D.Cli_2)./(AC.Wing2.Airfoil.Cl_max-D.Cli_2)).^2);
    
if false          
    figure() %#ok<UNRCH>
    hold on
    plot(cl_1,D.Cdp_min_1+D.deltaLCdp_1)
    plot(cl_2,D.Cdp_min_2+D.deltaLCdp_2)
    title('Profile Drag of an SC(3)-0712 Airfoil Section','interpreter','latex')
    xlabel('$C_l$','interpreter','latex')
    ylabel('$C_d$','interpreter','latex')
    legend('Lead Wing','Rear Wing','Location','southeast')
    legend('boxoff')
    saveFigure(ME.FiguresFolder,'AirfoilProfileDrag')
end
    
    
%WINGS
    %Three dimensional profile drag coefficient derived from data at MAC airfoil
    D.Cdp_wing_min_1 = 2*Cf_airfoil_1*(1+Phi_w_1*cosd(AC.Wing1.Sweep_12));
    D.Cdp_wing_min_2 = 2*Cf_airfoil_2*(1+Phi_w_2*cosd(AC.Wing2.Sweep_12));
    D.Cli_wing_1 = AC.Wing1.Airfoil.Cldesign*cosd(AC.Wing1.Sweep_12);
    D.Cli_wing_2 = AC.Wing2.Airfoil.Cldesign*cosd(AC.Wing2.Sweep_12);
    D.CDp_wing_1 = D.Cdp_wing_min_1 * (AC.Wing1.Snet/AC.Wing1.Sw) + 0.75*deltaLCdp_ref_1*((AC.Wing1.CL_wf-D.Cli_wing_1)/(AC.Wing1.CLmax*cosd(AC.Wing1.Sweep_12)-D.Cli_wing_1))^2;
    D.CDp_wing_2 = D.Cdp_wing_min_2 * (AC.Wing2.Snet/AC.Wing2.Sw) + 0.75*deltaLCdp_ref_2*((AC.Wing2.CL_wf-D.Cli_wing_2)/(AC.Wing2.CLmax*cosd(AC.Wing2.Sweep_12)-D.Cli_wing_2))^2;

   
    
    
%FUSELAGE
    fusReynolds = DP.CruiseSpeed * AC.Fuselage.fusLength / nu;
    Cf_fus = getCf(fusReynolds, 0.05);
    Df_eff     = sqrt(4*AC.Fuselage.frontArea/pi);
    lambda_eff = min([AC.Fuselage.fusLength/Df_eff,(AC.Fuselage.la+AC.Fuselage.ln)/Df_eff+2]);
    Phi_f      = 2.2/lambda_eff^1.5 + 3.8/lambda_eff^3;


    DS.fuselageProfileBasic    = Cf_fus*AC.Fuselage.Swet*(1+Phi_f); 
    DS.fuselageProfileDeltaAoA = AC.Fuselage.A_I * abs(sind(AC.Fuselage.fuselage_AoA)^3)+ AC.Fuselage.A_II * ...
                                 abs(sind(AC.Fuselage.fuselage_AoA - AC.Fuselage.beta)^3)/cosd(AC.Fuselage.beta);

%     alpha_f_prime = AC.Fuselage.beta * (sqrt(AC.Fuselage.A_I/AC.Fuselage.A_II)-1)/(AC.Fuselage.A_I/AC.Fuselage.A_II-1);
%     [Wf_Wto_fit,gof] = fit(parametro(indexPar)',Wf_Wto(indexPar)','a*x^2+b*x+c','StartPoint',[0.1 0.1 0],'Lower',[0 -Inf,0],'Upper',[Inf Inf 0]); %order-->[a,b,c]

%NACELLES
    %4 Contribuciones:  *Fan cowling & Gas Generator cowling & Plug & Pylon
    nacReynolds  = sqrt(DP.CruiseSpeed^2*Parameters.q1_qinf) * AC.Engine.Length / nu;
    Cf_nac       = getCf(nacReynolds,0);
    %FAN COWLING
        SwetNacelles = pi * AC.Engine.Diameter * 0.25*AC.Engine.Length;
        DS.nacellesProfileFan = 1.25*AC.Engine.Number*Cf_nac*SwetNacelles;
    %GAS GENERATOR COWLING
        %Included in effective thrust loss
    %PLUG
        %Ususally considered as a loss in engine fross trhust
    %PYLON
%         DS.nacellesProfilePylon = AC.Engine.Number*Cf_nac*(1+2.75*AC.Engine.Pylon_t_c*cosd(AC.Engine.Pylon_Sweep)^2)*AC.Engine.Pylon_Swet;
        DS.nacellesProfilePylon = 0;
        %https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19870009105.pdf
 

%VERTICAL TAIL PLANE
    VTPReynolds   = DP.CruiseSpeed * AC.VTP.CMG / nu;
	Cf_VTP        = getCf( VTPReynolds, 0);
	DS.VTPProfile = 2*Cf_VTP*(1+2.75*AC.VTP.t_c*cosd(AC.VTP.Sweep_12)^2)*AC.VTP.Sw;
%     DS.VTPProfile = 0;
%     if ~ismember('AC:notDefinedVTP',ME.errorList)
%         ME.errorList{end+1} = 'AC:notDefinedVTP';
%         warning('AC:notDefinedVTP','Cuando se defina el VTP y se metan los datos hay que descomentar su resistencia.')
%     end


%PLOT EACH CONTRIBUTION
if dragPlotFlag
    X = [Parameters.q1_qinf * AC.Wing1.Sw * D.CDp_wing_1, ...
         Parameters.q2_qinf * AC.Wing2.Sw * D.CDp_wing_2, ...
         DS.nacellesProfileFan, ...
         DS.fuselageProfileBasic, DS.fuselageProfileDeltaAoA,...
         DS.VTPProfile];
    labels =  {'Ala delantera: ';...
               'Ala trasera: ';...
               'Fricción de las góndolas: ';...
               'Fricción del fuselaje: ';'Ángulo de ataque del fuselaje: ';...
               'Estabilizador vertical: '};
    fProfile = drawCustomPieChart(X,labels,false,false);
    title('Contribuciones a la resistencia parasita','interpreter','latex')
    saveFigure(ME.FiguresFolder,'ProfileDragPie') 
end


%%INTERFERENCE EFFECTS
%WETTED AREA CORRECTIONS
    %Nada, al considerar gross wetted areas para todo, ya es mayor de la real y así se incorpora una "correccion"
%WING/FUSELAGE INTERFERENCE
    fusDiameter = (AC.Fuselage.fusHeight+AC.Fuselage.fusWidth)/2; %Fuselage diameter in [m]
    %Inducida
        eta_fus_1 = fusDiameter/AC.Wing1.WingSpan;
        eta_fus_2 = fusDiameter/AC.Wing2.WingSpan;
        D.deltaICDv_1 = (0.55*eta_fus_1*(2-pi*eta_fus_1)*AC.Wing1.CL_wf^2)/((1+AC.Wing1.TaperRatio)*pi*AC.Wing1.AspectRatio);
        D.deltaICDv_2 = (0.55*eta_fus_2*(2-pi*eta_fus_2)*AC.Wing2.CL_wf^2)/((1+AC.Wing2.TaperRatio)*pi*AC.Wing2.AspectRatio);
    %Viscosa debida al engordamiento de la capa limite en la interseccion
        Cci_1 = 4.5*AC.Wing1.RootChord; % Length for both wing halves of the wing/fuselage intersection line, aprox 4.5 times the root chord
        Cci_2 = 4.5*AC.Wing2.RootChord; % Length for both wing halves of the wing/fuselage intersection line, aprox 4.5 times the root chord
        RootReynolds_1 = sqrt(DP.CruiseSpeed^2*Parameters.q1_qinf) * AC.Wing1.RootChord / nu;
        RootReynolds_2 = sqrt(DP.CruiseSpeed^2*Parameters.q2_qinf) * AC.Wing2.RootChord / nu;
        Cf_root_1 = getCf( RootReynolds_1, 0);
        Cf_root_2 = getCf( RootReynolds_2, 0);
        DS.deltaIDp_boundary_1 = 1.5*Cf_root_1*AC.Wing1.Airfoil.t_c*AC.Wing1.RootChord*Cci_1*cosd(AC.Wing1.Sweep_12)^2;
        DS.deltaIDp_boundary_2 = 1.5*Cf_root_2*AC.Wing2.Airfoil.t_c*AC.Wing2.RootChord*Cci_2*cosd(AC.Wing2.Sweep_12)^2;
    %Viscosa debida al incremento de velocidad sobre la cara superior del ala frente a la parte de abajo debido a la sustentacion
        CMAReynolds_1 = sqrt(DP.CruiseSpeed^2*Parameters.q1_qinf) * AC.Wing1.RootChord / nu;
        CMAReynolds_2 = sqrt(DP.CruiseSpeed^2*Parameters.q2_qinf) * AC.Wing2.RootChord / nu;
        Cf_CMA_1 = getCf( CMAReynolds_1, 0);
        Cf_CMA_2 = getCf( CMAReynolds_2, 0);
        DS.deltaIDp_speed_1 = -0.81*Cf_CMA_1*AC.Wing1.CL_wf*AC.Wing1.RootChord*fusDiameter; %Valid for high wing
        DS.deltaIDp_speed_2 = -0.81*Cf_CMA_2*AC.Wing2.CL_wf*AC.Wing2.RootChord*fusDiameter; %Valid for high wing
%NACELLES/AIRFRAME INTERFERENCE
        %Jet engines
%         DS.deltaIDp_nacelles = 0.2*(DS.nacellesProfileFan + DS.nacellesProfilePylon);
        DS.deltaIDp_nacelles = 0;
%TAILPLANE/AIRFRAME INTERFERENCE --> Wing2 / Wing 1 interference
        
%PLOT EACH CONTRIBUTION
if dragPlotFlag
    X = [Parameters.q1_qinf * AC.Wing1.Sw * D.deltaICDv_1, ...
         Parameters.q2_qinf * AC.Wing2.Sw * D.deltaICDv_2, ...
         DS.deltaIDp_boundary_1, DS.deltaIDp_boundary_2, ...
         abs(DS.deltaIDp_speed_1), abs(DS.deltaIDp_speed_2)];
    labels =  {'Resistencia inducida del ala delantera: ';...
               'Resistencia inducida del ala trasera: ';...
               'Capa limite del ala delantera: ';'Capa limite del ala trasera: ';...
               'Diferencia de velocidades ala delantera: ';'Diferencia de velocidades ala trasera: '};
    fInterferences = drawCustomPieChart(X,labels,false,false);
    title('Contribuciones a la resistencia de las interferencias','interpreter','latex')
    saveFigure(ME.FiguresFolder,'InterferencesDragPie') 
end

%%PROTUBERANCES AND IMPERFECTIONS
%FIXED UNDERCARRIAGE
%CANOPIES AND WINDSHIELDS
    %Canopies (Las cubiertas de las cabinas que sobresalen, por ejemplo en los cazas) --> No hay
    %Windshields
        DS.windshieldsProtuberance = 0.02*DS.fuselageProfileBasic; % 2-3% del fuselage drag
%WHEELS
%STREAMLINE STRUTS
%POWERPLANT INSTALLATION
%EXTERNAL FUEL TANKS
%OTHER EFFECTS, SURFACE IMPERFECTIONS
    %Meter los porcentajes que dice
    DS.Imperfections = DS.windshieldsProtuberance + ...
					   0.06 * Parameters.q1_qinf * AC.Wing1.Sw * D.CDp_wing_1 + ...
                       0.06 * Parameters.q2_qinf * AC.Wing2.Sw * D.CDp_wing_2 + ...
                       0.07 * DS.fuselageProfileBasic + ...
                       0.15 * DS.nacellesProfileFan + ...
                       0.01 * (Parameters.q1_qinf * AC.Wing1.Sw * D.CDp_wing_1 + ...
                               Parameters.q2_qinf * AC.Wing2.Sw * D.CDp_wing_2 + ...
                               DS.fuselageProfileBasic + DS.fuselageProfileDeltaAoA + ...
                               DS.nacellesProfileFan + DS.nacellesProfilePylon + ...
                               DS.VTPProfile);

                           
                           
%PLOT EACH CONTRIBUTION
if dragPlotFlag
    X = [DS.windshieldsProtuberance, ...
         0.06 * Parameters.q1_qinf * AC.Wing1.Sw * D.CDp_wing_1, ...
         0.06 * Parameters.q2_qinf * AC.Wing2.Sw * D.CDp_wing_2, ...
         0.15 * DS.nacellesProfileFan, ...
         0.07 * DS.fuselageProfileBasic, ...
         0.01 * (Parameters.q1_qinf * AC.Wing1.Sw * D.CDp_wing_1 + ...
                 Parameters.q2_qinf * AC.Wing2.Sw * D.CDp_wing_2 + ...
                 DS.fuselageProfileBasic + DS.fuselageProfileDeltaAoA + ...
                 DS.nacellesProfileFan + DS.nacellesProfilePylon + ...
                 DS.VTPProfile)];
    labels =  {'Ventanillas de los pilotos: ';...
               'Ala delantera: ';...
               'Ala trasera: ';...
               'Montaje de los motores: ';...
               'Fuselaje y VTP: ';...
               'Otros sistemas: '};
    fImperfections = drawCustomPieChart(X,labels,false,false);
    title('Contribuciones a la resistencia de las imperfecciones','interpreter','latex')
    saveFigure(ME.FiguresFolder,'ImperfectionsDragPie') 
end




%%TOTAL DRAG
    DS.Induced = Parameters.q1_qinf*AC.Wing1.Sw*(D.CDv_1 + D.deltaEpsilonCDv_1) + ...
                 Parameters.q2_qinf*AC.Wing2.Sw*(D.CDv_2 + D.deltaEpsilonCDv_2) + ...
                 DS.fuselageInduced;

    DS.Profile = Parameters.q1_qinf * AC.Wing1.Sw * D.CDp_wing_1 + ...
                 Parameters.q2_qinf * AC.Wing2.Sw * D.CDp_wing_2 + ...
                 DS.fuselageProfileBasic + DS.fuselageProfileDeltaAoA + ...
                 DS.nacellesProfileFan + DS.nacellesProfilePylon + ...
                 DS.VTPProfile;

    DS.Interferences = Parameters.q1_qinf * AC.Wing1.Sw * D.deltaICDv_1 + ...
                       Parameters.q2_qinf * AC.Wing2.Sw * D.deltaICDv_2 + ...
                       DS.deltaIDp_boundary_1 + DS.deltaIDp_boundary_2 + ...
                       DS.deltaIDp_speed_1 + DS.deltaIDp_speed_2 + ...
                       DS.deltaIDp_nacelles;
                   
    CD = (DS.Induced + DS.Profile + DS.Interferences + DS.Imperfections)/AC.Wing.Sw;

    
if dragPlotFlag
    X = [DS.Induced, DS.Interferences, DS.Profile, DS.Imperfections];
    labels =  {'Resistencia Inducida: ';...
               'Resistencia debida a las interferencias: ';...
               'Resistencia Parásita: ';...
               'Resistencia debida a las imperfecciones: '};
    fDrag = drawCustomPieChart(X,labels,false,false);
    title('Contribuciones a la resistencia','interpreter','latex')
    saveFigure(ME.FiguresFolder,'DragPie') 
end
    
% clear D DS

    
%%CLEAR WORKSPACE
clear a T P rho nu cl_1 cl_2 Cci_1 Cci_2 Df_eff eta_cp_1 eta_cp_2 eta_fus_1 eta_fus_2
clear Cf_airfoil_1 Cf_airfoil_2 Cf_CMA_1 Cf_CMA_2 Cf_fus Cf_nac Cf_root_1 Cf_root_2 Cf_VTP
clear fusDiameter fusReynolds CMAReynolds_1 CMAReynolds_2 generalizedProfileDrag lambda_eff
clear W_WTO HighCLSpeed delta_1 delta_2 deltaLCdp_ref_1 deltaLCdp_ref_2 nacReynolds SwetNacelles
clear Phi_f Phi_w_1 Phi_w_2 AirfoilReynolds1 AirfoilReynolds2 RootReynolds_1 RootReynolds_2 v w 


end

function [f] = drawCustomPieChart(X,labels,explodeFlag,partialFlag)
    
    %Create figure
    f=figure();
    
    %Check if explode graph
    if explodeFlag
        explode = ones(1,length(X));
    else
        explode = zeros(1,length(X));
    end
    
    %Check if partial graph
    if ~partialFlag
        while sum(X)<1
            X=10.*X;
        end
    end
    
    %Create pie chart
    h=pie(X,explode);
    
    %Store Precalculated Percent Values
    hText = findobj(h,'Type','text'); % text object handles
    percentValues = get(hText,'String'); % percent values
    
    %Combine Percent Values and Additional Text
        %labels = {'Item A: ';'Item B: ';'Item C: '}; % strings
    combinedtxt = strcat(labels,percentValues); % strings and percent values
    
    %Store the text Extent property values for the current labels
    oldExtents_cell = get(hText,'Extent');  % cell array
    oldExtents = cell2mat(oldExtents_cell); % numeric array
    for i=1:length(X)
        hText(i).String = combinedtxt(i);
    end
    
    %Determine Horizontal Distance to Move Each Label
    newExtents_cell = get(hText,'Extent'); % cell array
    newExtents = cell2mat(newExtents_cell); % numeric array 
    width_change = newExtents(:,3)-oldExtents(:,3);
    signValues = sign(oldExtents(:,1));
    offset = signValues.*(width_change/2);
    
    %Position New Label
    textPositions_cell = get(hText,{'Position'}); % cell array
    textPositions = cell2mat(textPositions_cell); % numeric array
    textPositions(:,1) = textPositions(:,1) + offset; % add offset 
    for i=1:length(X)
        hText(i).Position = textPositions(i,:);
    end

end
