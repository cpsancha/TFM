%Apendice F

%% VORTEX INDUCED DRAG
%METHOD A --> A. Garner [Ref. F-39]
%Spanwise center of pressure    
    eta_cp_1 = trapz(AC.Wing1.eta,AC.Wing1.eta.*AC.Wing1.cla.*AC.Wing1.c./AC.Wing1.CL_wf./AC.Wing1.CMG);
    eta_cp_2 = trapz(AC.Wing2.eta,AC.Wing2.eta.*AC.Wing1.cla.*AC.Wing2.c./AC.Wing2.CL_wf./AC.Wing2.CMG);
    
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
    D.deltaEpsilonCDv_1 = AC.Wing1.CL_wf*(deg2rad(AC.Wing1.TipTwist)*AC.Wing1.CL_alpha_wf/E1)*v+(deg2rad(AC.Wing1.TipTwist)*AC.Wing1.CL_alpha_wf/E1)^2*w;
    D.deltaEpsilonCDv_2 = AC.Wing2.CL_wf*(deg2rad(AC.Wing2.TipTwist)*AC.Wing2.CL_alpha_wf/E2)*v+(deg2rad(AC.Wing2.TipTwist)*AC.Wing2.CL_alpha_wf/E2)^2*w;
    
%INDUCED DRAG DUE TO WING TIP CORRECTION
    %No se corrije porq no se sabe como
    
%INDUCED DRAG DUE TO FUSELAGE LIFT CORRECTION
    DS.fuselageInduced = 0.15*AC.Fuselage.fuselage_AoA^2*AC.Fuselage.Volume^(2/3);

%NACELLES CONTRIBUTION
    %No sense to estimate, interferences will be much larger
 
    
    
%% PROFILE DRAF OF SMOOTH, ISOLATED MAJOR COMPONENTS --> FLAT PLATE ANALOGY
%WING SECTIONS

%WINGS
%FUSELAGE
%NACELLES



%% INTERFERENCE EFFECTS
%WETTED AREA CORRECTIONS
%WING/FUSELAGE INTERFERENCE
%NACELLES/AIRFRAME INTERFERENCE
%TAILPLANE/AIRFRAME INTERFERENCE



%% PROTUBERANCES AND IMPERFECTIONS
%FIXED UNDERCARRIAGE
    %No tenemos en nuestro caso
%CANOPIES AND WINDSHIELDS (LAS VENTANILLAS DE LOS PILOTOS)
%POWERPLANT INSTALLATION
%EXTERNAL FUEL TANKS
    %No se aplica
%OTHER EFFECTS
    %Meter los porcentajes que dice















    
% clear C1_1 C1_2 C2_1 C2_2 C3_1 C3_2 f1 f2 c1 c2 y1 y2 eta1 eta2 Twist1 Twist2
% clear eta_cp_1 eta_cp_2 delta_1 delta_2