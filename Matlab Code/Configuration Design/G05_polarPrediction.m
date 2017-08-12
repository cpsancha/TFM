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
    [rho,a,T,P,nu,~] = atmos(DP.CruiseAltitude);
    
%Minimum Cd profile
    [D.Cdp_min_1,index_1] = min(AC.Wing1.Airfoil.Polar(2,:));
    [D.Cdp_min_2,index_2] = min(AC.Wing2.Airfoil.Polar(2,:));
    
%Cl for minimum Cd profile
    D.Cli_1 = AC.Wing1.Airfoil.Polar(1,index_1);
    D.Cli_2 = AC.Wing2.Airfoil.Polar(1,index_2);
    
%Extrapolated profile drag increment at stalling angle of attack
    HighCLSpeed = sqrt((2*DP.Wing1_Wing2*AC.Weight.Weight*CST.GravitySI)/(rho*AC.Wing1.Sw*AC.Wing1.Airfoil.Cl_max*0.75))/cosd(AC.Wing1.Sweep_14);
    Reynolds1  = HighCLSpeed * AC.Wing1.CMA / nu;
    Reynolds2  = HighCLSpeed * AC.Wing2.CMA / nu;
    if Reynolds1 < 1e7
        deltaLCdp_ref_1 = (67*AC.Wing1.Airfoil.Cl_max)/(log10(Reynolds1))^4.5 - 0.0046*(1+2.75*AC.Wing1.Airfoil.t_c);
    else
        deltaLCdp_ref_1 = 0.01*AC.Wing1.Airfoil.Cl_max - 0.0046*(1+2.75*AC.Wing1.Airfoil.t_c);
    end
    if Reynolds2 < 1e7
        deltaLCdp_ref_2 = (67*AC.Wing2.Airfoil.Cl_max)/(log10(Reynolds2))^4.5 - 0.0046*(1+2.75*AC.Wing2.Airfoil.t_c);
    else
        deltaLCdp_ref_2 = 0.01*AC.Wing2.Airfoil.Cl_max - 0.0046*(1+2.75*AC.Wing2.Airfoil.t_c);
    end
    
    run generalized_Profile_Drag
    cl_1 = linspace(-0.5,AC.Wing1.Airfoil.Cl_max,50);
    cl_2 = linspace(-0.5,AC.Wing2.Airfoil.Cl_max,50);
    D.deltaLCdp_1 = deltaLCdp_ref_1 .* interp1(generalizedProfileDrag(:,1),generalizedProfileDrag(:,2),...
                    ((cl_1-D.Cli_1)./(AC.Wing1.Airfoil.Cl_max-D.Cli_1)).^2);
    D.deltaLCdp_2 = deltaLCdp_ref_2 .* interp1(generalizedProfileDrag(:,1),generalizedProfileDrag(:,2),...
                    ((cl_2-D.Cli_2)./(AC.Wing2.Airfoil.Cl_max-D.Cli_2)).^2);
    
%     figure()
%     hold on
%     plot(cl_1,D.Cdp_min_1+D.deltaLCdp_1,'r')
%     plot(cl_2,D.Cdp_min_2+D.deltaLCdp_2,'b')
%     title('Profile Drag of an Airfoil Section','interpreter','latex')
%     xlabel('$C_l$','interpreter','latex')
%     ylabel('$C_d$','interpreter','latex')
%     saveFigure(ME.FiguresFolder,'AirfoilProfileDrag')

    
    
%WINGS
    %Three dimensional profile drag coefficient derived from data at MAC airfoil
    D.CDp_1 = D.Cdp_min_1 * (AC.Wing1.Snet/AC.Wing1.Sw) + 0.75*deltaLCdp_ref_1*((AC.Wing1.CL_wf-D.Cli_1)/(AC.Wing1.CLmax-D.Cli_1))^2;
    D.CDp_2 = D.Cdp_min_2 * (AC.Wing2.Snet/AC.Wing2.Sw) + 0.75*deltaLCdp_ref_2*((AC.Wing2.CL_wf-D.Cli_2)/(AC.Wing2.CLmax-D.Cli_2))^2;

    D.CD_1 = D.CDv_1 + D.CDp_1;
    D.CD_2 = D.CDv_2 + D.CDp_2;
    
%FUSELAGE
beta = 17; %grados
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