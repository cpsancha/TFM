% *************************************************************************
% *                                                                       *
% *                       DESIGN PARAMETERS                               *
% *                                                                       *
% *************************************************************************
% 
% This script is intended to be a compendium of all the design parameters
% to be manually decided, for simplicity.

switch ME.MissionType
    case 5  %5.Business Jets
        
    %% REQUIREMENTS
        DP.Payload     =  750; %[kg] <-- 6 pax (80+45 each)
        DP.Range       = 10e3; %[km]
        DP.MaxSpeed    =  250; %[m/s]  (Mach=0.85)
        DP.TOFL        = 1200; %[m]    Take-Off Field Length, from similar planes: max-->1972m, min-->956m, mean-->1488m
        DP.LFL         =  800; %[m]    Landing Field Length, from similar planes: max-->1015m, min-->631m, mean-->733m
        
        
    %% PLOTTING OPTIONS
        DP.ShowReportFigures      = false; %Show all the available figures for reports [true] or only the most relevant ones [false]
        DP.ShowAircraftLayout     = true;  %Show a layout of the aircraft and the wings
        DP.selectDesignPoint      = false; %Ask user to select design point [true] or use the saved value [false]
        DP.showRoskamRequirements = false; %Show the Take-Off and Landing requirements obtained with Roskam constants [true] or only the SP ones [false]
        
        
    %% CUSTOM SELECTED PARAMETERS
        %Take-Off
        DP.StallSpeed_TO       =    NaN; %[m/s] - Take-Off Stall Speed
           
        %Climb
        DP.ClimbRate           =   12.7; %[m/s] (2500 ft/min)
        DP.ClimbSpeed          =    140; %Climb speed in m/s (horizontal) (270 kts)
        DP.ClimbTime           =    NaN; %[s]   Time to reach cruise altitude, now is calculated from the ClimbRate and CruiseAltitude
        
        %Cruise
        DP.StallSpeed          =    NaN; %[m/s] -- Cruise Stall Speed
        DP.CruiseAltitude      =   12e3; %[m]
        DP.CruiseSpeed         =    243; %[m/s] (875 km/h  Mach=0.825)
        DP.CruiseEfficiency    =     14; %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.CruiseTSFC          =  0.597; %0.625 %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Loiter
        DP.LoiterTime          =  30*60;  %[s]
        DP.LoiterEfficiency    =     14;  %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.LoiterTSFC          =   0.35;  %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Low Height
        DP.LowHeightEfficiency =     12;  %Lower than CruiseEfficiency because of lower altitude (Typically lower than 10.000ft)
        DP.LowHeightTSFC       =    0.30;  %[lbm/(lbf·h)] Not sure if higher or lower than CruiseTSFC as the fly will probably be carried out below
                                           %              10.000ft and at a max of 250kts (128.611 m/s) in accordance with FAA regulations
        
        %Fuel Reserves & Alternate Airport
        DP.AlternateRange      = 370; %[km] Range to alternate airport (200 nautic miles --> 370km)
        
        %Landing
        DP.StallSpeed_L    = NaN; %[m/s] -- Landing Stall Speed
        
        %Wing - Airfoil
        DP.Wing1_Wing2     = 0.85; %[-] - Parte de la sustentación que se lleva el ala delantera,
        DP.Incidence_1     =    2;%2.85; %[º] - Degrees of angle of the wing/body incidence of wing 1 at root section
        DP.Incidence_2     =   -2;%1.10; %[º] - Degrees of angle of the wing/body incidence of wing 2 at root section
        DP.WingLoading     =  400; %[kg/m^2] - Default value if selection not active
        DP.AspectRatio     =    8; %[-] - Aspect ratio, from similar planes: max-->9.7166, min-->8.0139, mean-->9.0017
        DP.TaperRatio      =  0.6; %[-]
        DP.Dihedral        =  0.0; %[º]
        DP.TipTwist        = -5.0; %[º] - Positive twist: nose rotated upwards (Wash-in). Negative twist: nose rotated downwards (Wash-out)
        DP.Stagger         =  7.0; %[m]
        DP.VerticalGap     =    0; %[m]
        DP.Wing1LongPos    = 3.65; %[m] Longitudinal position of the first wing
        DP.Sweep_14        = 27.5; %[º] Flecha en la linea 1/4
        DP.lowSpeedFlag    =false; %Disable function @getAirfoilData and set Clalpha=2*pi 
        DP.CLmax           =  1.2; %Porque si, hay que calcularlo bien... los valores estimados en crucero son muy bajos por ser la velocidad muy alta
        DP.CLmax_TO        =  2.0; %From similar planes: max-->2.3447, min-->1.5414, mean-->2.0622
        DP.CLmax_L         =  2.8; %From similar planes: max-->3.7689, min-->2.2523, mean-->3.0764
        
        %Flaps
        DP.cf_c            = 0.35;
        DP.cle_c           = 0.10;
        DP.delta_f         =   40; %[deg]
        
        %Weight
        DP.chooseWeights   =  true;
        DP.EW              =  8950; %[kg]
        DP.MTOW            = 20655; %[kg]
        DP.MFW             = 10650; %[kg]
        DP.x_cg            = 10.30; %[m] Posición longitudinal del centro de gravedad, se debe calcular, solo es para que no pete el código.
        DP.y_cg            =     0; %[m] Posición lateral del centro de gravedad.
        DP.z_cg            =     1; %[m] Posición vertical del centro de gravedad, se debe calcular, solo es para que no pete el código.
        DP.MLW_MTOW        =  0.85; %From SP: min-->0.7900, max-->0.9267, mean-->0.8753
        DP.MRW_MTOW        = 1.005; %From SP: min-->0.9918, max-->1.0286, mean-->1.0051
        DP.EWnew_EWold     = 0.95;  %Weight reduction of the empty weight as being fully manufatured in composite materials
                                    %From: http://www.compositesworld.com/news/revolutionary-fuselage-concept-unveiled-by-mtorres --> Fuselage weight
                                    %reduction estimated between 10% and 30%, and from Roskam Part V, Chapter 2 page 11, average FuselageWeight/EW=0.2
                                    %so, total EW reduction estimated between 2% and 6%, we chose a confident average of 5% --> EWnew=0.95*EWold
        
        %Engines
        DP.EngineNumber    = 2;
        DP.EngineModel     = 'Snecma Silvercrest 2D';%'Rolls-Royce AE 3007A1E';
        DP.xEngine         =   14;
        DP.yEngine         = 1.40;
        DP.zEngine         = 2.25;
        DP.Pylon_t_c       = 0.12; %NACA 0012
        DP.Pylon_Swet      = 0.25;
        DP.Pylon_Sweep     =   30; %[º]
        
        %Crew
        DP.CrewNumber     = 2;
        
        % Fuselage Shape
        DP.fuselage_AoA   =     0; %[º] Angle of attack of the fuselage
        DP.fusWidth       =  2.50; %[m]
        DP.fusHeight      = 2.735; %[m]
        DP.fusLength      =  18.6; %[m] Total length
        DP.cabLength      = 8.575; %[m] Length of the cabin including the gallery and lavatory,  Gallery: 1.48m   Cabin: 5.05m
        DP.cabHeight      =  2.00; %[m]
        DP.cabWidth       =  2.25; %[m]
        DP.galleryLength  =  1.48; %[m]
        DP.lavatoryLength = 2.045; %[m]
        DP.ln             = 3.625; %[m] Length of the nose
        DP.la             = DP.fusLength-DP.ln-DP.cabLength; %[m]
        DP.frontArea      = pi*DP.fusWidth*DP.fusHeight/4;   %[m^2]
        DP.totalFusVolume =     DP.frontArea*DP.fusLength;   %[m^3]
        DP.Swet           = 120.57; %[m^2] Exterior wet surface
        DP.A_I            =  28.13; %[m^2]
        DP.A_II           =  11.09; %[m^2]
        DP.tailConeAngle  =     17; %[º]
        DP.depositsArea   =  0.754; %[m^2]
        DP.depositsUsage  =   0.80; %Percentage of usage of the deposits for fuel
         
  
        % VTP
        DP.VTP_X_ac        = 17.0;
        DP.VTP_AspectRatio = 1.40;
        DP.VTP_Sweep_LE    = 45.0; %[º] Leading edge VTP sweep
        DP.VTP_Sweep_r     = 25.0; %[º] Rudder sweep
        DP.VTP_Sr_Sv       = 0.30;
        DP.VTP_TaperRatio  = 0.40;
        DP.VTP_deltar_max  = 35.0;
        DP.VTP_t_c         = 0.15;
        DP.VTP_Airfoil     ='0015';
        DP.VTP_b_14        = 2.50;
        DP.VTP_b_34        = 2.43;
        DP.VTP_h_14        = 2.30;
        DP.VTP_h_34        = 2.28;
        DP.VTP_Svertical   = 36.476;
        
        
        % Undercarriage
        DP.UnderCarriageNose =  4; %[m]
        DP.UnderCarriageMain = 12; %[m]
        
        
        
    case 11  %11. Flying boats, amphibious, float airplanes
       
        
        %% PLOTTING OPTIONS
        DP.ShowReportFigures      = true; %Show all the available figures for reports [true] or only the most relevant ones [false]
        DP.selectDesignPoint      = false; %Ask user to select design point [true] or use the saved value [false]
        DP.showRoskamRequirements = false; %Show the Take-Off and Landing requirements obtained with Roskam constants [true] or only the SP ones [false]

        
        
        
end