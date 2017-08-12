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
        DP.Payload     =  800; %[kg] <-- 6 pax (80+50 each)
        DP.Range       = 10e3; %[km]
        DP.MaxSpeed    =  257; %[m/s]  (Mach=0.87)
        DP.TOFL        = 1200; %[m]    Take-Off Field Length, from similar planes: max-->1972m, min-->956m, mean-->1488m
        DP.LFL         =  800; %[m]    Landing Field Length, from similar planes: max-->1015m, min-->631m, mean-->733m
        
        
    %% PLOTTING OPTIONS
        DP.ShowReportFigures      = false; %Show all the available figures for reports [true] or only the most relevant ones [false]
        DP.ShowAircraftLayout     = false;  %Show a layout of the aircraft and the wings
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
        DP.CruiseSpeed         =    250; %[m/s] (Mach=0.85)
        DP.CruiseEfficiency    =     14; %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.CruiseTSFC          =  0.597; %0.625 %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Loiter
        DP.LoiterTime          =  30*60;  %[s]
        DP.LoiterEfficiency    =     14;  %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.LoiterTSFC          =    0.35;  %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Low Height
        DP.LowHeightEfficiency =     12;  %Lower than CruiseEfficiency because of lower altitude (Typically lower than 10.000ft)
        DP.LowHeightTSFC       =    0.30;  %[lbm/(lbf·h)] Not sure if higher or lower than CruiseTSFC as the fly will probably be carried out below
                                           %              10.000ft and at a max of 250kts (128.611 m/s) in accordance with FAA regulations
        
        %Fuel Reserves & Alternate Airport
        DP.AlternateRange      = 370; %[km] Range to alternate airport (200 nautic miles --> 370km)
        
        %Landing
        DP.StallSpeed_L    = NaN; %[m/s] -- Landing Stall Speed
        
        %Wing - Airfoil
        DP.Wing1_Wing2     = 0.75; %[-] - Parte de la sustentación que se lleva el ala delantera,
        DP.Incidence_1     = 2.85; %[º] - Degrees of angle of the wing/body incidence of wing 1 at root section
        DP.Incidence_2     = 1.10; %[º] - Degrees of angle of the wing/body incidence of wing 2 at root section
        DP.AspectRatio     =    8; %[-] - Aspect ratio, from similar planes: max-->9.7166, min-->8.0139, mean-->9.0017
        DP.TaperRatio      =  0.6; %[-]
        DP.Dihedral        =  0.0; %[º]
        DP.TipTwist        = -5.0; %[º] - Positive twist: nose rotated upwards (Wash-in). Negative twist: nose rotated downwards (Wash-out)
        DP.Stagger         = 6.75; %[m]
        DP.VerticalGap     =    0; %[m]
        DP.Wing1LongPos    =    4; %[m] Longitudinal position of the first wing
        DP.Sweep_14        =   30; %[º] Flecha en la linea 1/4 
        DP.CLmax           =  1.2; %Porque si, hay que calcularlo bien... los valores estimados en crucero son muy bajos por ser la velocidad muy alta
        DP.CLmax_TO        =  2.0; %From similar planes: max-->2.3447, min-->1.5414, mean-->2.0622
        DP.CLmax_L         =  2.8; %From similar planes: max-->3.7689, min-->2.2523, mean-->3.0764
        
        %Weight
        DP.x_cg            =    11; %[m] Posición longitudinal del centro de gravedad, se debe calcular, solo es para que no pete el código.
        DP.y_cg            =     0; %[m] Posición lateral del centro de gravedad.
        DP.z_cg            =     0; %[m] Posición vertical del centro de gravedad, se debe calcular, solo es para que no pete el código.
        DP.MLW_MTOW        =  0.85; %From SP: min-->0.7900, max-->0.9267, mean-->0.8753
        DP.MRW_MTOW        = 1.005; %From SP: min-->0.9918, max-->1.0286, mean-->1.0051
        DP.EWnew_EWold     = 0.95;  %Weight reduction of the empty weight as being fully manufatured in composite materials
                                    %From: http://www.compositesworld.com/news/revolutionary-fuselage-concept-unveiled-by-mtorres --> Fuselage weight
                                    %reduction estimated between 10% and 30%, and from Roskam Part V, Chapter 2 page 11, average FuselageWeight/EW=0.2
                                    %so, total EW reduction estimated between 2% and 6%, we chose a confident average of 5% --> EWnew=0.95*EWold
        
        %Engines
        DP.EngineNumber    = 2;
        DP.EngineModel     = 'Snecma Silvercrest 2D';%'Rolls-Royce AE 3007A1E';
        
        %Crew
        DP.CrewNumber     = 2;
        
        % Fuselage Shape
        DP.fuselage_AoA   =     0; %[º] Angle of attack of the fuselage
        DP.fusWidth       =  2.50; %[m]
        DP.fusHeight      = 2.735; %[m]
        DP.fusLength      =    23; %12.15 + 3*(DP.fusWidth+DP.fusHeight)/2;
        DP.totalFusVolume =  pi*DP.fusWidth*DP.fusHeight*DP.fusLength/4; %[m^3]

        
        
        
        
        
    case 11  %11. Flying boats, amphibious, float airplanes
       
        
        %% PLOTTING OPTIONS
        DP.ShowReportFigures      = true; %Show all the available figures for reports [true] or only the most relevant ones [false]
        DP.selectDesignPoint      = false; %Ask user to select design point [true] or use the saved value [false]
        DP.showRoskamRequirements = false; %Show the Take-Off and Landing requirements obtained with Roskam constants [true] or only the SP ones [false]
        
        
        
        
end