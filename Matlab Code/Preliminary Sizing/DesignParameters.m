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
        DP.Payload     =  800; %[kg] <-- 6 pax
        DP.Range       = 10e6; %[m]
        DP.MaxSpeed    =  257; %[m/s]  (Mach=0.87)
        DP.TOFL        = 1200; %[m]    Take-Off Field Length, from similar planes: max-->1972m, min-->956m, mean-->1488m
        DP.LFL         =  800; %[m]    Landing Field Length, from similar planes: max-->1015m, min-->631m, mean-->733m
        
        
    %% PLOTTING OPTIONS
        DP.ShowReportFigures      = true; %Show all the available figures for reports [true] or only the most relevant ones [false]
        DP.selectDesignPoint      = false; %Ask user to select design point [true] or use the saved value [false]
        DP.showRoskamRequirements = false; %Show the Take-Off and Landing requirements obtained with Roskam constants [true] or only the SP ones [false]
        warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
        
        
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
        DP.CruiseTSFC          =   0.60; %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Loiter
        DP.LoiterTime          =  30*60;  %[s]
        DP.LoiterEfficiency    =     14;  %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.LoiterTSFC          =    0.35;  %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Low Height
        DP.LowHeightEfficiency =     12;  %Lower than CruiseEfficiency because of lower altitude (Typically lower than 10.000ft)
        DP.LowHeightTSFC       =    0.30;  %[lbm/(lbf·h)] Not sure if higher or lower than CruiseTSFC as the fly will probably be carried out below
                                          %              10.000ft and at a max of 250kts (128.611 m/s) in accordance with FAA regulations
        
        %Fuel Reserves & Alternate Airport
        DP.AlternateRange      = 370*1e3; %[m] Range to alternate airport (200 nautic miles --> 370km)
        
        %Landing
        DP.StallSpeed_L    = NaN; %[m/s] -- Landing Stall Speed
        
        %Wing - Airfoil
        DP.AspectRatio     =   9; %[-] - Aspect ratio, from similar planes: max-->9.7166, min-->8.0139, mean-->9.0017
        DP.CLmax           = 1.4; %Porque si, hay que calcularlo bien... los valores estimados en crucero son muy bajos por ser la velocidad muy alta
        DP.CLmax_TO        = 2.0; %From similar planes: max-->2.3447, min-->1.5414, mean-->2.0622
        DP.CLmax_L         = 2.8; %From similar planes: max-->3.7689, min-->2.2523, mean-->3.0764
        
        %Weight
        DP.MLW_MTOW        = 0.85;  %From SP: min-->0.7900, max-->0.9267, mean-->0.8753
        DP.MRW_MTOW        = 1.005; %From SP: min-->0.9918, max-->1.0286, mean-->1.0051
        DP.EWnew_EWold     = 0.95;  %Weight reduction of the empty weight as being fully manufatured in composite materials
                                    %From: http://www.compositesworld.com/news/revolutionary-fuselage-concept-unveiled-by-mtorres --> Fuselage weight
                                    %reduction estimated between 10% and 30%, and from Roskam Part V, Chapter 2 page 11, average FuselageWeight/EW=0.2
                                    %so, total EW reduction estimated between 2% and 6%, we chose a confident average of 5% --> EWnew=0.95*EWold
        
        %Engines
        DP.NumberEngines    = 2;
        % [ ] 2*Rolls-Royce AE 3007A1E --> Embraer Legacy family
        DP.EngineOptions(1).Name   = 'Rolls-Royce AE 3007A1E';
        DP.EngineOptions(1).Thrust =    40; %[kN] at sea level
        DP.EngineOptions(1).Weight = 751.6; %[kg] dry weight
        DP.EngineOptions(1).TSFC   = 0.625; %[lb/(lbf·h)] at cruise 0.33 lb/(lbf·h) at TO
        % [ ] 2*Snecma Silvercrest --> Cessna Citation Hemisphere
        DP.EngineOptions(2).Name   = 'Snecma Silvercrest 2C';
        DP.EngineOptions(2).Thrust =  53.4; %[kN] at sea level
        DP.EngineOptions(2).Weight =  1040; %[kg] dry weight
        DP.EngineOptions(2).TSFC   = 0.628; %[lb/(lbf·h)] at cruise
        % [X] 2*Snecma Silvercrest --> Dassault Falcon 5X
        DP.EngineOptions(3).Name   = 'Snecma Silvercrest 2D';
        DP.EngineOptions(3).Thrust =  50.9; %[kN] at sea level
        DP.EngineOptions(3).Weight =  1040; %[kg] dry weight
        DP.EngineOptions(3).TSFC   = 0.628; %[lb/(lbf·h)] at cruise
        % [ ] 2*Rolls-Royce BR710A2-20 --> Bombardier Global 5000
        DP.EngineOptions(4).Name   = 'Rolls-Royce BR710A2-20';
        DP.EngineOptions(4).Thrust = 65.6; %[kN] at sea level
        DP.EngineOptions(4).Weight = 1633; %[kg] dry weight
        DP.EngineOptions(4).TSFC   = 0.63; %[lb/(lbf·h)] at cruise
        % [ ] 2*P&W Canada PW800 --> Gulfstream G500/G600 (67.36/69.75 kN)
        DP.EngineOptions(5).Name   = 'P&W Canada PW814GA';
        DP.EngineOptions(5).Thrust = 67.36; %[kN] at sea level
        DP.EngineOptions(5).Weight =  1422; %[kg] dry weight
        DP.EngineOptions(5).TSFC   =    []; %[lb/(lbf·h)] at cruise
        % [ ] 2*P&W PW1215G --> Mitsubishi Regional Jet
        DP.EngineOptions(6).Name   = 'P&W PW1215G';
        DP.EngineOptions(6).Thrust = 66.7; %[kN] at sea level
        DP.EngineOptions(6).Weight =   []; %[kg] dry weight
        DP.EngineOptions(6).TSFC   =   []; %[lb/(lbf·h)] at cruise
        % [ ] 2*GE Passport 20–17BB1A (44 to 89 kN)       -->  ~90kN ~kg ~lb/(lbf·h)--> Bombardier Global 7000/8000 (73.4 kN)
        DP.EngineOptions(7).Name   = 'GE Passport 20–17BB1A';
        DP.EngineOptions(7).Thrust = 74.8; %[kN] at sea level
        DP.EngineOptions(7).Weight = []; %[kg] dry weight
        DP.EngineOptions(7).TSFC   = []; %[lb/(lbf·h)] at cruise
        
        %Crew
        DP.CrewNumber     = 2;
        
        
        
        
    case 11  %11. Flying boats, amphibious, float airplanes
       
        
        %% PLOTTING OPTIONS
        DP.ShowReportFigures      = true; %Show all the available figures for reports [true] or only the most relevant ones [false]
        DP.selectDesignPoint      = false; %Ask user to select design point [true] or use the saved value [false]
        DP.showRoskamRequirements = false; %Show the Take-Off and Landing requirements obtained with Roskam constants [true] or only the SP ones [false]
        
        
        
        
end