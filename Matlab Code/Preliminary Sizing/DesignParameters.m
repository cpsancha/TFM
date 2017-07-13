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
        DP.Payload       = 1250; %[kg]
        DP.Range         = 10e6; %[m]
       
        
    %% CUSTOM SELECTED PARAMETERS
        %Take-Off
        DP.TOFL          = 1200; %[m]   - Take-Off Field Length, from similar planes: max-->1972m, min-->956m, mean-->1488m
        DP.StallSpeed_TO =  NaN; %[m/s] - Take-Off Stall Speed
           
        %Climb
        DP.ClimbRate           =   12.7; %[m/s] (2500 ft/min)
        DP.ClimbSpeed          =    140; %Climb speed in m/s (horizontal) (270 kts)
        DP.ClimbTime           =    NaN; %[s]   Time to reach cruise altitude, now is calculated from the ClimbRate and CruiseAltitude
        
        %Cruise
        DP.StallSpeed          =    NaN; %[m/s] -- Cruise Stall Speed
        DP.CruiseAltitude      = 12.5e3; %[m]
        DP.CruiseSpeed         =    225; %[m/s]
        DP.CruiseEfficiency    =     14; %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.CruiseTSFC          =   0.50; %[lbm/(lbf�h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Loiter
        DP.LoiterTime          = 30*60;  %[s]
        DP.LoiterEfficiency    = 16;     %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.LoiterTSFC          = 0.4;    %[lbm/(lbf�h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Low Height
        DP.LowHeightEfficiency = 12;     %Lower than CruiseEfficiency because of lower altitude (Typically lower than 10.000ft)
        DP.LowHeightTSFC       = 0.3;    %[lbm/(lbf�h)] Not sure if higher or lower than CruiseTSFC as the fly will probably be carried out below
                                         %              10.000ft and at a max of 250kts (128.611 m/s) in accordance with FAA regulations
        
        %Fuel Reserves & Alternate Airport
        DP.AlternateRange  = 370*1e3; %[m] Range to alternate airport (200 nautic miles --> 370km)
        
        %Landing
        DP.LFL             = 700; %[m]   -- Landing Field Length, from similar planes: max-->1015m, min-->631m, mean-->733m
        DP.StallSpeed_L    = NaN; %[m/s] -- Landing Stall Speed
        
        %Wing - Airfoil
        DP.AspectRatio     =   9; %[-] - Aspect ratio, from similar planes: max-->9.7166, min-->8.0139, mean-->9.0017
        DP.CLmax           = 1.4; %Porque si, hay que calcularlo bien... los valores estimados en crucero son muy bajos por ser la velocidad muy alta
        DP.CLmax_TO        = 2.0; %From similar planes: max-->2.3447, min-->1.5414, mean-->2.0622
        DP.CLmax_L         = 3.1; %From similar planes: max-->3.7689, min-->2.2523, mean-->3.0764
        
        %Weight
        DP.MLW_MTOW        = 0.85;  %From SP: min-->0.7900, max-->0.9267, mean-->0.8753
        DP.MRW_MTOW        = 1.005; %From SP: min-->0.9918, max-->1.0286, mean-->1.0051
        
        %Engines
        DP.NumberEngines  = 2; 
        
        
        
        
    case 11  %11. Flying boats, amphibious, float airplanes
        %% REQUIREMENTS
       
        
        %% CUSTOM SELECTED PARAMETERS
        
        
end