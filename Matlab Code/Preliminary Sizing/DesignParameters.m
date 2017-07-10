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
        DP.Payload     = 1250; %[kg]
        DP.Range       = 10e6; %[m]
        DP.CruiseSpeed =  225; %[m/s]
        DP.StallSpeed  =  NaN; %[m/s] 
        DP.TOFL        =  NaN; %[m]   Take-Off Field Length 
        DP.LFL         =  NaN; %[m]   Landing Field Length
        DP.ClimbRate   = 12.7; %[m/s] (2500 ft/min)
        DP.ClimbTime   =  NaN; %[s]   Time to reach cruise altitude, now is calculated from the ClimbRate and CruiseAltitude
        
    %% CUSTOM SELECTED PARAMETERS
        %Climb
        DP.ClimbSpeed          = 140;  %Climb speed in m/s (horizontal) (270 kts)
        
        %Cruise
        DP.CruiseAltitude      = 12.5e3; %[m]
        DP.CruiseEfficiency    = 14;     %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.CruiseTSFC          = 0.50;   %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Loiter
        DP.LoiterTime          = 30*60;  %[s]
        DP.LoiterEfficiency    = 16;     %[-] mean(loadFields(SP,'Actuations.L_D'),'omitnan') --> 12.25
        DP.LoiterTSFC          = 0.4;    %[lbm/(lbf·h)] mean(loadFields(SP,'Engine.TSFC'),'omitnan') --> 0.661
        
        %Low Height
        DP.LowHeightEfficiency = 12;   %Lower than CruiseEfficiency because of lower altitude (Typically lower than 10.000ft)
        DP.LowHeightTSFC       = 0.3;  %[lbm/(lbf·h)] Not sure if higher or lower than CruiseTSFC as the fly will probably be carried out below
                                       %              10.000ft and at a max of 250kts (128.611 m/s) in accordance with FAA regulations
        
        %Fuel Reserves & Alternate Airport
        DP.AlternateRange      = 370*1e3; %[m] Range to alternate airport (200 nautic miles --> 370km)
        
    case 11  %11. Flying boats, amphibious, float airplanes
end