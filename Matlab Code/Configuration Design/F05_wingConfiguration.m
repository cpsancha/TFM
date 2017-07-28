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
    AC.Wing1.Sweep_12 = rad2deg(atan(tan(deg2rad(AC.Wing1.Sweep_14))+(4/AC.Wing1.AspectRatio)*((1-AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*(0.25-0.5)));
    AC.Wing2.Sweep_12 = rad2deg(atan(tan(deg2rad(AC.Wing2.Sweep_14))+(4/AC.Wing2.AspectRatio)*((1-AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*(0.25-0.5)));

%Sweep_LE
    AC.Wing1.Sweep_LE = rad2deg(atan(tan(deg2rad(AC.Wing1.Sweep_14))+(4/AC.Wing1.AspectRatio)*((1-AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*(0.25-0)));
    AC.Wing2.Sweep_LE = rad2deg(atan(tan(deg2rad(AC.Wing2.Sweep_14))+(4/AC.Wing2.AspectRatio)*((1-AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*(0.25-0)));

%Dihedral
    AC.Wing1.Dihedral = DP.Dihedral;
    AC.Wing2.Dihedral = DP.Dihedral;
 
%Longitudinal Position of the leading edge at root
    AC.Wing1.LongPos = 4;
    AC.Wing2.LongPos = AC.Wing1.LongPos+AC.Wing1.RootChord+DP.Stagger;
    
%Y_CMA
    AC.Wing1.CMA_b = (AC.Wing1.WingSpan/6)*((1+2*AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio));
    AC.Wing2.CMA_b = (AC.Wing2.WingSpan/6)*((1+2*AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio));

%X_CMA_LE
    AC.Wing1.CMA_LE = (AC.Wing1.WingSpan/6)*((1+2*AC.Wing1.TaperRatio)/(1+AC.Wing1.TaperRatio))*tan(deg2rad(AC.Wing1.Sweep_LE));
    AC.Wing2.CMA_LE = (AC.Wing2.WingSpan/6)*((1+2*AC.Wing2.TaperRatio)/(1+AC.Wing2.TaperRatio))*tan(deg2rad(AC.Wing2.Sweep_LE));

%X_Tip_LE
    AC.Wing1.TipSweep = (AC.Wing1.WingSpan/2)*tan(deg2rad(AC.Wing1.Sweep_LE));
    AC.Wing2.TipSweep = (AC.Wing2.WingSpan/2)*tan(deg2rad(AC.Wing2.Sweep_LE));


    
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
   

%PLOT WING LAYOUT  
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
    clear Xfus Yfus
     
    
    
    
 % Espesor en la raíz
 
 AC.Wing1.RootWidth = AC.Wing1.WingSpan/2/22; %overhang ratio of 22




