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
    
    
AC.Wing.WingSpan      = sqrt(AC.Wing.AspectRatio*AC.Wing.Sw);
AC.Wing.CMG           = AC.Wing.Sw / AC.Wing.AspectRatio;


AC.Wing1.WingSpan     = AC.Wing.WingSpan/sqrt(2);  
AC.Wing2.WingSpan     = AC.Wing.WingSpan/sqrt(2);  

