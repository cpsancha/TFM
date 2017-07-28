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
    
%1: Relate obtained values (Sw, AR) to a virtual global wing
%2: Relate global wing to each of the wings
AC.Wing1.Sw          = AC.Wing.Sw/2;
AC.Wing1.AspectRatio = AC.Wing.AspectRatio;


% Wing extra parameters

% La combinación de estos parametros gobernará el punto de entrada en
% pérdidad:
phi = 0.3; % location of prismoidal section in fraction of semispan 
AC.Wing1.TaperRatio = 0.38;
AC.Wing1.Dihedral = 0;

%3: Obtain derived values from geometry
AC.Wing1.WingSpan   = sqrt(AC.Wing1.AspectRatio*AC.Wing1.Sw);
AC.Wing1.CMG        = AC.Wing1.WingSpan / AC.Wing1.AspectRatio; %equals to AC.Wing1.Sw/AC.Wing1.WingSpan 
AC.Wing1.CMA        = AC.Wing1.CMG*(1 + (1+3*phi)/(3*(1-phi))*((1-AC.Wing1.TaperRatio)/((1+phi)/(1-phi)+AC.Wing1.TaperRatio))^2);

AC.Wing1.RootChord  = (AC.Wing1.Sw/(2*AC.Wing1.WingSpan))/(phi+(1-phi)*(1+AC.Wing1.TaperRatio)/2);
AC.Wing1.TipChord   = AC.Wing1.RootChord * AC.Wing1.TaperRatio;



%
y = linspace(0,AC.Wing1.WingSpan/2,100);
for i=1:length(y)
    if y(i)<phi*AC.Wing1.WingSpan/2
        c(i) = AC.Wing1.RootChord;
        x_leadingEdge(i) = 0;
        x_trailingEdge(i) = AC.Wing1.RootChord;
    else
        c(i) = AC.Wing1.RootChord + (y(i)-phi*0.5*AC.Wing1.WingSpan)*(AC.Wing1.TipChord - AC.Wing1.RootChord)*2/(1-phi)/AC.Wing1.WingSpan;
        x_leadingEdge(i)  = (AC.Wing1.RootChord-c(i))/2;
        x_trailingEdge(i) = AC.Wing1.RootChord - x_leadingEdge(i);
    end
end

flecha_1_4 =(180/2/pi)*atan((x_leadingEdge(end)+0.25*c(end)-0.25*AC.Wing1.RootChord)/(0.5*AC.Wing1.WingSpan*(1-phi)));
figure()
hold on
axis equal
    plot(y,x_leadingEdge(:),'k')
    plot(y,x_trailingEdge(:),'k')
% 

 x_ac = 0.25 * AC.Wing1.RootChord;
 
 % Espesor en la raíz
 
 AC.Wing1.RootWidth = AC.Wing1.WingSpan/2/22; %overhang ratio of 22




