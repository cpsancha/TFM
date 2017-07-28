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


%% Wing extra parameters
DP.Stagger = 5;
AC.Wing1.LongPos = 4;
AC.Wing2.LongPos = AC.Wing1.LongPos+AC.Wing1.RootChord+DP.Stagger;


% La combinaci�n de estos parametros gobernar� el punto de entrada en
% p�rdidad:
phi = 0.3; % location of prismoidal section in fraction of semispan 
AC.Wing1.TaperRatio = 0.38; %Menor que 0.4, argumentando que el taper ratio "efectivo" es diferente respecto a un ala de simple estrechamiento
AC.Wing1.Dihedral = 0;

%% 3: Obtain derived values from geometry
AC.Wing1.WingSpan   = sqrt(AC.Wing1.AspectRatio*AC.Wing1.Sw);
AC.Wing1.CMG        = AC.Wing1.WingSpan / AC.Wing1.AspectRatio; %equals to AC.Wing1.Sw/AC.Wing1.WingSpan 
AC.Wing1.CMA        = AC.Wing1.CMG*(1 + (1+3*phi)/(3*(1-phi))*((1-AC.Wing1.TaperRatio)/((1+phi)/(1-phi)+AC.Wing1.TaperRatio))^2);

%Root and tip chord
AC.Wing1.RootChord  = (AC.Wing1.Sw/(2*AC.Wing1.WingSpan))/(phi+(1-phi)*(1+AC.Wing1.TaperRatio)/2);
AC.Wing1.TipChord   = AC.Wing1.RootChord * AC.Wing1.TaperRatio;

% Chord, leading edge and trailing edge vector coordinates
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

%% Plotting
%Longitudinal Position of the leading edge at root
   

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
figure()
    hold on
    axis equal
    %COCKPIT
        plot(Xfus, Yfus,'k')
        plot(Xfus,-Yfus,'k')
    %CABIN
        plot([Xfus(end),AC.Fuselage.fusLength],[ Yfus(end), Yfus(end)],'k')
        plot([Xfus(end),AC.Fuselage.fusLength],[-Yfus(end),-Yfus(end)],'k')
    % Wing 
        plot(AC.Wing1.LongPos+x_leadingEdge,y,'b')
        plot(AC.Wing1.LongPos+x_leadingEdge,-y,'b')
        plot(AC.Wing1.LongPos+x_trailingEdge,y,'b')
        plot(AC.Wing1.LongPos+x_trailingEdge,-y,'b')
        
        plot(AC.Wing2.LongPos+x_leadingEdge,y,'r')
        plot(AC.Wing2.LongPos+x_trailingEdge,y,'r')
        plot(AC.Wing2.LongPos+x_leadingEdge,-y,'r')       
        plot(AC.Wing2.LongPos+x_trailingEdge,-y,'r')


%%

% Sweep at each chord location 
AC.Wing1.Sweep_14 =(180/2/pi)*atan((x_leadingEdge(end)+0.25*c(end)-0.25*AC.Wing1.RootChord)/(0.5*AC.Wing1.WingSpan*(1-phi)));
AC.Wing1.Sweep_12 = 0;
AC.Wing1.Sweep_LE = (180/2/pi)*atan((AC.Wing1.RootChord - c(end))/2/(0.5*AC.Wing1.WingSpan*(1-phi)));

% Preassure center:
 x_ac = 0.25 * AC.Wing1.RootChord;
 
 
 % Espesor en la ra�z
 % Quiz�s deber�a escoger primero t/c ?
 % De esta forma y teniendo la expresi�n de la cuerda tendr�a el espesor
 % m�ximo en cada secci�n. Creo que el t/c viene fijado por el n�mero de
 % perfil que elija. Habr�a que establecer congruencia entre este criterio
 % y el overhang ratio. Por ej:
 % 1: Escojo perfil
 % 2: Compruebo el espesor en la ra�z tr= c_root*(t/c)
 % 3: Calculo que el overhang ratio est� entre los l�mites [18,22]
 AC.Wing1.RootWidth = AC.Wing1.WingSpan/2/22; %overhang ratio of 22
 
 %% Airfoil characteristics
 % Get Reynolds number at each chord location:
 [~,~,~,~,nu,~] = atmos(ME.Cruise.Altitude);
 Reynolds = ME.Cruise.Speed .* c / nu;

 
cl_alpha =(180/2/pi)*(0.991-0.15)/(8-0);  %1/rad,  NACA 63A210 pag 557 theory of wing sections
 
%% Wing lift slope
l = 0.5*AC.Wing1.WingSpan*(1-phi)/cos(AC.Wing1.Sweep_LE*2*pi/180);
E = (4*0.5*AC.Wing1.WingSpan*phi+4*l+2*AC.Wing1.TipChord)/2/AC.Wing1.WingSpan; %planform semiperimeter / wingspan

CL_alpha = (0.995 * cl_alpha / (E + cl_alpha/(pi*ME.Cruise.beta*AC.Wing1.AspectRatio))) / ME.Cruise.beta; clear E l 

%% Lift distribution
eta = y./(AC.Wing1.WingSpan/2);
x = 2*pi*AC.Wing1.AspectRatio/(cl_alpha*cos(AC.Wing1.Sweep_14*2*pi/180))







