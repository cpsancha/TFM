% %% WING CONFIGURATION --> Chapter 6 & 7 Roskam; Chapter 7 Torenbeek
% 
% % Selection of wing configuration and geometric characteristics
% 
%     % Cantilever wing (without braces)
%     % Leading wing --> High[ ] / Low[ ]
%     % Rear wing --> High[ ] / Low[ ]
%     % Zero sweep[ ] / Positive sweep[X] / Negative sweep[ ]
%     % Aspect ratio
%     % Thickness ratio
%     % Airfoils
%     % Taper ratio
%     % Twist
%     % Incidence angle
%     % Dihedral angle
%     % High lift and control surface requirements
%     % Winglets
% 
 %% Airfoil characteristics
%  % Get Reynolds number at each chord location: (para qué?)0
%  [~,~,~,~,nu,~] = atmos(ME.Cruise.Altitude);
%  Reynolds = ME.Cruise.Speed .* c / nu;

 
AF.cl_alpha =1*(180/pi)*(0.991-0.15)/(8-0); %1/rad,  NACA 63A210 pag 557 theory of wing sections
AF.alpha_l0 = -1.5*pi/180;
AF.cl_max = 1.41;
AF.c_mac = -0.25;
AF.t_c   = 0.1;
AF.cli   = 0.005; %cl for minimum drag

%% Wing extra parameters
%Only for Wing1
AC.Wing1.Root_LE = 4;

%Common parameters
AC.Wing1.Torsion = -6*pi/180; 
AC.Wing1.TaperRatio = 0.45; %Menor que 0.4, argumentando que el taper ratio "efectivo" es diferente respecto a un ala de simple estrechamiento
AC.Wing1.Dihedral = 0;

DP.Stagger = 10;
DP.VerticalGap = 0;

%%
plotFlag = 0;
options = optimoptions('fsolve','FunctionTolerance',1e-12,...
                           'StepTolerance',1e-9,...
                           'Display','none');
%Get medium cruise weight:
CruiseFF = 1; %Mission Fuel Fraction
for i=1:4
    CruiseFF = CruiseFF*Parameters.fuelFraction(i).value;
end
W =CST.GravitySI*AC.Weight.MTOW*(CruiseFF + CruiseFF*Parameters.fuelFraction(5).value)/2;

incidences = fsolve(@(incidences)getWings(incidences,AC, DP, ME, Parameters, AF, plotFlag,W), [0.1,0.1,10],options);   

%Call function again for plotting
plotFlag = 1;
getWings(incidences,AC, DP, ME, Parameters, AF, plotFlag, W);
plotFlag = 0;



