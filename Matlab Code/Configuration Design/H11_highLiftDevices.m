% HighLift Devices

%Most restrictive case: Landing with MLW -> now we get Vl


Parameters.CL_max_L = 3.6;

% !IMPORTANTE
% Sobreescribimos momentaneamente las condiciones de crucero por las de
% aterrizaje para poder utilizar la ecuacion del trimado 
%Guardo las cond de crucero en una nueva variable:
cruiseConditions = ME.Cruise;
%Sobreescribo:
ME.Cruise.Altitude = 0;     % in m
[~, a, ~, rho] = atmosisa(ME.Cruise.Altitude);
ME.Cruise.Density = rho;
Vl = sqrt(2*CST.GravitySI *  AC.Weight.MLW/( ME.Cruise.Density*AC.Wing.Sw * Parameters.CL_max_L));
ME.Cruise.Speed = Vl; 
ME.Cruise.Mach = ME.Cruise.Speed/a;
ME.Cruise.beta = sqrt(1-ME.Cruise.Mach^2);
ME.Cruise.Density = rho;
ME.Cruise.q = 0.5*ME.Cruise.Density*ME.Cruise.Speed^2;

clear a rho

% Calculo el peso equivalente que correspondería a la situación de CLmaxL,
% para poder ser introducido en la misma ecuacion : 
% terminoA = Parameters.CL_max_L*AC.Wing.Sw/AC.Wing1.Sw = Weq/(ME.Cruise.q*AC.Wing1.Sw)
Weq = ME.Cruise.q*AC.Wing.Sw*Parameters.CL_max_L ;

options = optimoptions('fsolve','FunctionTolerance',1e-12,...
'StepTolerance',1e-9,...
'Display','none');


%Redifine cruise conditions
x = fsolve(@(x)getTrim( x, AC, DP, ME, Parameters, AF,plotFlag , Weq), [0.1,20],options); 

% Hallemos los flaps necesarios para el ala delantera
deltaCLmax = AC.Wing1.CL_wf - AC.Wing1.CLmax;

% AC.Wing2.CL_wf


%% Vuelvo a cargar las condiciones de crucero
ME.Cruise = cruiseConditions;
