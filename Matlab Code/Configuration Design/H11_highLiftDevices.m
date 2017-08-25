% HighLift Devices

%Most restrictive case: Landing with MLW -> now we get Vl



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
ME.Landing.Speed = Vl;
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


%Redifined cruise conditions
x = fsolve(@(x)getTrim( x, AC, DP, ME, Parameters, AF,plotFlag , Weq), [0.1,20],options); 

% Hallamos el incremento de sustentación máxima necesaria para el ala delantera
deltaCLmax = AC.Wing1.CL_wf - AC.Wing1.CLmax;

%% Flap design parameters
% Valores de los parámetros que proporcionan un incremento de sustentacion maxima
% , obtenidos tras optimizar:
delta_f_L = 44.2597;
cf_c = 0.3297;

%Otros parámetros, leading edge devices:
cle_c = 1.1;

% Obtención de la superficie de flaps necesaria para proporcionar el
% incremento de CLmax calculado anteriormente:
Swf_S = getSwf(deltaCLmax , delta_f_L, cf_c, cle_c, AF.cl_alpha);

%% Obtención del porcentaje de envergadura que corresponde a la superficie de
%flaps
AC.Wing1.Swf = Swf_S * AC.Wing1.Sw;

index =find(AC.Wing1.eta<0.5*AC.Fuselage.fusWidth/(AC.Wing1.WingSpan/2),1,'last');
eta_inboard = AC.Wing1.eta(index);
AC.Wing1.eta_inboard = eta_inboard;

j=1;
for i =(index+1):length(AC.Wing1.eta)
   %integral desde bf/2 hasta b/2
   integrando = AC.Wing1.c(index:i);
   Swf1(j) = 2*trapz(AC.Wing1.eta(index:i).*AC.Wing1.WingSpan./2, integrando);
   j=j+1;
end

AC.Wing1.eta_outboard = AC.Wing1.eta(find(Swf1<AC.Wing1.Swf,1,'last'));
   
%% Wing 2
% Hallamos el incremento de sustentación máxima necesaria para el ala
% trasera
deltaCLmax = AC.Wing2.CL_wf - AC.Wing2.CLmax;

%% Flap design parameters
% Valores de los parámetros que proporcionan un incremento de sustentacion maxima
% , obtenidos tras optimizar:
delta_f_L = 44.2597;
cf_c = 0.3297;

%Otros parámetros, leading edge devices:
cle_c = 1.1;

% Obtención de la superficie de flaps necesaria para proporcionar el
% incremento de CLmax calculado anteriormente:
Swf_S = getSwf(deltaCLmax , delta_f_L, cf_c, cle_c, AF.cl_alpha);

%% Obtención del porcentaje de envergadura que corresponde a la superficie de
%flaps
AC.Wing2.Swf = Swf_S * AC.Wing2.Sw;

index =find(AC.Wing2.eta<0.5*AC.Fuselage.fusWidth/(AC.Wing2.WingSpan/2),1,'last');
eta_inboard = AC.Wing2.eta(index);
AC.Wing2.eta_inboard = eta_inboard;

j=1;
for i =(index+1):length(AC.Wing2.eta)
   %integral desde bf/2 hasta b/2
   integrando = AC.Wing2.c(index:i);
   Swf1(j) = 2*trapz(AC.Wing2.eta(index:i).*AC.Wing2.WingSpan./2, integrando);
   j=j+1;
end

AC.Wing2.eta_outboard = AC.Wing2.eta(find(Swf1<AC.Wing2.Swf,1,'last'));
   
   %% Vuelvo a cargar las condiciones de crucero
ME.Cruise = cruiseConditions;






% F = getEtaFlaps(0.82, AC, AC.Wing1.Swf)
%     options = optimoptions('fsolve','FunctionTolerance',10,...
%                            'StepTolerance',1e-1);
% fsolve(@(x)getEtaFlaps(x, AC, AC.Wing1.Swf),0.5,options)
%% Auxiliar functions

function Swf_S = getSwf(deltaCLmax , delta_f, cf_c, cle_c, cl_alpha)
%Basada en la función de abajo, 'getDeltaCLmax'
load('alphaDeltaf.mat')
load('alphaDeltaf15.mat')
a15 = interp1(alphaDeltaf15(:,1),alphaDeltaf15(:,2),delta_f,'linear','extrap');
a40 = interp1(alphaDeltaf(:,1),alphaDeltaf(:,2),delta_f,'linear','extrap');
alpha_delta_f = interp1([0.15 , 0.4] , [a15 , a40], cf_c);

load('K_Clmax_Cl.mat')
K = interp1( K_Clmax_Cl(:,1),K_Clmax_Cl(:,2), cf_c,'linear','extrap');
Swf_S = deltaCLmax/(cle_c * cl_alpha * (1 + cf_c)* alpha_delta_f * (delta_f*pi/180) * K);

end

function deltaCLmax = getDeltaCLmax(Swf_S , delta_f, cf_c, cle_c, cl_alpha)
%This function obtains the total increment of CLmax in the wing using the
%following parameters:
% Swf_S: flap area versus wing area (ver figura)
% cf_c : chord flap versus chord
% delta_f : flap deflection [º]
% cle_c: chord leading edge devices
% cl_alpha: airfoil lift slope

load('alphaDeltaf.mat')
load('alphaDeltaf15.mat')
a15 = interp1(alphaDeltaf15(:,1),alphaDeltaf15(:,2),delta_f,'linear','extrap');
a40 = interp1(alphaDeltaf(:,1),alphaDeltaf(:,2),delta_f,'linear','extrap');
alpha_delta_f = interp1([0.15 , 0.4] , [a15 , a40], cf_c);

load('K_Clmax_Cl.mat')
K = interp1( K_Clmax_Cl(:,1),K_Clmax_Cl(:,2), cf_c,'linear','extrap');

deltaCLmax = cle_c * cl_alpha * (1 + cf_c)* alpha_delta_f * (delta_f*pi/180) * K * Swf_S;

end

function F = getEtaFlaps(x, AC, Swf)
    %indice correspondiente a bf/2
   index =find(AC.Wing1.eta<0.5*AC.Fuselage.fusWidth/(AC.Wing1.WingSpan/2),1,'last');
   eta_inboard = AC.Wing1.eta(index);
   
   eta_outboard = x;
   index_out = find(AC.Wing1.eta<eta_outboard,1,'last');
   
   %integral desde bf/2 hasta b/2
   integrando = AC.Wing1.c(index:index_out);
   Swf1 = 2*trapz(AC.Wing1.eta(index:index_out).*AC.Wing1.WingSpan./2, integrando);
    
   
   F = Swf - Swf1;
   AC.Wing1.eta_inboard = eta_inboard;
   AC.Wing1.eta_outboard = eta_outboard;
end    



