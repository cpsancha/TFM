%Apéndice C del Torenbeek


%Store cruise conditions
Backup.CruiseSpeed    = DP.CruiseSpeed;
Backup.CruiseAltitude = DP.CruiseAltitude;

%Define approach conditions
DP.CruiseAltitude = 0;
[~, a, ~, rho]    = atmosisa(ME.Cruise.Altitude);
Vapp_SP   = sqrt((loadFields(SP,'Actuations.Sl')'\(loadFields(SP,'Actuations.Vapprox').^2)').*DP.LFL);
VStall_L  = Vapp_SP/1.3; %[m/s]
% DP.CruiseSpeed    = sqrt(( 2*CST.GravitySI*AC.Weight.MLW )/( rho*AC.Wing.Sw * DP.CLmax_L))/1.3;
DP.CruiseSpeed  = VStall_L;
DP.lowSpeedFlag = true;


% Calculo el peso equivalente que correspondería a la situación de CLmaxL,
% para poder ser introducido en la misma ecuacion : 
% terminoA = Parameters.CL_max_L*AC.Wing.Sw/AC.Wing1.Sw = Weq/(ME.Cruise.q*AC.Wing1.Sw)
% Weq = 0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw*DP.CLmax_L;


%Trim aircraft
options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
[~,~,exitflag,~] = fsolve(@(X)trimAircraft(X, AC.Weight.MLW, AC, ME, DP, Parameters, CST, CF),[0, 0],options);
if ~isequal(exitflag,1)
    error('El solver del trimado no ha logrado converger correctamente. Se debería revisar el resultado.')
else
    clear exitflag
end

%Get maximum delta CL provided by flaps
deltaCLmax_1 = 1.05*(AC.Wing1.CL_wf - AC.Wing1.CLmax); %Factor 1.05 due to aditional penalyties due to trim with flags
deltaCLmax_2 = 1.05*(AC.Wing2.CL_wf - AC.Wing2.CLmax - AC.Wing2.deltaCLdeltaE); %Factor 1.05 due to aditional penalyties due to trim with flags

%Define flaps design parameters
cf_c    = DP.cf_c;    %[-]
cle_c   = DP.cle_c;   %[-]
delta_f = DP.delta_f; %[deg]

%Cociente de la superficie de flaps necesaria
Swf_Sw_1 = getSwf(deltaCLmax_1 , delta_f, cf_c, cle_c, AC.Wing1.Airfoil.Cl_alpha, AC.Wing1.Sweep_14);
[~,~,exitflag1,~] = fsolve(@(X)getEtaFlaps(X, AC, AC.Wing1, Swf_Sw_1),0.7,options);
if ~isequal(exitflag1,1)
    error('El solver del calculo de eta_inboard y eta_outboard no ha logrado converger correctamente. Se debería revisar el resultado.')
else
    clear exitflag1
end

if deltaCLmax_2 > 0
    disp('El ala 2 necesita flaps')
    Swf_Sw_2 = getSwf(deltaCLmax_2 , delta_f, cf_c, cle_c, AC.Wing2.Airfoil.Cl_alpha, AC.Wing2.Sweep_14);
    [~,~,exitflag2,~] = fsolve(@(X)getEtaFlaps(X, AC, AC.Wing2, Swf_Sw_2),0.7,options);
    if ~isequal(exitflag2,1)
        error('El solver del calculo de eta_inboard y eta_outboard no ha logrado converger correctamente. Se debería revisar el resultado.')
    else
        clear exitflag2
    end
else
    AC.Wing2.eta_inboard  = AC.Wing1.eta_inboard;
    AC.Wing2.eta_outboard = AC.Wing1.eta_inboard;
    AC.Wing2.Swf = 0;
end


%Diseño preliminar flaps --> pg 400 (Es el usado)
%Flaps in airfoils --> 8.1.2 pg 1916
%Flaps in wings    --> 8.1.4 pg 1949



%UNDO CHANGES
DP.CruiseAltitude = Backup.CruiseAltitude;
DP.CruiseSpeed    = Backup.CruiseSpeed;
DP.lowSpeedFlag   = false;

clear Backup cf_c cle_c delta_f Vapp_SP VStall_L a rho deltaCLmax_1 deltaCLmax_2 Swf_Sw_1 Swf_Sw_2 options

%% USEFUL FUNCIONS
function [Error, ME] = trimAircraft(X, AircraftWeight, AC, ME, DP, Parameters, CST, CF)
  %Parse inputs
    AC.Fuselage.fuselage_AoA = X(1);
    AC.Wing2.deltaCLdeltaE   = X(2);
    
  %Run wing's script
    ME = wingsDesign(AC, ME, DP, Parameters, CST, CF);
    
  %Necessary calculation
	[rho,~,~,~,~,~] = atmos(DP.CruiseAltitude);
 	CL0             = AircraftWeight*CST.GravitySI / (0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw); 
    
  %Parse outputs
    Error(1) =  AC.Wing.CL_wf - CL0;
    Error(2) =  AC.Wing.Cm_wf;   
end

function [Error] = getEtaFlaps(X, AC, Wing, Swf_Sw)
        %X = eta out
        
        %Expected flaps surface
        Wing.Swf = Wing.Sw * Swf_Sw;
        
        %Limits of estimated flaps surface
        y_in  = AC.Fuselage.fusWidth/2;
        y_out = X*Wing.WingSpan/2;
        Wing.eta_inboard  = y_in/(Wing.WingSpan/2);
        Wing.eta_outboard = y_out/(Wing.WingSpan/2);
        
        %Create arrays
        YArray     = linspace(y_in, y_out, 1e3);
        ChordArray = getChord(YArray, Wing.WingSpan, Wing.TaperRatio, Wing.Sw);
        
        %Estimated flaps surface
        Swf = trapz(YArray, ChordArray);
        
        Error = Wing.Swf - Swf;
end 

function [Swf_Sw] = getSwf(deltaCLmax , delta_f, cf_c, cle_c, cl_alpha, Sweep_14)
%This function obtains the total increment of CLmax in the wing using the
%following parameters: --> de hecho saca Swf_Sw (cambio de planes)
% Swf_Sw  : flap area versus wing area (ver figura)
% cf_c    : chord flap versus chord
% delta_f : flap deflection [º]
% cle_c   : chord leading edge devices
% cl_alpha: airfoil lift slope
%Basada en la función 'getDeltaCLmax'
load('alphaDeltaf.mat')
load('alphaDeltaf15.mat')
a15 = interp1(alphaDeltaf15(:,1),alphaDeltaf15(:,2),delta_f,'linear','extrap'); %#ok<NODEF>
a40 = interp1(alphaDeltaf(:,1),  alphaDeltaf(:,2),  delta_f,'linear','extrap'); %#ok<NODEF>
alpha_delta_f = interp1([0.15 , 0.4], [a15 , a40], cf_c);

load('K_Clmax_Cl.mat')
K = interp1( K_Clmax_Cl(:,1),K_Clmax_Cl(:,2), cf_c,'linear','extrap'); %#ok<NODEF>

K_Lambda    = (1-0.08*cosd(Sweep_14)^2)*cosd(Sweep_14)^(3/4);
cl_alpha_f  = cl_alpha * (1 + cf_c) * (1 + cle_c);
deltaCl     = cl_alpha_f*alpha_delta_f*deg2rad(delta_f);
deltaCl_max = K * deltaCl;
Swf_Sw      = (deltaCLmax*K_Lambda)/deltaCl_max;

end

function chord = getChord(y, SpanWidth, TaperRatio, Sw)
    chord = (2*Sw/((1+TaperRatio)*SpanWidth)).*(1-(2*(1-TaperRatio)/SpanWidth).*abs(y));
end