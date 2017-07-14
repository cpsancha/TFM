% Estimating take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

% Get a first guess of W_TO from similar airplanes:
switch ME.MissionType
    case 5
        W_TO_guess = 30000; %[kg]
    case 11
        W_TO_guess = 40000; %[kg]
end


% Solve the iterative process:
options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
[W_TO,~,exitflag,~] = fsolve(@(x)getWeights(x,ME,CST,CF,Parameters),W_TO_guess,options);
if ~isequal(exitflag,1)
    disp('El solver no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear options W_TO_guess exitflag
end


% Get the other weights:
[~, AC.Weight.MTOW, AC.Weight.EW, AC.Weight.MFW] = getWeights( W_TO, ME, CST, CF, Parameters );

%Weight reduction of the empty weight as being fully manufatured in composite materials
switch ME.MissionType
    case 5
        AC.Weight.MTOW = AC.Weight.MTOW - (1-DP.EWnew_EWold)*AC.Weight.EW;     
        AC.Weight.EW   = DP.EWnew_EWold * AC.Weight.EW;
    case 11
end

% Define more AC weights;
switch ME.MissionType
    case 5
        AC.Weight.MRW = DP.MRW_MTOW * AC.Weight.MTOW;
        AC.Weight.MLW = DP.MLW_MTOW * AC.Weight.MTOW;
    case 11
end

clear W_TO

%Create figure showing the convergence
switch ME.MissionType
    case 5
        run Convergencia
    case 11
end


%% CONTINUE...

% run D_airplaneDesignParameters.m

%% Function definitions:

function [ F, W_TO_guess, W_E, W_F ] = getWeights( x, ME, CST, CF, Parameters )
%GETWEIGHTS: Gets the estimation of take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

%% 1. Determine the mission payload weight (W_PL)
W_PL = ME.Payload; %[kg]


%% 2. Guessing a likely value of W_TO_guess:
%An initial guess is obtained by comparing the mission specification of the
%airplane with the mission capabilities of similar airplanes.
W_TO_guess = x; %[kg]


%% 3. Determination of mission fuel weight:
 %Eq 2.13
M_ff = 1; %Mission Fuel Fraction
for i=1:length(Parameters.fuelFraction(:))
    M_ff = M_ff*Parameters.fuelFraction(i).value;
end

W_F_res = 0; %NEEDS TO BE ESTABLISHED (FAR?), en mi caso ya están incluidas, así que avisa si las vas a añadir aquí
W_F = (1 - M_ff)*W_TO_guess + W_F_res; %Eq 2.15 [kg]


%% Step 4. Calculate a tentative value for W_OE from:
W_OE_tent = W_TO_guess - W_F - W_PL;  %Eq 2.4 [kg]


%% Step 5. Calculate a tentative value for W_E from:
W_tfo = 0.005*W_TO_guess; % 0.5% of MTOW, taken from example pag.52. Note that W_tfo (trapped fuel-oil) is often neglected in this stage (page 7)
W_E_tent = W_OE_tent - W_tfo - ME.CrewWeight;   %Eq 2.4.  [kg]

%% 4. Finding the allowable value for W_E
W_E = 10^((log10(W_TO_guess*CST.GravitySI*CF.N2lbf)-Parameters.Table_2_15.a)/Parameters.Table_2_15.b); %In lbf, remember that Roskam correlation is in lbf
W_E = W_E*CF.lbf2N/CST.GravitySI; %W_E in kg

% W_E = 10^((log10(W_TO_guess)-Parameters.Table_2_15.a)/Parameters.Table_2_15.b);


% W_E = 10^((log10(W_TO_guess/CF.lbs_to_kg)-parameters.Table_2_15.a)/parameters.Table_2_15.b); %lb
% W_E = W_E*CF.lbs_to_kg;

% disp(['WTO: ',num2str(W_TO_guess)])
% disp(['WF: ',num2str(W_F)])
% disp(['WE_tent: ',num2str(W_E_tent)])
% disp(['WE: ',num2str(W_E)])
F = W_E_tent-W_E;

end