

x1 = linspace(1000,50000, 100);
for i = 1:length(x1)
    x=x1(i);
    [ F, W_TO_guess, W_E, W_F, W_E_tent ] = getWeights( x, ME, Parameters, CF );
    y1(i)= W_E_tent;
    y2(i)=W_E;
end
plot(x1,y1,'r') ;
hold all
plot(x1,y2,'b') %regresion
%% Function definitions:

function [ F, W_TO_guess, W_E, W_F, W_E_tent ] = getWeights( x, ME, Parameters, CF )
%GETWEIGHTS: Gets the estimation of take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

%% 1. Determine the mission payload weight (W_PL)
W_PL = ME.Payload;


%% 2. Guessing a likely value of W_TO_guess:
%An initial guess is obtained by comparing the mission specification of the
%airplane with the mission capabilities of similar airplanes.
W_TO_guess = x; 


%% 3. Determination of mission fuel weight:
 %Eq 2.13
M_ff = 1; %Mission Fuel Fraction
for i=1:length(Parameters.fuelFraction(:))
    M_ff = M_ff*Parameters.fuelFraction(i).value;
end

W_F_res = 0; %NEEDS TO BE ESTABLISHED (FAR?)
W_F = (1 - M_ff)*W_TO_guess + W_F_res; %Eq 2.15


%% Step 4. Calculate a tentative value for W_OE from:
W_OE_tent = W_TO_guess - W_F - W_PL;  %Eq 2.4


%% Step 5. Calculate a tentative value for W_E from:
W_tfo = 0.005*W_TO_guess; % 0.5% of MTOW, taken from example pag.52. Note that W_tfo (trapped fuel-oil) is often neglected in this stage (page 7)
W_E_tent = W_OE_tent - W_tfo - ME.CrewWeight;   %Eq 2.4. 

%% 4. Finding the allowable value for W_E
W_E = 10^((log10(W_TO_guess)-Parameters.Table_2_15.a)/Parameters.Table_2_15.b);

% W_E = 10^((log10(W_TO_guess/CF.lbs_to_kg)-parameters.Table_2_15.a)/parameters.Table_2_15.b); %lb
% W_E = W_E*CF.lbs_to_kg;
% 
% disp(['WTO: ',num2str(W_TO_guess)])
% disp(['WF: ',num2str(W_F)])
% disp(['WE_tent: ',num2str(W_E_tent)])
% disp(['WE: ',num2str(W_E)])
F = W_E_tent-W_E;

end