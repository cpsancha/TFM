% Estimating take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

%% 

% Get a first guess of W_TO from similar airplanes
 W_TO_guess = 7900; 

% Solve the iterative process:
 W_TO = fsolve(@(x)getWeights(x,ME,parameters,CF),W_TO_guess);

 % Get the other weights:
% [F, W_TO, W_E, W_F] = getWeights( W_TO_guess,ME, parameters,CF )



%% Function definitions:

function [ F, W_TO_guess, W_E, W_F ] = getWeights( x,ME, parameters,CF )
%GETWEIGHTS Summary of this function goes here
%   Detailed explanation goes here

%% 1. Determine the mission payload weight (W_PL)

W_PL = ME.Payload


%% 2. Guessing a likely value of W_TO_guess:
%An initial guess is obtained by comparing the mission specification of the
%airplane with the mission capabilities of similar airplanes.

W_TO_guess = x; 


%% 3. Determination of mission fuel weight:
 %Eq 2.13
M_ff = 1;
for i=1:length(parameters.fuelFraction(:))
    M_ff = M_ff*parameters.fuelFraction(i).value;
end
M_ff
W_F = (1 - M_ff)*W_TO_guess *1.25 %Eq 2.15


%% Step 4. Calculate a tentative value for W_OE from:

W_OE_tent = W_TO_guess - W_F - W_PL;  %Eq 2.4
W_OE_tent
%% Step 5. Calculate a tentative value for W_E from:

W_tfo = 0.005*W_TO_guess; %Note that W_tfo is often neglected in this stage (page 7)

W_E_tent = W_OE_tent-W_tfo; %Eq 2.4. 

W_E_tent

%% 4. Finding the allowable value for W_E

W_E = 10^((log10(W_TO_guess)-parameters.A)/parameters.B);
W_E
F = W_E_tent-W_E;

end