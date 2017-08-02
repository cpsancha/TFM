


%% Fuselage
% Colebrook Equation
%   f = Darcy-Weisbach friction factor
%   R = Reynolds number
%   r = relative roughness

%k = 20e-6m , 25e-6m
% r = lf/k ????? no creo

%R=lf*u/nu 

 [~,~,~,~,nu,~] = atmos(ME.Cruise.Altitude);
 R = ME.Cruise.Speed .* AC.Fuselage.fusLength / nu;

Cf=colebrook(R,r)