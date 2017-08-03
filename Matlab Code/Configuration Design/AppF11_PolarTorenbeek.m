% CL vector definition
CL = linspace(0.0001,20,100)

%% INDUCED DRAG
% Wing 1
eta_cp = trapz(eta, cl.*c.*eta./CL./AC.Wing1.CMG);
delta = 46.264*(eta_cp - 4/3/pi)^2;
CD.i.wing = (1 + delta).*CL.^2./(pi*AC.Wing1.AspectRatio);
CD.i.wing = CD.v.wing + 3.7e-5.*AC.Wing1.Torsion.^2;  %increment due to torsion

plot(CL,CD.v.wing)

%Wing 2
%!!!!!!! ESTUDIAR!!!!! 

% Fuselage
alpha_f = (CL-AC.Wing1.CL0)./AC.Wing1.Clalpha; %Hay que meter el ángulo de ataque nulo del avión completa -> hallar incidencias primero
CD.i.fuselage = 0.15.*alpha_f.^2.*AC.Fuselage.Volume^(2/3) ;


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