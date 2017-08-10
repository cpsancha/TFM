% Appendix F PolarTorenbeek

% CL vector definition
CL = linspace(0.0001,20,100);

%% INDUCED DRAG
for i=1:length(CL)
    W(i) = ME.Cruise.q*AC.Wing1.Sw*CL(i);
    options = optimoptions('fsolve','FunctionTolerance',1e-12,...
                           'StepTolerance',1e-9,...
                           'Display','none');
    x = fsolve(@(x)getTrim( x, AC, DP, ME, Parameters, AF,plotFlag , W(i)), [0.1,0.1],options);  

alpha_f(i) = x(1);
% Wing 1
c   = AC.Wing1.c;
cl  = AC.Wing1.cl;
eta = AC.Wing1.eta;
eta_cp = trapz(eta, cl.*c.*eta./AC.Wing1.CL_wf./AC.Wing1.CMG);
delta = 46.264*(eta_cp - 4/3/pi)^2;
CD.i.wing1(i) = (1 + delta).*AC.Wing1.CL_wf.^2./(pi*AC.Wing1.AspectRatio) + 3.7e-5.*AC.Wing1.Torsion.^2;  %increment due to torsion


% Wing 2
c   = AC.Wing2.c;
cl  = AC.Wing2.cl;
eta = AC.Wing2.eta;

eta_cp = trapz(eta, cl.*c.*eta./AC.Wing2.CL_wf./AC.Wing1.CMG);
delta = 46.264*(eta_cp - 4/3/pi)^2;
CD.i.wing2(i) = (1 + delta).*AC.Wing2.CL_wf.^2./(pi*AC.Wing1.AspectRatio) + 3.7e-5.*AC.Wing2.Torsion.^2; %increment due to torsion


% Fuselage

CD.i.fuselage(i) = 0.15.*alpha_f(i).^2.*AC.Fuselage.Volume^(2/3) ;%<---Falta una superficie
end

figure()
hold on
plot(CL, CD.i.wing1,'b')
plot(CL, CD.i.wing2,'r')
plot(CL, CD.i.fuselage,'k')
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