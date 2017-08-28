

% Cruise range
W1 = CST.GravitySI * AC.Weight.MTOW; %Cruise initial weight
W2 = CST.GravitySI * (AC.Weight.EW+ME.Payload);   %Cruise final weight

Weight_cruise = linspace(W1,W2,100);
CL_cruise = Weight_cruise./(ME.Cruise.q*AC.Wing.Sw);
CD_cruise = interp1(CL,CD.Total,CL_cruise);
cp = AC.Engine.SFC*CF.c_p2SI;
% Integration (Breguet? pls.)
R = -trapz(Weight_cruise,Parameters.Cruise.n_p/cp*(CL_cruise./CD_cruise)./Weight_cruise)
E = R/ME.Cruise.Speed

%% Envolvente
%Load Factor for limit maneuvering
n = 2.1+(24e3/(10e3+AC.Weight.MTOW*CF.kg2lbm));
if n<2.5
 n = 2.5;
 elseif n>3.8
 n = 3.8;
end
%     %Load gust Factor 
mu = 2*(AC.Weight.MTOW/AC.Wing.Sw)/(ME.Cruise.Density*AC.Wing.CL_alpha_wf*AC.Wing.CMA*CST.GravitySI);
kg = 0.88*mu/(5.3+mu);
[~, asound, P, ~] = atmosisa(ME.Cruise.Altitude);
V_EAS =  correctairspeed(ME.Cruise.Speed, asound, P, 'TAS', 'EAS');
    if ME.Cruise.Altitude*CF.m2ft<20000
    U_EAS = 50*CF.ft2m;
    elseif ME.Cruise.Altitude*CF.m2ft>=20000
    U_EAS = ((50+ (ME.Cruise.Altitude*CF.m2ft-20000)*(25-50)/(50000-20000)))*CF.ft2m;
    end
ngust_V_C = 1+ kg* AC.Wing1.CL_alpha_wf*0.5*U_EAS*V_EAS*AC.Wing.Sw/(AC.Weight.MTOW*CST.GravitySI);

[~, asound, P, ~] = atmosisa(ME.Cruise.Altitude);
V_EAS =  correctairspeed(ME.Cruise.Speed/0.8, asound, P, 'TAS', 'EAS');
    if ME.Cruise.Altitude*CF.m2ft<20000
    U_EAS = 50*CF.ft2m;
    elseif ME.Cruise.Altitude*CF.m2ft>=20000
    U_EAS = ((50+ (ME.Cruise.Altitude*CF.m2ft-20000)*(25-50)/(50000-20000)))*CF.ft2m;
    end
ngust_V_D = 1+ kg* AC.Wing1.CL_alpha_wf*0.5*U_EAS*V_EAS*AC.Wing.Sw/(AC.Weight.MTOW*CST.GravitySI);


V_C = correctairspeed(ME.Cruise.Speed, asound, P, 'TAS', 'EAS');
V_D = correctairspeed(ME.Cruise.Speed/0.8, asound, P, 'TAS', 'EAS');
[~, ~, ~, rho0] = atmosisa(0);
n_vec = linspace(0,n,100);
n_vec_ = linspace(0,-1,100);
VA = (2*AC.Weight.MTOW*CST.GravitySI.*n_vec./(rho0*AC.Wing.Sw*AC.Wing1.CLmax)).^0.5;
VA_ = (2*AC.Weight.MTOW*CST.GravitySI.*1./(rho0*AC.Wing.Sw*AC.Wing1.CLmax)).^0.5;
VA_vec = linspace(0, VA_,100);
n_vec_ = -(VA_vec).^2.*rho0*AC.Wing.Sw*AC.Wing1.CLmax/(2*AC.Weight.MTOW*CST.GravitySI);
figure()
hold on
% Maniobra
plot(VA,n_vec,'k')
plot(linspace(VA(end),V_D,100),n*ones(1,100),'k')
plot(VA_vec,n_vec_,'k')
plot(linspace(VA_,V_C,100),-1*ones(1,100),'k')
plot(V_D.*ones(1,100),linspace(0,n,100),'k')
plot(linspace(V_C,V_D,100),-1+(1./(V_D-V_C)).*(linspace(V_C,V_D,100)-V_C),'k')
% Ráfaga
plot(linspace(0,V_C,100),1+(ngust_V_C-1)/V_C.*(linspace(0,V_C,100)),'k--')
plot(linspace(0,V_C,100),(1-(ngust_V_C-1)/V_C.*(linspace(0,V_C,100))),'k--')
plot(linspace(0,V_D,100),1+(ngust_V_D-1)/V_D.*(linspace(0,V_D,100)),'k--')
% plot(linspace(0,V_D,100),(1-(ngust_V_D-1)/V_D.*(linspace(0,V_D,100))),'k--')
