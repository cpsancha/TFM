% Appendix F PolarTorenbeek

% CL vector definition
CL = linspace(0.03,0.03,1);

%% INDUCED DRAG
for i=1:length(CL)
    %!!!!!!!!!!!!!!!!!!!!!!!! ojo con esta adimensionalizacion!
    W(i) = ME.Cruise.q*AC.Wing.Sw*CL(i)
    options = optimoptions('fsolve','FunctionTolerance',1e-12,...
                           'StepTolerance',1e-9,...
                           'Display','none');
    x = fsolve(@(x)getTrim( x, AC, DP, ME, Parameters, AF,plotFlag , W(i)), [-0.1,0.1],options);  

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


clfront(i) = AC.Wing1.CL_wf;
clrear(i) = AC.Wing2.CL_wf;

% Fuselage
CDS.i.fuselage(i) = 0.15.*alpha_f(i).^2.*AC.Fuselage.Volume^(2/3) ;%<---Falta una superficie


% figure()
% hold on
% plot(CL, CD.i.wing1,'b')
% plot(CL, CD.i.wing2,'r')
% 
% figure()
% hold on
% plot(CL, clfront, 'b')
% plot(CL, clrear, 'r')
% plot(CL, CD.i.fuselage,'k')




%% PROFILE DRAG. WING1
    
    phi_w = 2.7*AF.t_c + 100* AF.t_c^4;
    Reynolds1 = AC.Wing1.Reynolds;
    xt_c= 0.15; %Separation point
    Cf = getCf ( Reynolds1, xt_c);   
    Cdp_min = 2*Cf*(1+phi_w);
    

    for j=1:length(Reynolds1)
        if Reynolds1(j) < 1e7
            deltaLCdp_ref(j) = (67*AF.cl_max/(log10(Reynolds1(j)))^4.5) - 0.0046*(1+2.75*AF.t_c);
        else
            deltaLCdp_ref(j) = 0.01*AF.cl_max - 0.0046*(1+2.75*AF.t_c);
        end
    end

    run generalized_Profile_Drag
    deltaLCdp = deltaLCdp_ref.* interp1(generalizedProfileDrag(:,1),generalizedProfileDrag(:,2),...
                    ((AC.Wing1.cl-AF.cli)./(AF.cl_max-AF.cli)).^2);
    cdp= deltaLCdp + Cdp_min;    
    
    %indice correspondiente a bf/2
   index =find(AC.Wing1.eta<0.5*AC.Fuselage.fusWidth/(AC.Wing1.WingSpan/2));
   index = index(end);
   
   %integral desde bf/2 hasta b/2
   integrando = AC.Wing1.c(index:end).*cdp(index:end);
    CD.p.wing1(i) = 2/AC.Wing1.Sw * trapz(AC.Wing1.eta(index:end).*AC.Wing1.WingSpan./2, integrando);
    
    
    %% Wing 2
    Reynolds1 = AC.Wing2.Reynolds;
    Cf = getCf ( Reynolds1, xt_c);   
    Cdp_min = 2*Cf*(1+phi_w);
    

    for j=1:length(Reynolds1)
        if Reynolds1(j) < 1e7
            deltaLCdp_ref(j) = (67*AF.cl_max/(log10(Reynolds1(j)))^4.5) - 0.0046*(1+2.75*AF.t_c);
        else
            deltaLCdp_ref(j) = 0.01*AF.cl_max - 0.0046*(1+2.75*AF.t_c);
        end
    end

    run generalized_Profile_Drag
    deltaLCdp = deltaLCdp_ref.* interp1(generalizedProfileDrag(:,1),generalizedProfileDrag(:,2),...
                    ((AC.Wing2.cl-AF.cli)./(AF.cl_max-AF.cli)).^2);
    cdp= deltaLCdp + Cdp_min;    
    
    %indice correspondiente a bf/2
   index =find(AC.Wing2.eta<0.5*AC.Fuselage.fusWidth/(AC.Wing2.WingSpan/2));
   index = index(end);
   
   %integral desde bf/2 hasta b/2
   integrando = AC.Wing2.c(index:end).*cdp(index:end);
    CD.p.wing2(i) = 2/AC.Wing2.Sw * trapz(AC.Wing2.eta(index:end).*AC.Wing2.WingSpan./2, integrando);
    

%% Fuselage
    [~,~,~,~,nu,~] = atmos(ME.Cruise.Altitude);
    R = ME.Cruise.Speed * AC.Fuselage.fusLength / nu;
    Cf = getCf(R, 0.05);
    Df_eff = sqrt(4/pi*AC.Fuselage.frontArea);
    lambda_eff = min([AC.Fuselage.fusLength / Df_eff, (AC.Fuselage.ln + AC.Fuselage.la) / Df_eff + 2]);
    phi_f = 2.2/lambda_eff^1.5 + 3.8/lambda_eff^3;

    CD_S_basic = Cf*AC.Fuselage.Swet*(1+phi_f); 
    deltaCD_S = AC.Fuselage.A_I * abs(sin(alpha_f(i))^3)+ AC.Fuselage.A_II*abs(sin(alpha_f(i)-AC.Fuselage.beta*pi/180)^3)/cos(AC.Fuselage.beta*pi/180);

CDS.p.fuselage(i) = CD_S_basic + deltaCD_S;


%% Engine Nacelles
R = ME.Cruise.Speed * AC.Engine.Length/ nu;
Cf = getCf ( R,0);
lambdan_eff = AC.Engine.Length / AC.Engine.Width;

%Acordarse de que hay 2 motores !
CDS.p.nacelles(i) = Cf*(1+2.2/lambdan_eff^1.5 + 3.8/lambdan_eff^3)*AC.Engine.Swet;

%% Vertical Tailplane
Cf = getCf ( ME.Cruise.Speed * AC.VTP.CMA/ nu, 0);
CDS.p.VTP(i) = 2*Cf*(1+2.75*AC.VTP.t_c*cos(AC.VTP.Sweep_12*pi/180)^2)*AC.VTP.Swet;


%% INTERFERENCE EFFECTS
%WETTED AREA CORRECTIONS
%WING/FUSELAGE INTERFERENCE
% Vortex induced:
% No attempt is made to present
% corrections for the effect of low-wing
% or high-wing positions.
% Viscous interference: F-66 and F-67, F-69
% F-66 
% CDS.p.wing_fuselage = -0.81* Cf*AC.Wing1.CL_wf*AC.Wing1.RootChord*Df; Df?
% F-69: no está indentificado D2 en la forma anterior


%NACELLES/AIRFRAME INTERFERENCE
% Negligible drag is experienced if the naoelle
% centerline coincides with the local
% chord. -> Pues perfecto
% Drag due to interference from the boundary
% layers of the nacelle and the wing can be
% taken into account by ignoring the wetted
% area correction required according to Section
% F-4.1 -> Ya está, no se ha restado el trozo de contacto con el ala

%TAILPLANE/AIRFRAME INTERFERENCE
% Coger gross área en lugar de wetted

%FIXED UNDERCARRIAGE
    %No tenemos en nuestro caso
%CANOPIES AND WINDSHIELDS (LAS VENTANILLAS DE LOS PILOTOS)
% D = 0.02* Df;
%POWERPLANT INSTALLATION
%EXTERNAL FUEL TANKS
    %No se aplica
%OTHER EFFECTS
    %Meter los porcentajes que dice
 
 
    
    
%% Streamlined Struts
% D_st = (0.015*(1+t_c)+t_c^2)*chord*length;



end
%% AUXILIAR FUNCTIONS
   function [Cf] = getCf( Reynolds, xt_c)
x=[1e6,    1e7];
y=[0.0045, 0.003];
CF0=10.^polyval(polyfit(log10(x),log10(y),1),log10(Reynolds));
   
x=[1e6,    1e7];
y = [ 0.0013, 0.00041];
CF1=10.^polyval(polyfit(log10(x),log10(y),1),log10(Reynolds));

Cf = ones(1,length(Reynolds));
for ii=1:length(Reynolds)
Cf(ii)=interp1([0,1],[CF0(ii),CF1(ii)],xt_c);
end
% loglog(linspace(1e6,1e7,length(Reynolds)),Cf)
   end