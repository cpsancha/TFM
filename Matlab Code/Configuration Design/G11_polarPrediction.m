% Appendix F PolarTorenbeek

% CL vector definition
Wvec = CST.GravitySI.*linspace(AC.Weight.MTOW, 0.6*AC.Weight.MTOW,10);
CL = Wvec./(ME.Cruise.q*AC.Wing.Sw);
CL = linspace(0,AC.Wing.CLdesign*1.5,40);
% CL = linspace(0.6631/2,0.6631/2,1);
% CL = 2.0930e+05/(ME.Cruise.q*AC.Wing.Sw);



%% INDUCED DRAG
for i=1:length(CL)
    %!!!!!!!!!!!!!!!!!!!!!!!! ojo con esta adimensionalizacion!
    W(i) = ME.Cruise.q*AC.Wing.Sw*CL(i);
    options = optimoptions('fsolve','FunctionTolerance',1e-12,...
                           'StepTolerance',1e-9,...
                       'Display','None');
    x = fsolve(@(x)getTrim( x, AC, DP, ME, Parameters, AF,plotFlag , W(i)), [0.1,20],options);  
 
clfront(i) = AC.Wing1.CL_wf;
clrear(i) = AC.Wing2.CL_wf;
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
CDS.i.fuselage(i) = 0.15.*alpha_f(i).^2.*AC.Fuselage.Volume^(2/3) ;%<---Falta una superficie

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
% ((AC.Wing1.cl-AF.cli)./(AF.cl_max-AF.cli)).^2
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
CDS.p.nacelles(i) = ME.Powerplant.Number * Cf*(1+2.2/lambdan_eff^1.5 + 3.8/lambdan_eff^3)*AC.Engine.Swet;

% %% Vertical Tailplane
Cf = getCf ( ME.Cruise.Speed * AC.VTP.CMA/ nu, 0);
CDS.p.VTP(i) = Cf*(1+2.75*AC.VTP.t_c*cos(AC.VTP.Sweep_12*pi/180)^2)*AC.VTP.Swet;


%% INTERFERENCE EFFECTS
%WETTED AREA CORRECTIONS
%WING/FUSELAGE INTERFERENCE
 fusDiameter = (AC.Fuselage.fusHeight+AC.Fuselage.fusWidth)/2; %Fuselage diameter in [m]
% Vortex induced:
% No attempt is made to present
% corrections for the effect of low-wing
% or high-wing positions.
 eta_fus_1 = fusDiameter/AC.Wing1.WingSpan;
eta_fus_2 = fusDiameter/AC.Wing2.WingSpan;
CD.inter.deltaICDv_1(i) = (0.55*eta_fus_1*(2-pi*eta_fus_1)*AC.Wing1.CL_wf^2)/((1+AC.Wing1.TaperRatio)*pi*AC.Wing1.AspectRatio);
CD.inter.deltaICDv_2(i) = (0.55*eta_fus_2*(2-pi*eta_fus_2)*AC.Wing2.CL_wf^2)/((1+AC.Wing2.TaperRatio)*pi*AC.Wing2.AspectRatio);
 
% Viscous interference: F-66 and F-67, F-69
% F-66 
 %Viscosa debida al engordamiento de la capa limite en la interseccion
        Cci_1 = 4.5*AC.Wing1.RootChord; % Length for both wing halves of the wing/fuselage intersection line, aprox 4.5 times the root chord
        Cci_2 = 4.5*AC.Wing2.RootChord; % Length for both wing halves of the wing/fuselage intersection line, aprox 4.5 times the root chord
        RootReynolds_1 = ME.Cruise.Speed * AC.Wing1.RootChord / nu;
        RootReynolds_2 = sqrt(ME.Cruise.Speed^2*Parameters.q2_qinf) * AC.Wing2.RootChord / nu;
        Cf_root_1 = getCf( RootReynolds_1, 0);
        Cf_root_2 = getCf( RootReynolds_2, 0);
        CDS.inter.deltaIDp_boundary_1(i) = 1.5*Cf_root_1*AF.t_c*AC.Wing1.RootChord*Cci_1*cosd(AC.Wing1.Sweep_12)^2;
        CDS.inter.deltaIDp_boundary_2(i) = 1.5*Cf_root_2*AF.t_c*AC.Wing2.RootChord*Cci_2*cosd(AC.Wing2.Sweep_12)^2;
%Viscosa debida al incremento de velocidad sobre la cara superior del ala frente a la parte de abajo debido a la sustentacion
        CMAReynolds_1 = ME.Cruise.Speed * AC.Wing1.RootChord / nu;
        CMAReynolds_2 = sqrt(ME.Cruise.Speed^2*Parameters.q2_qinf) * AC.Wing2.RootChord / nu;
        Cf_CMA_1 = getCf( CMAReynolds_1, 0);
        Cf_CMA_2 = getCf( CMAReynolds_2, 0);
        CDS.inter.deltaIDp_speed_1(i) = -0.81*Cf_CMA_1*AC.Wing1.CL_wf*AC.Wing1.RootChord*fusDiameter; %Valid for high wing
        CDS.inter.deltaIDp_speed_2(i) = -0.81*Cf_CMA_2*AC.Wing2.CL_wf*AC.Wing2.RootChord*fusDiameter; %Valid for high wing
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
CDS.windshieldsProtuberance(i) = 0.02*CD_S_basic;
%POWERPLANT INSTALLATION
%EXTERNAL FUEL TANKS
    %No se aplica
%OTHER EFFECTS
    %Meter los porcentajes que dice
    
%% Flat plate drag area of the Boat Hull
%Boat Hull Geometry
%b stands for Boat Hull
b =AC.Hull.Beam;
Lh = AC.Hull.Length;
rb = (b/2); %Radius of Boat Hull [m]
KA = 0.7; %Proportionality Coefficient
AHull = KA*Lh*b; %Area of Load Water Plane of Hull [m^2]
Swetb = 0.5*((pi*rb^2)+AHull+(pi*rb*Lh));%Wetted area of Boat Hull [m^2]
Qb = 1.25; %Interference Factor
Reb = (ME.Cruise.Speed*Lh*ME.Cruise.Density)/nu; %Reynolds Number
Cfb = 0.455./(log10(Reb)).^2.58;%Friction coefficient
Amaxb = (pi*(rb/2)^2)/4; %Boat Hull Cross Area [m^2]
ldb = Lh/sqrt((4/pi)*Amaxb); %Fineness ratio
Fb = 1+(60/(ldb^3))+(ldb/400); %Boat Hull Form Factor
CDS.p.Hull(i) = Cfb.*Fb.*Qb.*Swetb; %Flat plate drag area [m^2]
%% Flat Plate Drag Area of Floats
% %Float Geometry
% %f stands for Float
% nf = 2; %Number of Outriggers
% if nf == 0;
% ff = 0;
% else
bo = AC.Hull.bstabWT;
Lo = AC.Hull.LstabWT;
ro = (bo/2); %Radius [m]
AFloat = KA*Lo*bo; %Area of Load Water Plane Float [m^2]
Sexpf = (0.5*pi*ME.Cruise.Density^2)+AFloat+(pi*ME.Cruise.Density*Lo);%Float Exposed Area [m^2]
Qf = 1.5; %Interference factor
gf = Lo/bo; %Effective Fineness ratio
Ref = (ME.Cruise.Speed*Lo*ME.Cruise.Density)/nu; %Reynolds Number
Cff = 0.455./(log10(Ref)).^2.58; %Friction Coefficient
Ff = 1+(0.35/gf); %Form Factor
CDS.p.Floats(i) = Cff.*Ff.*Qf.*Sexpf.*2; %Floats Drag Area [m^2]
 
    
    
%% Streamlined Struts
% D_st = (0.015*(1+t_c)+t_c^2)*chord*length;

%% imperfections
 CDS.Imperfections(i) = 0.06 * AC.Wing1.Sw * CD.p.wing1(i) + ...
                       0.06 * Parameters.q2_qinf * AC.Wing2.Sw * CD.p.wing2(i) + ...
                       0.07 * CD_S_basic + ...
                       0.01 * (1 * AC.Wing1.Sw * CD.p.wing1(i) + ...
                               Parameters.q2_qinf * AC.Wing2.Sw * CD.p.wing2(i) + ...
                               CDS.p.fuselage(i) + ...
                               CDS.p.nacelles(i) + ...
                               CDS.p.VTP(i));

end
%% TOTAL DRAG
    DS.Induced = AC.Wing1.Sw * CD.i.wing1 + ...
                 Parameters.q2_qinf * AC.Wing2.Sw* CD.i.wing2 + ...
                 CDS.i.fuselage;
    DS.Profile = AC.Wing1.Sw * CD.p.wing1 + ...
                 Parameters.q2_qinf * AC.Wing2.Sw * CD.p.wing2 + ...
                 CDS.p.fuselage + ...
                 CDS.p.nacelles + ...
                 CDS.p.Hull + CDS.p.Floats + ...
                 CDS.p.VTP;

    DS.Interferences =  AC.Wing1.Sw * CD.inter.deltaICDv_1 + ...
                       Parameters.q2_qinf * AC.Wing2.Sw * CD.inter.deltaICDv_2 + ...
                       CDS.inter.deltaIDp_boundary_1 + CDS.inter.deltaIDp_boundary_2 + ...
                       CDS.inter.deltaIDp_speed_1 + CDS.inter.deltaIDp_speed_1+...
                       +CDS.windshieldsProtuberance;
                   
    CD.Total = (DS.Induced + DS.Profile + DS.Interferences + CDS.Imperfections)/AC.Wing.Sw;



%% Plotting

% figure()
% hold on
% plot(CL, CD.i.wing1,'b')
% plot(CL, CD.i.wing2,'r')
% % 
% figure()
% hold on
% plot(CL, clfront, 'b')
% plot(CL, clrear, 'r')
% plot(CL, CD.i.fuselage,'k')

figure()  
hold on
plot(CL, CD.Total,'k')
plot(CL, DS.Induced/AC.Wing.Sw,'b')
plot(CL, DS.Profile/AC.Wing.Sw,'r')

legend('Total Drag','Induced Drag','Profile Drag')
xlabel('C_{L} [-]')
ylabel('C_{D} [-]')
grafWidth   = 16;
grafAR      = 0.6;
set(gcf,'DefaultLineLineWidth',1.5);
set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
set(gca,'FontSize',10,'FontName','Times new Roman','box','on')


% Drag Breakdown
     clear i   
     i =  find(abs(CL-AC.Wing.CLdesign) < 0.001,1,'last');
     if isempty(i)
         i =  find(abs(CL-AC.Wing.CLdesign) < 0.01,'last');
     end

    figure()       

        
    Plot.Induced = [AC.Wing1.Sw * CD.i.wing1(i), Parameters.q2_qinf * AC.Wing2.Sw* CD.i.wing2(i), CDS.i.fuselage(i)];
    labels = {'Front Wing','Rear Wing','Fuselage'};
    ax1 = subplot(1,2,1);
    pie3(ax1,[Plot.Induced],labels)
    title(ax1,'Induced Drag');
    saveFigure(ME.FiguresFolder,'Total Drag Breakdown')
    
    ax2 = subplot(1,2,2);
    Plot.Profile = [AC.Wing1.Sw * CD.p.wing1(i), Parameters.q2_qinf * AC.Wing2.Sw * CD.p.wing2(i),  CDS.p.fuselage(i), CDS.p.nacelles(i),CDS.p.Hull(i), CDS.p.Floats(i)...
        CDS.p.VTP(i)];
     labels = {'Front Wing','Rear Wing','Fuselage','Engines','Hull','Floats','Vertical Tailplane'};
     pie3(ax2,[Plot.Profile],labels)
     title(ax2,'Profile Drag')
     
     saveFigure(ME.FiguresFolder,'Total Drag Breakdown')
     
     figure()
     Plot.Total = [DS.Induced(i), DS.Profile(i), DS.Interferences(i),CDS.Imperfections(i)];
     labels = {'Induced','Profile','Interferences', 'Imperfections'};
     pie3(Plot.Total, labels)
     saveFigure(ME.FiguresFolder,'Total Drag Breakdown')
     
     
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