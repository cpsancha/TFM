%% WING CONFIGURATION --> Chapter 6 & 7 Roskam; Chapter 7 Torenbeek

% Selection of wing configuration and geometric characteristics

    % Cantilever wing (without braces)
    % Leading wing --> High[ ] / Low[ ]
    % Rear wing --> High[ ] / Low[ ]
    % Zero sweep[ ] / Positive sweep[X] / Negative sweep[ ]
    % Aspect ratio
    % Thickness ratio
    % Airfoils
    % Taper ratio
    % Twist
    % Incidence angle
    % Dihedral angle
    % High lift and control surface requirements
    % Winglets


%1: Relate obtained values (Sw, AR) to a virtual global wing
%2: Relate global wing to each of the wings
AC.Wing1.Sw          = AC.Wing.Sw/2;
AC.Wing1.AspectRatio = AC.Wing.AspectRatio;


%% Wing extra parameters
%Only for Wing1
AC.Wing1.Incidence = 0;
AC.Wing1.Root_LE = 4;

%Common parameters
AC.Wing1.Torsion = -6*pi/180; 
AC.Wing1.TaperRatio = 0.45; %Menor que 0.4, argumentando que el taper ratio "efectivo" es diferente respecto a un ala de simple estrechamiento
AC.Wing1.Dihedral = 0;

%% 3: Obtain derived values from geometry
AC.Wing1.WingSpan   = sqrt(AC.Wing1.AspectRatio*AC.Wing1.Sw);
phi = AC.Engine.Position(2)/(AC.Wing1.WingSpan/2); % location of prismoidal section in fraction of semispan 
AC.Wing1.CMG        = AC.Wing1.WingSpan / AC.Wing1.AspectRatio; %equals to AC.Wing1.Sw/AC.Wing1.WingSpan 
AC.Wing1.CMA        = AC.Wing1.CMG*(1 + (1+3*phi)/(3*(1-phi))*((1-AC.Wing1.TaperRatio)/((1+phi)/(1-phi)+AC.Wing1.TaperRatio))^2);


%Root and tip chord
AC.Wing1.RootChord  = (AC.Wing1.Sw/(AC.Wing1.WingSpan))/(phi+(1-phi)*(1+AC.Wing1.TaperRatio)/2);
AC.Wing1.TipChord   = AC.Wing1.RootChord * AC.Wing1.TaperRatio;

% Chord, leading edge and trailing edge vector coordinates
y = linspace(0,AC.Wing1.WingSpan/2,100);
for i=1:length(y)
    if y(i)<phi*AC.Wing1.WingSpan/2
        c(i) = AC.Wing1.RootChord;
        x_leadingEdge(i) = 0;
        x_trailingEdge(i) = AC.Wing1.RootChord;
    else
        c(i) = AC.Wing1.RootChord + (y(i)-phi*0.5*AC.Wing1.WingSpan)*(AC.Wing1.TipChord - AC.Wing1.RootChord)*2/(1-phi)/AC.Wing1.WingSpan;
        x_leadingEdge(i)  = (AC.Wing1.RootChord-c(i))/2;
        x_trailingEdge(i) = AC.Wing1.RootChord - x_leadingEdge(i);
    end
end
AC.Wing1.TipSweep =  (AC.Wing1.RootChord-c(end))/2;


%%

% Sweep at each chord location 
AC.Wing1.Sweep_14 =(180/pi)*atan((x_leadingEdge(end)+0.25*c(end)-0.25*AC.Wing1.RootChord)/(0.5*AC.Wing1.WingSpan*(1-phi)));
AC.Wing1.Sweep_12 = 0;
AC.Wing1.Sweep_LE = (180/pi)*atan((AC.Wing1.RootChord - c(end))/2/(0.5*AC.Wing1.WingSpan*(1-phi)));

% Preassure center:
 x_ac = 0.25 * AC.Wing1.RootChord;
 
 
 % Espesor en la ra�z
 % Quiz�s deber�a escoger primero t/c ?
 % De esta forma y teniendo la expresi�n de la cuerda tendr�a el espesor
 % m�ximo en cada secci�n. Creo que el t/c viene fijado por el n�mero de
 % perfil que elija. Habr�a que establecer congruencia entre este criterio
 % y el overhang ratio. Por ej:
 % 1: Escojo perfil
 % 2: Compruebo el espesor en la ra�z tr= c_root*(t/c)
 % 3: Calculo que el overhang ratio est� entre los l�mites [18,22]
 AC.Wing1.RootWidth = AC.Wing1.WingSpan/2/22; %overhang ratio of 22
 
 %% Airfoil characteristics
 % Get Reynolds number at each chord location: (para qu�?)0
 [~,~,~,~,nu,~] = atmos(ME.Cruise.Altitude);
 Reynolds = ME.Cruise.Speed .* c / nu;

 
cl_alpha =1*(180/pi)*(0.991-0.15)/(8-0); %1/rad,  NACA 63A210 pag 557 theory of wing sections
alpha_l0 = -1.5*pi/180;
cl_max = 1.41;
c_mac = -0.25;
 
%% Wing lift slope
l = 0.5*AC.Wing1.WingSpan*(1-phi)/cos(AC.Wing1.Sweep_LE*pi/180);
E = (4*0.5*AC.Wing1.WingSpan*phi+4*l+2*AC.Wing1.TipChord)/2/AC.Wing1.WingSpan; %planform semiperimeter / wingspan

CL_alpha = (0.995 * cl_alpha / (E + cl_alpha/(pi*ME.Cruise.beta*AC.Wing1.AspectRatio))) / ME.Cruise.beta; clear l 

%% Lift distribution
eta = y./(AC.Wing1.WingSpan/2);

%x value in Fig E-5:
x = 2*pi*AC.Wing1.AspectRatio/(cl_alpha*cos(AC.Wing1.Sweep_14*pi/180)); %<--Te sobra un 2 al pasar a radianes, no? por qu� no usas cosd() y as� no te equivocas nunca?

%Load digitized points, make a 3rd order polynomic regression (better
%extrapolation than higher order(situational but who cares)) and evaluate in the point of interest:
run ('C1C2C3C4f0f30.m')

% xt = linspace(0,20,100)
C1 = polyval( polyfit(C1x,C1y,3) , x);
C2 = polyval( polyfit(C2x,C2y,3) , x);
C3 = polyval( polyfit(C3x,C3y,3) , x);
C4 = polyval( polyfit(C4x,C4y,3) , x);

% %Plot for checking regression:
% figure()
% hold on
% plot(xt, C4)
% plot(C4x,C4y)

f0 = polyval( polyfit(f0x,f0y,5) , eta); %high order, no need for extrapolation but for precission
f30 = polyval( polyfit(f30x,f30y,5) , eta);
f = f0+(f30-f0)./30.*AC.Wing1.Sweep_14;  %perform interpolation for Sweep 14
% figure()
% hold on
% plot(eta, f0)
% plot(f0x,f0y)
% plot(eta, f30)
% plot(f30x,f30y)
% plot(eta,f)

%Aditional lift distribution calculation La
La = C1.*c./AC.Wing1.CMG+(C2).*4./pi.*sqrt(1-eta.^2)+C3.*f; %<--Est� mal, no? sobra el �ltimo t�rmino, creo que has mezclado las f�rmulas E13 y E14

%Basic lift distribution Lb. Si no tienes torsion no hace falta calcularla
%porque luego se multiplica por la torsion en la punta pa k kieres saber
%eso jaja salu2
epsilon_epsilont = eta; %Torsi�n lineal desde 0 hasta la torsion en punta epsilon=epsilont*eta;
alpha_0_1 = -trapz(eta,epsilon_epsilont.*La); %Eq E-16
lambda_beta = atan((tan(AC.Wing1.Sweep_12)/ME.Cruise.beta)); %Valor nulo, flecha 1/2 nula
Lb = La.*C4.*cos(lambda_beta).*(epsilon_epsilont+alpha_0_1).*ME.Cruise.beta*E;

% figure()
% hold on
% plot(eta, La)
% plot(eta,Lb)

%cl distribution
% CL = 2*AC.Weight.MTOW*CST.GravitySI/(ME.Cruise.Density*ME.Cruise.Speed^2*AC.Wing.Sw); %Coeficiente de sustentacion orientativo en crucero
CL = 1;
cl = (La.*CL+AC.Wing1.Torsion.*cl_alpha.*Lb./E)./c.*AC.Wing1.CMG;
cla = CL.*AC.Wing1.CMG.*La./c;
clb = AC.Wing1.Torsion .* cl_alpha .* AC.Wing1.CMG .* Lb ./ (E.*c);

%% Zero lift angle
alpha_L_0_r = alpha_l0 + alpha_0_1 * AC.Wing1.Torsion; %No confundir los alphas; el primero corresponde al angulo de ataque del perfil tal que cl=0, 
% el segundo corresponde a la integral hallada en el apartado anterior
%% Maximum lift

[CLmax,index] = min((cl_max.*ones(length(eta),1)-clb')./cla');
AC.Wing1.CLmax = CLmax;

figure()
hold on

plot(eta, cla, 'b')
plot(eta, clb,'k')
plot(eta, cl, 'y','DisplayName','cl')
plot(eta, cl_max*ones(length(eta),1),'r' )
% plot(eta, (cl_max.*ones(length(eta))-clb)./cla , 'g')
plot(eta, CLmax.*ones(length(eta),1),'DisplayName','CLmax wing')

        legend('Aditional lift distribution','Basic lift distribution','Total lift distribution','clmax airfoil','CLmax wing','Location','best')
        legend('boxoff')
        xlabel('$\frac{y}{b/2}$','interpreter','latex')
        ylabel('$C_l$','interpreter','latex')
        title(['Spanwise Lift Distribution for $C_{Lmax}=',num2str(CLmax),'$ and $\varepsilon_t=',num2str(AC.Wing1.Torsion),'$'],'interpreter','latex')
        saveFigure(ME.FiguresFolder,'SpanwiseLiftDistribution')


legend('show')


%% Pitching moment

%Eq E-25
for i=1:length(y)
    if y(i)<phi*AC.Wing1.WingSpan/2
        integrando(i) = clb(i)*c(i)*y(i)*tan(0); %En el primer tramo la flecha es nula

    else
        integrando(i) = clb(i)*c(i)*y(i)*tan(AC.Wing1.Sweep_14*pi/180);
    end
end

% figure()
% plot(y,integrando)

deltaEpsilonCmac = (-2/(AC.Wing1.Sw*AC.Wing1.CMG))*trapz(y,integrando);
Cmac_w = c_mac + deltaEpsilonCmac ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From this point, values may change for the rear wing,as it doesnt share
% incidence and longitudinal position:
%% Correction1: Wing-fuselage interference
AC.Wing1.Snet = AC.Wing1.Sw - AC.Wing1.RootChord * AC.Fuselage.fusWidth;
K_I = (1 + 2.15 * AC.Fuselage.fusWidth / AC.Wing1.WingSpan ) * AC.Wing1.Snet/AC.Wing1.Sw ...
    + pi * AC.Fuselage.fusWidth^2 / (2 * CL_alpha * AC.Wing1.Sw);
K_II = (1 + 0.7 * AC.Fuselage.fusWidth / AC.Wing1.WingSpan ) * AC.Wing1.Snet/AC.Wing1.Sw;
deltazCL = -0.1 * AC.Wing1.RootChord * AC.Fuselage.fusWidth / AC.Wing1.Sw;
CL_alpha_wf = K_I * CL_alpha;

%The following expression gives the total CL of wing+ fus in function of
%fuselage angle of attack alpha_f:
% CL_wf = CL_alpha_wf *( (alpha_f - alpha_0_1*AC.Wing1.Torsion)+ K_I/K_II*(AC.Wing1.Incidence - alpha_l0))+deltazCL;

%% Correction2: Pressure center correction
deltaf1 = -1.8*AC.Wing1.Root_LE * AC.Fuselage.fusHeight * AC.Fuselage.fusWidth / (CL_alpha_wf * AC.Wing1.Sw * AC.Wing1.CMG);
%f2 correction es de dudosa aplicaci�n, no esta definida para un ala
%generica.
deltaf2 = tan(AC.Wing1.Sweep_14*pi/180)*0.273*AC.Fuselage.fusWidth*AC.Wing1.CMG*(AC.Wing1.WingSpan-AC.Fuselage.fusWidth) ...
    /((1+AC.Wing1.TaperRatio)*AC.Wing1.CMG^2*(AC.Wing1.WingSpan+2.15*AC.Fuselage.fusWidth));

x_ac_wf  = AC.Wing1.CMG * (x_ac/AC.Wing1.CMG + deltaf1+ deltaf2);
%% Correction3: Pitching moment Cmac_wf
CL0 = CL_alpha_wf *( (0 - alpha_0_1*AC.Wing1.Torsion)+ K_I/K_II*(AC.Wing1.Incidence - alpha_l0))+deltazCL;
deltaCmac = -1.8*(1 - 2.5 * AC.Fuselage.fusWidth/AC.Fuselage.fusLength)...
    *pi*AC.Fuselage.fusWidth*AC.Fuselage.fusHeight*AC.Fuselage.fusLength*CL0...
    /(4*AC.Wing1.Sw*AC.Wing1.CMG*CL_alpha_wf);   %Falta meter correcci�n por �rea no circular

%% Final values of Wing1
AC.Wing1.alpha_zeroLift = alpha_L_0_r;
AC.Wing1.x_ac_wf = x_ac_wf;
AC.Wing1.CL_wf = CL_alpha_wf *( (0 - alpha_0_1*AC.Wing1.Torsion)+ K_I/K_II*(AC.Wing1.Incidence - alpha_l0))+deltazCL;
AC.Wing1.Cm_ac_wf = Cmac_w + deltaCmac;
AC.Wing1.CL_alpha_wf = CL_alpha_wf;


%% WING2: parameters
for fn = fieldnames(AC.Wing1)'
   AC.Wing2.(fn{1}) = AC.Wing1.(fn{1});
end

DP.Stagger = 5;
AC.Wing2.Root_LE = AC.Wing1.Root_LE + AC.Wing1.RootChord + DP.Stagger;
AC.Wing2.Incidence = 0.1;

AC.Wing2.Snet = AC.Wing2.Sw - AC.Wing2.RootChord * AC.Fuselage.fusWidth;
K_I = (1 + 2.15 * AC.Fuselage.fusWidth / AC.Wing2.WingSpan ) * AC.Wing2.Snet/AC.Wing2.Sw ...
    + pi * AC.Fuselage.fusWidth^2 / (2 * CL_alpha * AC.Wing2.Sw);
K_II = (1 + 0.7 * AC.Fuselage.fusWidth / AC.Wing2.WingSpan ) * AC.Wing2.Snet/AC.Wing2.Sw;
deltazCL = -0.1 * AC.Wing2.RootChord * AC.Fuselage.fusWidth / AC.Wing2.Sw;
CL_alpha_wf = K_I * CL_alpha;

deltaf1 = -1.8*AC.Wing2.Root_LE * AC.Fuselage.fusHeight * AC.Fuselage.fusWidth / (CL_alpha_wf * AC.Wing2.Sw * AC.Wing2.CMG);
deltaf2 = tan(AC.Wing2.Sweep_14*pi/180)*0.273*AC.Fuselage.fusWidth*AC.Wing2.CMG*(AC.Wing2.WingSpan-AC.Fuselage.fusWidth) ...
    /((1+AC.Wing2.TaperRatio)*AC.Wing2.CMG^2*(AC.Wing2.WingSpan+2.15*AC.Fuselage.fusWidth));

x_ac_wf  = AC.Wing2.CMG * (x_ac/AC.Wing2.CMG + deltaf1+ deltaf2);

CL0 = CL_alpha_wf *( (0 - alpha_0_1*AC.Wing2.Torsion)+ K_I/K_II*(AC.Wing2.Incidence - alpha_l0))+deltazCL;
deltaCmac = -1.8*(1 - 2.5 * AC.Fuselage.fusWidth/AC.Fuselage.fusLength)...
    *pi*AC.Fuselage.fusWidth*AC.Fuselage.fusHeight*AC.Fuselage.fusLength*CL0...
    /(4*AC.Wing2.Sw*AC.Wing2.CMG*CL_alpha_wf);   %Falta meter correcci�n por �rea no circular

AC.Wing2.alpha_zeroLift = alpha_L_0_r;
AC.Wing2.x_ac_wf = x_ac_wf;
AC.Wing2.CL_wf = CL_alpha_wf *( (0 - alpha_0_1*AC.Wing1.Torsion)+ K_I/K_II*(AC.Wing1.Incidence - alpha_l0))+deltazCL;
AC.Wing2.Cm_ac_wf = Cmac_w + deltaCmac;
AC.Wing2.CL_alpha_wf = CL_alpha_wf;

%% Plotting
%Longitudinal Position of the leading edge at root
   

    %Import fuselage layout from data file
sr=which(mfilename);
i=max(strfind(lower(sr),lower('MTORRES')))+6;
if i>0
  sr=sr(1:i);
else
  error('Cannot locate MTORRES directory. You must add the path to the MTorres directory.')
end
FuselageFile = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'fuselage.dat');
[Xfus,Yfus]  = importFuselage(FuselageFile);
clear sr i FuselageFile
figure()
    hold on
    axis equal
    %COCKPIT
        plot(Xfus, Yfus,'k')
        plot(Xfus,-Yfus,'k')
    %CABIN
        plot([Xfus(end),AC.Fuselage.fusLength],[ Yfus(end), Yfus(end)],'k')
        plot([Xfus(end),AC.Fuselage.fusLength],[-Yfus(end),-Yfus(end)],'k')
    % Wing 
        plot(AC.Wing1.Root_LE+x_leadingEdge,y,'b')
        plot(AC.Wing1.Root_LE+x_leadingEdge,-y,'b')
        plot(AC.Wing1.Root_LE+x_trailingEdge,y,'b')
        plot(AC.Wing1.Root_LE+x_trailingEdge,-y,'b')

        %root chord
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.RootChord],[ 0, 0],'b')
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.RootChord],[-0,-0],'b')
        %tip chord
            plot([AC.Wing1.Root_LE+AC.Wing1.TipSweep,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ AC.Wing1.WingSpan/2, AC.Wing1.WingSpan/2],'b')
            plot([AC.Wing1.Root_LE+AC.Wing1.TipSweep,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-AC.Wing1.WingSpan/2,-AC.Wing1.WingSpan/2],'b')
  
fsolve(@(x)getWings(x,AC, DP, ME, Parameters), [0,0])    


