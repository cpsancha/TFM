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
    
    
%% CALCULATE WINGS INCIDENCE AND STAGGER
%Create solver options
options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
                   
%First run, to remove error warnings
ME = wingsDesign(AC, ME, DP, Parameters, CST, CF);

%Solve wing incidence and stagger
[X,~,exitflag,~] = fsolve(@(X)getWingsIncidence(X, AC, ME, DP, Parameters, CST, CF),[DP.Incidence_1, DP.Incidence_2, DP.Stagger],options);
DP.Incidence_1 = X(1);
DP.Incidence_2 = X(2);
DP.Stagger     = X(3);
if ~isequal(exitflag,1)
    error('El solver que calcula las incidencias y el stagger no ha logrado converger correctamente. Se debería revisar el resultado.')
else
    clear exitflag X options
end



%% CALCULATE NECESSARY VTP
%Volume coefficient = Sv * lv / (Sw * bw)
% coefVolume = 0.065; % From similar planes, NOT USED --> 0.073


%DEFINE VTP PARAMETERS
AC.VTP.AspectRatio = DP.VTP_AspectRatio;
AC.VTP.Sweep_LE    = DP.VTP_Sweep_LE;
AC.VTP.TaperRatio  = DP.VTP_TaperRatio;
AC.VTP.t_c         = DP.VTP_t_c;
%Define airfoil
AC.VTP.Airfoil.designation = DP.VTP_Airfoil;
AC.VTP.Airfoil.wantFile    = 0;
AC.VTP.Airfoil.n           = 30;
AC.VTP.Airfoil.is_finiteTE = 0;
AC.VTP.Airfoil.HalfCosineSpacing = 1;
AC.VTP.Airfoil.Data = naca4gen(AC.VTP.Airfoil);   
   

%CRITICAL ENGINE CRITERION
run rudderGraphs;
lv = DP.VTP_X_ac - AC.Weight.x_cg;
Ye = AC.Engine.Position(2) - AC.Weight.y_cg;
kv = 1; %por no ser cola en T
kdr = interp1(kdr_drmax(:,1),kdr_drmax(:,2),DP.VTP_deltar_max);

%Torenbeek graph 9-23 pag 336
X_parameter = Ye/lv * (AC.Engine.Thrust*1e3 * DP.CLmax_TO) / ((AC.Weight.MTOW - ME.Payload)*CST.GravitySI);
Y_parameter = interp1(rudderEstimationGraph(:,1), rudderEstimationGraph(:,2), X_parameter);
Sv_S = Y_parameter / (kdr*kv*(DP.VTP_Sr_Sv*AC.VTP.AspectRatio*cosd(DP.VTP_Sweep_r))^(1/3));
S_criticalEngine = Sv_S * AC.Wing.Sw;


%CROSSWIND CRITERION
k_beta = 0.3*(AC.Weight.x_cg/AC.Fuselage.fusLength) + 0.75*(AC.Fuselage.fusHeight/AC.Fuselage.fusLength) - 0.105;
Cn_beta_f = -k_beta*((DP.VTP_Svertical*AC.Fuselage.fusLength)/(AC.Wing.Sw*AC.Wing.WingSpan))*...
             (DP.VTP_h_14/DP.VTP_h_34)^(1/2)*(DP.VTP_b_34/DP.VTP_b_14)^(1/3);
Cn_beta_i = -0.017; %High wing
VTP_volumeParameter = interp1(rudderVolumeGraph(:,1), rudderVolumeGraph(:,2), Cn_beta_f+Cn_beta_i);
S_crossWind = VTP_volumeParameter*AC.Wing.Sw*AC.Wing.WingSpan/lv;


%Decide criterion
AC.VTP.Sw        = max([S_criticalEngine,S_crossWind]);
AC.VTP.Swet      = 2*AC.VTP.Sw;
AC.VTP.WingSpan  = sqrt(AC.VTP.Sw * AC.VTP.AspectRatio);
AC.VTP.RootChord = (2/(1+AC.VTP.TaperRatio))*sqrt(AC.VTP.Sw/AC.VTP.AspectRatio);
AC.VTP.TipChord  = AC.VTP.TaperRatio*AC.VTP.RootChord;
AC.VTP.CMG       = AC.VTP.RootChord*((1+AC.VTP.TaperRatio)/2);
AC.VTP.CMA       = (2/3)*AC.VTP.RootChord*((1+AC.VTP.TaperRatio+AC.VTP.TaperRatio^2)/(1+AC.VTP.TaperRatio));
AC.VTP.Sweep_14  = atand(tand(AC.VTP.Sweep_LE)+(4/AC.VTP.AspectRatio)*((1-AC.VTP.TaperRatio)/(1+AC.VTP.TaperRatio))*(0-0.25));
AC.VTP.Sweep_12  = atand(tand(AC.VTP.Sweep_LE)+(4/AC.VTP.AspectRatio)*((1-AC.VTP.TaperRatio)/(1+AC.VTP.TaperRatio))*(0-0.50));
AC.VTP.Sweep_RE  = atand(tand(AC.VTP.Sweep_LE)+(4/AC.VTP.AspectRatio)*((1-AC.VTP.TaperRatio)/(1+AC.VTP.TaperRatio))*(0-1.00));
AC.VTP.CMA_b     = (AC.VTP.WingSpan/6)*((1+2*AC.VTP.TaperRatio)/(1+AC.VTP.TaperRatio));
AC.VTP.CMA_14    = DP.VTP_X_ac;
AC.VTP.CMA_LE    = AC.VTP.CMA_14 - AC.VTP.CMA/4;
AC.VTP.Root_LE   = AC.VTP.CMA_LE - AC.VTP.CMA_b*tand(AC.VTP.Sweep_LE);
AC.VTP.TipSweep  = AC.VTP.WingSpan/2 * tand(AC.VTP.Sweep_LE);
AC.VTP.Airfoil.rootCoordinates.xU = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.xU;
AC.VTP.Airfoil.rootCoordinates.zU = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.zU;
AC.VTP.Airfoil.rootCoordinates.xL = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.xL;
AC.VTP.Airfoil.rootCoordinates.zL = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.zL;
AC.VTP.Airfoil.tipCoordinates.xU  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.xU;
AC.VTP.Airfoil.tipCoordinates.zU  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.zU;
AC.VTP.Airfoil.tipCoordinates.xL  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.xL;
AC.VTP.Airfoil.tipCoordinates.zL  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.zL;




clear lv Ye kdr_drmax rudderEstimationGraph rudderVolumeGraph kv kdr X_parameter Y_parameter Sv_S
clear k_beta Cn_beta_f Cn_beta_i VTP_volumeParameter S_criticalEngine S_crossWind



% y_vtp = linspace(0, AC.VTP.WingSpan, 100);
% c_vtp = AC.VTP.RootChord + (AC.VTP.TipChord-AC.VTP.RootChord)/AC.VTP.WingSpan .*y_vtp;
% AC.VTP.x_ac_wf = 0.25* AC.VTP.RootChord + tan(AC.VTP.Sweep_14*pi/180)/AC.VTP.Sw *trapz(y_vtp,c_vtp.*y_vtp);
% 
% AC.VTP.Root_LE = DP.x_cg + l_v - AC.VTP.x_ac_wf;
% AC.VTP.t_c = 0.12; %NACA 0012
% 
% % Airfoil coordinates for fancy plotting in getWings
% AC.VTP.Airfoil.designation='0012';
% AC.VTP.Airfoil.wantFile = 0;
% AC.VTP.Airfoil.n=30;
% AC.VTP.Airfoil.HalfCosineSpacing=1;
% AC.VTP.Airfoil.is_finiteTE=0;
% AC.VTP.Airfoil.Points = naca4gen(AC.VTP.Airfoil);


%% PLOT LAYOUT
wingConfigurationPlotting(AC, ME, DP, Parameters, CST, CF)
   

   
    
%% USEFUL FUNCTIONS
function [] = wingConfigurationPlotting(AC, ME, DP, Parameters, CST, CF) %#ok<INUSD>
% Show Divergence Mach depending on sweep
if DP.ShowReportFigures
    sweepArray = [15,20,25,30,32.5,35,37.5];
    t_c = linspace(10,18,5);
    LegendStr = cell(0);
    options = optimoptions('fsolve',...
                           'StepTolerance',1e-9,...
                           'Display','none');
    figure()
    hold on
    for i=1:length(sweepArray)
        for j=1:length(t_c)
            [MachDiv(j),~,exitFlag,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, t_c(j)./100, sweepArray(i),AC.Wing1.CLdesign,'Supercritical'),0.8,options); %#ok<AGROW,SAGROW>
            if ~isequal(exitFlag,1)
                disp('El solver del mach de divergencia al generar la figura no ha logrado converger correctamente. Se debería revisar el resultado.')
            end
        end
        plot(t_c,MachDiv,'LineWidth',1.25,'Color',Parameters.Colors(i,:))
        LegendStr{end+1} = ['$\Lambda_{1/4}=',num2str(sweepArray(i)),'^o$']; %#ok<AGROW>
    end
    plot(AC.Wing1.Airfoil.t_c*100,AC.Wing1.MachDiv,'o','LineWidth',1.25,'Color',Parameters.Colors(i+1,:))
    LegendStr{end+1}='Design Point';
    title('$M_{dd}\ en\ funcion\ de\ la\ flecha\ y\ el\ espesor\ relativo\ del\ perfil$','interpreter','latex')
    xlabel('$t/c\ [-]$','interpreter','latex')
    ylabel('$M_{dd}\ [-]$','interpreter','latex')
    legend(LegendStr,'Location','northeast','interpreter','latex')
    legend('boxoff')
    saveFigure(ME.FiguresFolder,'SweepDecision')
    clear sweepArray t_c i j exitFlag LegendStr MachDiv options
end    
    
    


% Display Wing Lift Distribution    
if DP.ShowReportFigures
    %Wing 1
        figure()
        hold on
        [~,index] = max(AC.Wing1.clb+AC.Wing1.CLmax./AC.Wing1.CL_wf.*AC.Wing1.cla);
        plot(AC.Wing1.eta(index),AC.Wing1.clb(index)+AC.Wing1.CLmax./AC.Wing1.CL_wf.*AC.Wing1.cla(index),'o','LineWidth',1.25,'Color',Parameters.Colors(1,:))
        plot(AC.Wing1.eta,ones(1,length(AC.Wing1.eta)).*AC.Wing1.Airfoil.Cl_max*cosd(AC.Wing1.Sweep_14),'--','LineWidth',1.25,'Color',Parameters.Colors(2,:));
        plot(AC.Wing1.eta,AC.Wing1.clb+AC.Wing1.CLmax./AC.Wing1.CL_wf.*AC.Wing1.cla,'LineWidth',1.25,'Color',Parameters.Colors(3,:));
        plot(AC.Wing1.eta,AC.Wing1.clb,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
        plot(AC.Wing1.eta,AC.Wing1.CLmax./AC.Wing1.CL_wf.*AC.Wing1.cla,'LineWidth',1.25,'Color',Parameters.Colors(8,:));
        plot(AC.Wing1.eta,AC.Wing1.cla./AC.Wing1.CL_wf,'LineWidth',1.25,'Color',Parameters.Colors(4,:));
        legend('First point of stall','Maximum lift of the airfoil','Total lift distribution','Basic lift distribution','Aditional lift distribution','Aditional lift distribution for C_{L_{max}}=1','Location','southwest')
        legend('boxoff')
        xlabel('$\frac{y}{b/2}$','interpreter','latex')
        ylabel('$C_l$','interpreter','latex')
        title(['Spanwise Lift Distribution of wing 1 for $C_{Lmax}=',num2str(AC.Wing1.CLmax),'$ and $\varepsilon_t=',num2str(AC.Wing1.TipTwist),'^o$'],'interpreter','latex')
        saveFigure(ME.FiguresFolder,'SpanwiseLiftDistribution_1')
    %Wing 2
        figure()
        hold on
        [~,index] = max(AC.Wing2.clb+AC.Wing2.CLmax./AC.Wing2.CL_wf.*AC.Wing2.cla);
        plot(AC.Wing2.eta(index),AC.Wing2.clb(index)+AC.Wing2.CLmax./AC.Wing2.CL_wf.*AC.Wing2.cla(index),'o','LineWidth',1.25,'Color',Parameters.Colors(1,:))
        plot(AC.Wing2.eta,ones(1,length(AC.Wing2.eta)).*AC.Wing2.Airfoil.Cl_max*cosd(AC.Wing2.Sweep_14),'--','LineWidth',1.25,'Color',Parameters.Colors(2,:));
        plot(AC.Wing2.eta,AC.Wing2.clb+AC.Wing2.CLmax./AC.Wing2.CL_wf.*AC.Wing2.cla,'LineWidth',1.25,'Color',Parameters.Colors(3,:));
        plot(AC.Wing2.eta,AC.Wing2.clb,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
        plot(AC.Wing2.eta,AC.Wing2.CLmax./AC.Wing2.CL_wf.*AC.Wing2.cla,'LineWidth',1.25,'Color',Parameters.Colors(8,:));
        plot(AC.Wing2.eta,AC.Wing2.cla./AC.Wing2.CL_wf,'LineWidth',1.25,'Color',Parameters.Colors(4,:));
        legend('First point of stall','Maximum lift of the airfoil','Total lift distribution','Basic lift distribution','Aditional lift distribution','Aditional lift distribution for C_{L_{max}}=1','Location','southwest')
        legend('boxoff')
        xlabel('$\frac{y}{b/2}$','interpreter','latex')
        ylabel('$C_l$','interpreter','latex')
        title(['Spanwise Lift Distribution of wing 2 for $C_{Lmax}=',num2str(AC.Wing2.CLmax),'$ and $\varepsilon_t=',num2str(AC.Wing2.TipTwist),'^o$'],'interpreter','latex')
        saveFigure(ME.FiguresFolder,'SpanwiseLiftDistribution_2')
end




    
% PLOT WING LAYOUT
if DP.ShowAircraftLayout
%Import fuselage layout and tail cone layout from data file
sr=which(mfilename);
i=max(strfind(lower(sr),lower('MTORRES')))+6;
if i>0
  sr=sr(1:i);
else
  error('Cannot locate MTORRES directory. You must add the path to the MTorres directory.')
end
FuselageFile = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'fuselage.dat');
AfterbodyFile = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'tailCoordinates.dat');
[Xfus,Yfus]  = importFuselage(FuselageFile);
[Xtail,Ytail]  = importFuselage(AfterbodyFile);
clear sr i FuselageFile AfterbodyFile


%Create figure and plotting  
figure()
    hold on
    axis equal
    %COCKPIT
        plot(Xfus, Yfus,'k')
        plot(Xfus,-Yfus,'k')
    %CABIN
        plot([AC.Fuselage.ln,AC.Fuselage.ln+AC.Fuselage.cabLength],[ AC.Fuselage.fusWidth/2, AC.Fuselage.fusWidth/2],'k')
        plot([AC.Fuselage.ln,AC.Fuselage.ln+AC.Fuselage.cabLength],[-AC.Fuselage.fusWidth/2,-AC.Fuselage.fusWidth/2],'k')
    %TAIL CONE
        plot(AC.Fuselage.fusLength-AC.Fuselage.la+Xtail, Ytail,'k')
        plot(AC.Fuselage.fusLength-AC.Fuselage.la+Xtail,-Ytail,'k')
    %WING1
        %root chord
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.RootChord],[ 0, 0],'r')
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.RootChord],[-0,-0],'r')
        %tip chord
            plot([AC.Wing1.Root_LE+AC.Wing1.TipSweep,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ AC.Wing1.WingSpan/2, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.Root_LE+AC.Wing1.TipSweep,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-AC.Wing1.WingSpan/2,-AC.Wing1.WingSpan/2],'r')
        %leading edge
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.TipSweep],[ 0, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.Root_LE,AC.Wing1.Root_LE+AC.Wing1.TipSweep],[-0,-AC.Wing1.WingSpan/2],'r')
        %trailing edge
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[ 0, AC.Wing1.WingSpan/2],'r')
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord],[-0,-AC.Wing1.WingSpan/2],'r')
        %position c/4
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord/4,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[ 0, AC.Wing1.WingSpan/2],'g')
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord/4,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[-0,-AC.Wing1.WingSpan/2],'g')
        %MAC
            plot([AC.Wing1.CMA_LE,AC.Wing1.CMA_LE+AC.Wing1.CMA],[ AC.Wing1.CMA_b, AC.Wing1.CMA_b],'m')
            plot([AC.Wing1.CMA_LE,AC.Wing1.CMA_LE+AC.Wing1.CMA],[-AC.Wing1.CMA_b,-AC.Wing1.CMA_b],'m')
    %WING2
        %root chord
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.RootChord],[ 0, 0],'b')
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.RootChord],[-0,-0],'b')
        %tip chord
            plot([AC.Wing2.Root_LE+AC.Wing2.TipSweep,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[ AC.Wing2.WingSpan/2, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.Root_LE+AC.Wing2.TipSweep,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[-AC.Wing2.WingSpan/2,-AC.Wing2.WingSpan/2],'b')
        %leading edge
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.TipSweep],[ 0, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.Root_LE,AC.Wing2.Root_LE+AC.Wing2.TipSweep],[-0,-AC.Wing2.WingSpan/2],'b')
        %trailing edge
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[ 0, AC.Wing2.WingSpan/2],'b')
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord],[-0,-AC.Wing2.WingSpan/2],'b')
        %position c/4
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord/4,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[ 0, AC.Wing2.WingSpan/2],'g')
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord/4,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[-0,-AC.Wing2.WingSpan/2],'g')
        %MAC
            plot([AC.Wing2.CMA_LE,AC.Wing2.CMA_LE+AC.Wing2.CMA],[ AC.Wing2.CMA_b, AC.Wing2.CMA_b],'m')
            plot([AC.Wing2.CMA_LE,AC.Wing2.CMA_LE+AC.Wing2.CMA],[-AC.Wing2.CMA_b,-AC.Wing2.CMA_b],'m')
    %VTP
        %Root Airfoil
            plot(AC.VTP.Root_LE+AC.VTP.Airfoil.rootCoordinates.xU,AC.VTP.Airfoil.rootCoordinates.zU,'y')
            plot(AC.VTP.Root_LE+AC.VTP.Airfoil.rootCoordinates.xL,AC.VTP.Airfoil.rootCoordinates.zL,'y')
        %Tip Airfoil
            plot(AC.VTP.Root_LE+AC.VTP.TipSweep+AC.VTP.Airfoil.tipCoordinates.xU,AC.VTP.Airfoil.tipCoordinates.zU,'y')
            plot(AC.VTP.Root_LE+AC.VTP.TipSweep+AC.VTP.Airfoil.tipCoordinates.xL,AC.VTP.Airfoil.tipCoordinates.zL,'y')
    %WEIGHT
        %Center of gravity
            plot(AC.Weight.x_cg,AC.Weight.y_cg,'*')
    %ENGINES
        %Engine front
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)-AC.Engine.Length/2],...
                 [ AC.Engine.Position(2)-AC.Engine.Diameter/2,  AC.Engine.Position(2)+AC.Engine.Diameter/2],'c')
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)-AC.Engine.Length/2],...
                 [-AC.Engine.Position(2)-AC.Engine.Diameter/2, -AC.Engine.Position(2)+AC.Engine.Diameter/2],'c')
        %Engine back
            plot([ AC.Engine.Position(1)+AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [ AC.Engine.Position(2)-AC.Engine.Diameter/2,  AC.Engine.Position(2)+AC.Engine.Diameter/2],'c')
            plot([ AC.Engine.Position(1)+AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [-AC.Engine.Position(2)-AC.Engine.Diameter/2, -AC.Engine.Position(2)+AC.Engine.Diameter/2],'c')
        %Engine inside
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [ AC.Engine.Position(2)-AC.Engine.Diameter/2,  AC.Engine.Position(2)-AC.Engine.Diameter/2],'c')
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [-AC.Engine.Position(2)+AC.Engine.Diameter/2, -AC.Engine.Position(2)+AC.Engine.Diameter/2],'c')
        %Engine outside
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [ AC.Engine.Position(2)+AC.Engine.Diameter/2,  AC.Engine.Position(2)+AC.Engine.Diameter/2],'c')
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [-AC.Engine.Position(2)-AC.Engine.Diameter/2, -AC.Engine.Position(2)-AC.Engine.Diameter/2],'c')
        
    clear Xfus Yfus Xtail Ytail
end
end

function Error = getDivergenceMach(DivergenceMach,t_c_airfoil,Sweep,CLwing,AirfoilType)
    switch lower(AirfoilType)
        case 'conventional'
            Mstar = 1.00 - 0.25*CLwing*(cosd(Sweep))^-2;
        case 'high-speed'
            Mstar = 1.05 - 0.25*CLwing*(cosd(Sweep))^-2;
        case 'supercritical'
            Mstar = 1.15 - 0.25*CLwing*(cosd(Sweep))^-2;
    end
    maximum_t_c_wing = (0.3/DivergenceMach) * ...
                       ((1/(DivergenceMach*cosd(Sweep))-DivergenceMach*cosd(Sweep))^(1/3)) * ...
                       ((1-((5+(DivergenceMach*cosd(Sweep))^2)/(5+Mstar^2))^3.5)^(2/3));
    Error = t_c_airfoil - maximum_t_c_wing / cosd(Sweep);
end

function [Error, ME] = getWingsIncidence(X, AC, ME, DP, Parameters, CST, CF)
  %INPUTS:  
    %X(1) = Wing1 Incidence [º]
    %X(2) = Wing2 Incidence [º]
    %X(3) = Stagger [m]
  %OUTPUTS:  
    %Error(1) = LiftCoeff - CL0 [-]
    %Error(2) = MomentCoeff [-]
    %Error(3) = Lift wing1 - 0.7Weight [kN]
    
  %Parse inputs
    DP.Incidence_1 = X(1);
    DP.Incidence_2 = X(2);
    DP.Stagger     = X(3);
  
  %Run wing's script
    ME = wingsDesign(AC, ME, DP, Parameters, CST, CF);
    
  %Necessary calculation
    designWeight    = AC.Weight.EW + ME.Payload + AC.Weight.MFW/2;
    [rho,~,~,~,~,~] = atmos(DP.CruiseAltitude);
    CL0             = designWeight*CST.GravitySI / (0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw);
    
  %Parse outputs
    Error(1) =  AC.Wing.CL_wf - CL0;
    Error(2) =  AC.Wing.Cm_wf;
    Error(3) =  DP.Wing1_Wing2 - (AC.Wing1.Sw / AC.Wing.Sw * AC.Wing1.CL_w / CL0);
%     Error(3) = (DP.Wing1_Wing2*designWeight*CST.GravitySI - 0.5*rho*DP.CruiseSpeed^2*AC.Wing1.Sw*AC.Wing1.CL_wf)*1e-3;
    
end
