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

%Make sure x_ac_2 > x_cg
solveFlag = false;

%Solve wing incidence and stagger
while ~solveFlag
    [X,~,exitflag,~] = fsolve(@(X)getWingsIncidence(X, AC, ME, DP, Parameters, CST, CF),rand(1,3).*10,options);
    DP.Incidence_1 = X(1);
    DP.Incidence_2 = X(2);
    DP.Wing1_Wing2 = X(3);
    if ~isequal(exitflag,1)
        error('El solver que calcula las incidencias y el stagger no ha logrado converger correctamente. Se debería revisar el resultado.')
    else
        clear exitflag X
    end
    if AC.Wing2.x_ac_wf>DP.x_cg
        solveFlag = true;
    end
end

clear solveFlag


%% CALCULATE NECESSARY VTP
getVTP(AC, ME, DP, Parameters, CST, CF);



%% GET SOME COEFFICIENTS GRAPHS
    alphaRange = linspace(-5,10,25);
    CLRange    = zeros(1,length(alphaRange));
    CmRange    = zeros(1,length(alphaRange));
    for i=1:length(alphaRange)
        AC.Fuselage.fuselage_AoA = alphaRange(i);
        ME = wingsDesign(AC, ME, DP, Parameters, CST, CF);
        CLRange(i) = AC.Wing.CL_wf;
        CmRange(i) = AC.Wing.Cm_wf;
    end
    AC.Wing.CL_alpha_wf = deg2rad(alphaRange)'\CLRange';
    AC.Wing.Cm_alpha    = deg2rad(alphaRange)'\CmRange';
  
if DP.ShowReportFigures    
    %Show figure
    figure()
        yyaxis left
        plot(alphaRange,CLRange)
        title('Total aircraft lift and moment coefficients','interpreter','latex')
        xlabel('Fuselage angle of attack ($\alpha_{f}$) [$^o$]','interpreter','latex')
        ylabel('Lift coefficient ($C_L$) [-]','interpreter','latex')
        yyaxis right
        plot(alphaRange,CmRange)
        ylabel('Moment coefficient ($C_m$) [-]','interpreter','latex')
        saveFigure(ME.FiguresFolder,'AircraftCoefficients')
end

clear i alphaRange CLRange CmRange


%% PLOT LAYOUT
wingConfigurationPlotting(AC, ME, DP, Parameters, CST, CF);
   


%% GET STICK-FIXED NEUTRAL POINT
[AC.Weight.x_n,~,exitflag,~] = fsolve(@(Xn)getNeutralPoint(Xn, AC, ME, DP, Parameters, CST, CF),DP.x_cg,options);
if ~isequal(exitflag,1)
    error('El solver que calcula las incidencias y el stagger no ha logrado converger correctamente. Se debería revisar el resultado.')
else
    clear exitflag options
end


   
    
%% USEFUL FUNCTIONS
function [] = wingConfigurationPlotting(AC, ME, DP, Parameters, CST, CF) %#ok<INUSD>
% Show Divergence Mach depending on sweep
if DP.ShowReportFigures
    sweepArray = [15,20,25,27.5,30,32.5,35];
    t_c = linspace(10,16,5);
    LegendStr = cell(0);
    options = optimoptions('fsolve',...
                           'StepTolerance',1e-9,...
                           'Display','none');
    figure()
    hold on
    plot(AC.Wing1.Airfoil.t_c*100,AC.Wing1.MachDiv,'o','LineWidth',1.25,'Color',Parameters.Colors(1,:))
    LegendStr{end+1}='Design Point';
    for i=1:length(sweepArray)
        for j=1:length(t_c)
            [MachDiv(j),~,exitFlag,~] = fsolve(@(DivergenceMach)getDivergenceMach(DivergenceMach, t_c(j)./100, sweepArray(i),AC.Wing1.CLdesign,'Supercritical'),0.8,options); %#ok<AGROW>
            if ~isequal(exitFlag,1)
                disp('El solver del mach de divergencia al generar la figura no ha logrado converger correctamente. Se debería revisar el resultado.')
            end
        end
        plot(t_c,MachDiv,'LineWidth',1.25,'Color',Parameters.Colors(i+1,:))
        LegendStr{end+1} = ['$\Lambda_{1/4}=',num2str(sweepArray(i)),'^o$']; %#ok<AGROW>
    end
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
FuselageFile    = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'fuselage.dat');
AfterbodyFile   = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'tailCoordinates.dat');
VertNoseFile    = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'verticalNose.dat');
VertAftFile     = fullfile(sr,'Matlab Code',filesep,'Digitalized Data',filesep,'verticalAfterbody.dat');
[Xfus,Yfus]     = importFuselage(FuselageFile);
[Xtail,Ytail]   = importFuselage(AfterbodyFile);
[XvNose,YvNose] = importFuselage(VertNoseFile);
[XvAft,YvAft]   = importFuselage(VertAftFile);


clear sr i FuselageFile AfterbodyFile VertNoseFile VertAftFile


%Create layout figure and plotting  
figure()
subplot(3,1,[1 2])
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
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord/4,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[ 0, AC.Wing1.WingSpan/2],'y')
            plot([AC.Wing1.Root_LE+AC.Wing1.RootChord/4,AC.Wing1.Root_LE+AC.Wing1.TipSweep+AC.Wing1.TipChord/4],[-0,-AC.Wing1.WingSpan/2],'y')
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
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord/4,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[ 0, AC.Wing2.WingSpan/2],'y')
            plot([AC.Wing2.Root_LE+AC.Wing2.RootChord/4,AC.Wing2.Root_LE+AC.Wing2.TipSweep+AC.Wing2.TipChord/4],[-0,-AC.Wing2.WingSpan/2],'y')
        %MAC
            plot([AC.Wing2.CMA_LE,AC.Wing2.CMA_LE+AC.Wing2.CMA],[ AC.Wing2.CMA_b, AC.Wing2.CMA_b],'m')
            plot([AC.Wing2.CMA_LE,AC.Wing2.CMA_LE+AC.Wing2.CMA],[-AC.Wing2.CMA_b,-AC.Wing2.CMA_b],'m')
    %VTP
        %Root Airfoil
            plot(AC.VTP.Root_LE+AC.VTP.Airfoil.rootCoordinates.xU,AC.VTP.Airfoil.rootCoordinates.zU,'g')
            plot(AC.VTP.Root_LE+AC.VTP.Airfoil.rootCoordinates.xL,AC.VTP.Airfoil.rootCoordinates.zL,'g')
        %Tip Airfoil
            plot(AC.VTP.Root_LE+AC.VTP.TipSweep+AC.VTP.Airfoil.tipCoordinates.xU,AC.VTP.Airfoil.tipCoordinates.zU,'g')
            plot(AC.VTP.Root_LE+AC.VTP.TipSweep+AC.VTP.Airfoil.tipCoordinates.xL,AC.VTP.Airfoil.tipCoordinates.zL,'g')
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
    
    YL = ylim;
    ylim([YL(1)-0.25,YL(2)+0.25]);
    XL = xlim;
    clear Xfus Yfus Xtail Ytail YL
    
    
%Create vertical figure and plotting  
% figure()
subplot(3,1,3)
    hold on
    axis equal
    %VERTICAL NOSE
        plot(XvNose, YvNose,'k')
    %VERTICAL CABIN
        plot([AC.Fuselage.ln,AC.Fuselage.fusLength-AC.Fuselage.la],[AC.Fuselage.fusHeight,AC.Fuselage.fusHeight],'k')
        plot([AC.Fuselage.ln,AC.Fuselage.fusLength-AC.Fuselage.la],[0,0],'k')
    %VERTICAL AFTERBODY
        plot(XvAft, YvAft,'k')
     %ENGINES
        %Engine front
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)-AC.Engine.Length/2],...
                 [ AC.Engine.Position(3)-AC.Engine.Diameter/2,  AC.Engine.Position(3)+AC.Engine.Diameter/2],'c')
        %Engine back
            plot([ AC.Engine.Position(1)+AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [ AC.Engine.Position(3)-AC.Engine.Diameter/2,  AC.Engine.Position(3)+AC.Engine.Diameter/2],'c')
        %Engine up
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [ AC.Engine.Position(3)+AC.Engine.Diameter/2,  AC.Engine.Position(3)+AC.Engine.Diameter/2],'c')
        %Engine down
            plot([ AC.Engine.Position(1)-AC.Engine.Length/2,    AC.Engine.Position(1)+AC.Engine.Length/2],...
                 [ AC.Engine.Position(3)-AC.Engine.Diameter/2,  AC.Engine.Position(3)-AC.Engine.Diameter/2],'c')
    %VERTICAL VTP
        h_LE = interp1(XvAft(1:18),YvAft(1:18),AC.VTP.Root_LE);
        h_RE = interp1(XvAft(1:18),YvAft(1:18),AC.VTP.Root_LE+AC.VTP.RootChord);
        plot([AC.VTP.Root_LE,AC.VTP.Root_LE+AC.VTP.TipSweep],[h_LE,h_LE+AC.VTP.WingSpan],'g')
        plot([AC.VTP.Root_LE+AC.VTP.TipSweep,AC.VTP.Root_LE+AC.VTP.TipSweep+AC.VTP.TipChord],...
             [h_LE+AC.VTP.WingSpan,h_LE+AC.VTP.WingSpan],'g')
        plot([AC.VTP.Root_LE+AC.VTP.TipSweep+AC.VTP.TipChord,AC.VTP.Root_LE+AC.VTP.RootChord],...
             [h_LE+AC.VTP.WingSpan,h_RE],'g')
        plot([AC.VTP.CMA_LE,AC.VTP.CMA_LE+AC.VTP.CMA],[h_LE+AC.VTP.CMA_b,h_LE+AC.VTP.CMA_b],'m')
        plot(AC.VTP.CMA_14,h_LE+AC.VTP.CMA_b,'*g')
    %WING 1
        %Root
        plot(AC.Wing1.Root_LE + AC.Wing1.RootChord.*AC.Wing1.Airfoil.Coordinates.xU,...
             AC.Fuselage.fusHeight+AC.Wing1.RootChord.*AC.Wing1.Airfoil.Coordinates.zU,'r')
        plot(AC.Wing1.Root_LE + AC.Wing1.RootChord.*AC.Wing1.Airfoil.Coordinates.xL,...
             AC.Fuselage.fusHeight+AC.Wing1.RootChord.*AC.Wing1.Airfoil.Coordinates.zL,'r')
        %Tip
        RotationMatrix = [cosd(-AC.Wing1.TipTwist) -sind(-AC.Wing1.TipTwist); sind(-AC.Wing1.TipTwist) cosd(-AC.Wing1.TipTwist)];
        points=RotationMatrix*[AC.Wing1.Airfoil.Coordinates.xU-0.5 AC.Wing1.Airfoil.Coordinates.zU]';
        AC.Wing1.Airfoil.Coordinates.xU = points(1,:)'+0.5;
        AC.Wing1.Airfoil.Coordinates.zU = points(2,:)';
        points=RotationMatrix*[AC.Wing1.Airfoil.Coordinates.xL-0.5 AC.Wing1.Airfoil.Coordinates.zL]';
        AC.Wing1.Airfoil.Coordinates.xL = points(1,:)'+0.5;
        AC.Wing1.Airfoil.Coordinates.zL = points(2,:)';
        plot(AC.Wing1.Root_LE + AC.Wing1.TipSweep + AC.Wing1.TipChord.*AC.Wing1.Airfoil.Coordinates.xU,...
             AC.Fuselage.fusHeight+AC.Wing1.TipChord.*AC.Wing1.Airfoil.Coordinates.zU,'r')
        plot(AC.Wing1.Root_LE + AC.Wing1.TipSweep + AC.Wing1.TipChord.*AC.Wing1.Airfoil.Coordinates.xL,...
             AC.Fuselage.fusHeight+AC.Wing1.TipChord.*AC.Wing1.Airfoil.Coordinates.zL,'r')
    %WING 2
        %Root
        plot(AC.Wing2.Root_LE + AC.Wing2.RootChord.*AC.Wing2.Airfoil.Coordinates.xU,...
             AC.Fuselage.fusHeight+DP.VerticalGap+AC.Wing2.RootChord.*AC.Wing2.Airfoil.Coordinates.zU,'b')
        plot(AC.Wing2.Root_LE + AC.Wing2.RootChord.*AC.Wing2.Airfoil.Coordinates.xL,...
             AC.Fuselage.fusHeight+DP.VerticalGap+AC.Wing2.RootChord.*AC.Wing2.Airfoil.Coordinates.zL,'b')
        %Tip
        RotationMatrix = [cosd(-AC.Wing2.TipTwist) -sind(-AC.Wing2.TipTwist); sind(-AC.Wing2.TipTwist) cosd(-AC.Wing2.TipTwist)];
        points=RotationMatrix*[AC.Wing2.Airfoil.Coordinates.xU-0.5 AC.Wing2.Airfoil.Coordinates.zU]';
        AC.Wing2.Airfoil.Coordinates.xU = points(1,:)'+0.5;
        AC.Wing2.Airfoil.Coordinates.zU = points(2,:)';
        points=RotationMatrix*[AC.Wing2.Airfoil.Coordinates.xL-0.5 AC.Wing2.Airfoil.Coordinates.zL]';
        AC.Wing2.Airfoil.Coordinates.xL = points(1,:)'+0.5;
        AC.Wing2.Airfoil.Coordinates.zL = points(2,:)';
        plot(AC.Wing2.Root_LE + AC.Wing2.TipSweep + AC.Wing2.TipChord.*AC.Wing2.Airfoil.Coordinates.xU,...
             AC.Fuselage.fusHeight+AC.Wing2.TipChord.*AC.Wing2.Airfoil.Coordinates.zU,'b')
        plot(AC.Wing2.Root_LE + AC.Wing2.TipSweep + AC.Wing2.TipChord.*AC.Wing1.Airfoil.Coordinates.xL,...
             AC.Fuselage.fusHeight+AC.Wing2.TipChord.*AC.Wing2.Airfoil.Coordinates.zL,'b')

   
   xlim(XL);
   saveFigure(ME.FiguresFolder,'AircraftLayout')
   clear XvNose YvNose XvAft YvAft h_LE h_RE RotationMatrix points XL
end
end

function [] = getVTP(AC, ME, DP, Parameters, CST, CF) %#ok<INUSL,INUSD>

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
kdr = interp1(kdr_drmax(:,1),kdr_drmax(:,2),DP.VTP_deltar_max); %#ok<NODEF>

%Torenbeek graph 9-23 pag 336
X_parameter = Ye/lv * (AC.Engine.Thrust*1e3 * DP.CLmax_TO) / ((AC.Weight.MTOW - ME.Payload)*CST.GravitySI);
Y_parameter = interp1(rudderEstimationGraph(:,1), rudderEstimationGraph(:,2), X_parameter); %#ok<NODEF>
Sv_S = Y_parameter / (kdr*kv*(DP.VTP_Sr_Sv*AC.VTP.AspectRatio*cosd(DP.VTP_Sweep_r))^(1/3));
S_criticalEngine = Sv_S * AC.Wing.Sw;


%STABILITY CRITERION
k_beta = 0.3*(AC.Weight.x_cg/AC.Fuselage.fusLength) + 0.75*(AC.Fuselage.fusHeight/AC.Fuselage.fusLength) - 0.105;
Cn_beta_f = -k_beta*((DP.VTP_Svertical*AC.Fuselage.fusLength)/(AC.Wing.Sw*AC.Wing.WingSpan))*...
             (DP.VTP_h_14/DP.VTP_h_34)^(1/2)*(DP.VTP_b_34/DP.VTP_b_14)^(1/3);
Cn_beta_i = -0.017; %High wing
VTP_volumeParameter = interp1(rudderVolumeGraph(:,1), rudderVolumeGraph(:,2), Cn_beta_f+Cn_beta_i); %#ok<NODEF>
S_stability = VTP_volumeParameter*AC.Wing.Sw*AC.Wing.WingSpan/lv;
S_stability = 0.55 * S_stability; %Relaxation due to FCS

%Decide criterion
AC.VTP.Sw        = 0.75*max([S_criticalEngine,S_stability]);  %AREA
AC.VTP.Swet      = 2*AC.VTP.Sw;                          %WET AREA
AC.VTP.WingSpan  = sqrt(AC.VTP.Sw * AC.VTP.AspectRatio); %HEIGHT
AC.VTP.CMG       = AC.VTP.Sw/AC.VTP.WingSpan;            %CMG
AC.VTP.RootChord = 2*AC.VTP.CMG/(1+AC.VTP.TaperRatio);   %C_r
AC.VTP.TipChord  = AC.VTP.TaperRatio*AC.VTP.RootChord;   %C_t
AC.VTP.CMA       = (2/3)*AC.VTP.RootChord*((1+AC.VTP.TaperRatio+AC.VTP.TaperRatio^2)/(1+AC.VTP.TaperRatio)); %CMA
AC.VTP.CMA_b     = (AC.VTP.WingSpan/3)*((1+2*AC.VTP.TaperRatio)/(1+AC.VTP.TaperRatio));
AC.VTP.CMA_14    = DP.VTP_X_ac;
AC.VTP.CMA_LE    = AC.VTP.CMA_14 - AC.VTP.CMA/4;
AC.VTP.Root_LE   = AC.VTP.CMA_LE - AC.VTP.CMA_b*tand(AC.VTP.Sweep_LE);
AC.VTP.TipSweep  = AC.VTP.WingSpan * tand(AC.VTP.Sweep_LE);
AC.VTP.Sweep_14  = atand((AC.VTP.TipSweep+(AC.VTP.TipChord-AC.VTP.RootChord)/4)/AC.VTP.WingSpan);
AC.VTP.Sweep_12  = atand((AC.VTP.TipSweep+(AC.VTP.TipChord-AC.VTP.RootChord)/2)/AC.VTP.WingSpan);
AC.VTP.Sweep_RE  = atand((AC.VTP.TipSweep+(AC.VTP.TipChord-AC.VTP.RootChord)/1)/AC.VTP.WingSpan);



AC.VTP.Airfoil.rootCoordinates.xU = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.xU;
AC.VTP.Airfoil.rootCoordinates.zU = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.zU;
AC.VTP.Airfoil.rootCoordinates.xL = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.xL;
AC.VTP.Airfoil.rootCoordinates.zL = AC.VTP.RootChord .* AC.VTP.Airfoil.Data.zL;
AC.VTP.Airfoil.tipCoordinates.xU  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.xU;
AC.VTP.Airfoil.tipCoordinates.zU  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.zU;
AC.VTP.Airfoil.tipCoordinates.xL  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.xL;
AC.VTP.Airfoil.tipCoordinates.zL  = AC.VTP.TipChord  .* AC.VTP.Airfoil.Data.zL;




clear lv Ye kdr_drmax rudderEstimationGraph rudderVolumeGraph kv kdr X_parameter Y_parameter Sv_S
clear k_beta Cn_beta_f Cn_beta_i VTP_volumeParameter S_criticalEngine S_stability

end

function [Error] = getDivergenceMach(DivergenceMach,t_c_airfoil,Sweep,CLwing,AirfoilType)
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

function [Error] = getNeutralPoint(Xn, AC, ME, DP, Parameters, CST, CF)

    AC.Fuselage.fuselage_AoA = 0;
    AC.Wing2.deltaCLdeltaE   = 0;
    
    %Run wing's script
    wingsDesign(AC, ME, DP, Parameters, CST, CF);
    
    Error = Parameters.q1_qinf*AC.Wing1.Sw/AC.Wing.Sw*(Xn-AC.Wing1.x_ac_wf)/AC.Wing.CMA*AC.Wing1.CL_alpha_wf + ...
            Parameters.q2_qinf*AC.Wing2.Sw/AC.Wing.Sw*(Xn-AC.Wing2.x_ac_wf)/AC.Wing.CMA*AC.Wing2.CL_alpha_wf * ...
            (1-AC.Wing2.deltaE_deltaAlpha);
    
end

function [Error, ME] = getWingsIncidence(X, AC, ME, DP, Parameters, CST, CF)
  %INPUTS:  
    %X(1) = Wing1 Incidence [º]
    %X(2) = Wing2 Incidence [º]
    %X(3) = Wing1_Wing2 [-]
  %OUTPUTS:  
    %Error(1) = LiftCoeff - CL0 [-]
    %Error(2) = MomentCoeff [-]
    %Error(3) = Lift wing1 - 0.7Weight [kN]
    
  %Parse inputs
    DP.Incidence_1 = X(1);
    DP.Incidence_2 = X(2);
    DP.Wing1_Wing2 = X(3);
  
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
