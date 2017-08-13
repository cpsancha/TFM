% Estimating take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

% Get a first guess of W_TO from similar airplanes:
W_TO_guess = 30000; %[kg]


% Solve the iterative process without weight reduction due to application of composites on fuselage:
options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
[MTOW.oldRoskam,~,exitflag,~] = fsolve(@(x)getWeights(x,ME,CST,CF,Parameters,'EW_old'),W_TO_guess,options);
[~, MTOW.oldRoskam, EW.oldRoskam, MFW.oldRoskam] = getWeights( MTOW.oldRoskam, ME, CST, CF, Parameters, 'EW_old');
if ~isequal(exitflag,1)
    disp('El solver del MTOW por Roskam no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear W_TO_guess exitflag
end

% Solve the iterative process with weight reduction due to application of composites on fuselage:
[MTOW.newRoskam,~,exitflag,~] = fsolve(@(x)getWeights(x,ME,CST,CF,Parameters,'EW_new'),MTOW.oldRoskam,options);
[~, MTOW.newRoskam, EW.newRoskam, MFW.newRoskam] = getWeights( MTOW.newRoskam, ME, CST, CF, Parameters, 'EW_new');
if ~isequal(exitflag,1)
    disp('El solver del MTOW por Roskam no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear exitflag
end


%Create figure showing the convergence
if DP.ShowReportFigures
    % EXPECTED TO BE RUN AFTER B_loadParameters, IF NOT, COMMENT THIS:
    h1=gcf;
    h2=figure();
    objects=allchild(h1); %#ok<NASGU>
    copyobj(get(h1,'children'),h2);
    
    %Define linspace values of MTOW
    x=linspace(18.6e3,45e3,30);
    for i = 1:length(x)
        [ ~, ~, W_E, ~, W_E_tent] = getWeights(x(i), ME, CST, CF, Parameters, 'EW_new');
        y1(i) = W_E_tent; %#ok<SAGROW> %Fuel fraction method
        y2(i) = W_E;      %#ok<SAGROW> %New EW regresion, with the weight reduction corrections
        
    end
    plot(y2,x,'LineWidth',1.25,'Color',Parameters.Colors(3,:)); %<-- EW regresion with weight reduction corrections
    plot(y1,x,'LineWidth',1.25,'Color',Parameters.Colors(6,:)); %<-- EW calculated with the fuel fraction method
    %EW_old intersection point
    plot(EW.oldRoskam,MTOW.oldRoskam,'+','LineWidth',2,'Color',Parameters.Colors(5,:));
    txt=['$$EW_{old}:\ $$',num2str(round(EW.oldRoskam)),'kg\ \ $$\rightarrow\ \ \ $$'];
    text(EW.oldRoskam,MTOW.oldRoskam,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
    txt=['$$MTOW_{old}:\ $$',num2str(round(MTOW.oldRoskam)),'kg\ \ \ '];
    text(EW.oldRoskam,MTOW.oldRoskam+2e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
    %EW_new intersection point
    plot(EW.newRoskam,MTOW.newRoskam,'+','LineWidth',2,'Color',Parameters.Colors(4,:));
    txt=['EW:\ ',num2str(round(EW.newRoskam)),'kg\ \ $$\rightarrow\ \ \ $$'];
    text(EW.newRoskam,MTOW.newRoskam,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
    txt=['MTOW:\ ',num2str(round(MTOW.newRoskam)),'kg\ \ \ '];
    text(EW.newRoskam,MTOW.newRoskam+1.5e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
    p=findobj(gca,'Type','text');
    for i=1:length(p)
        uistack(p(i), 'top')
    end
    %Legend
    h = findobj(gca,'Type','line');
    legend([h(5),h(4),h(3)],{'Similar Planes','Similar Planes (Weight correction)','Fuel-Fraction Method'},'Location','southeast')
    legend('boxoff')
    saveFigure(ME.FiguresFolder,'MTOW definition')
    clear h h1 h2 x y1 y2 W_E W_E_tent i objects txt
end
clear MTOW.oldRoskam MTOW.newRoskam EW.oldRoskam p



%% TORENBEEK GUESSTIMATION
Wfix = 500;

%Load fields from similar plane and create regresions
Wf_Wto      =  loadFields(SP,'Weight.MFW')./loadFields(SP,'Weight.MTOW');
[T,~,P,~]   =  atmosisa(loadFields(SP,'Actuations.Hcruise')); %#ok<ASGLU>
[T0,a0,~,~] =  atmosisa(0); %#ok<ASGLU>
parentesis  =  (1./(loadFields(SP,'Actuations.Mcruise').*sqrt(loadFields(SP,'Wing.AspectRatio')))) + ...
               (0.068.*P.*loadFields(SP,'Actuations.Mcruise').*loadFields(SP,'Fuselage.fusLength').*...
               (loadFields(SP,'Fuselage.fusWidth')+loadFields(SP,'Fuselage.fusHeight'))./(2.*CST.GravitySI.*loadFields(SP,'Weight.MTOW')));
reserveRange = 0.75.*CF.hour2sec.*loadFields(SP,'Actuations.Vcruise'); %Range for reserves in meters     
parametro    = (loadFields(SP,'Actuations.Range').*1e3+reserveRange).*loadFields(SP,'Engine.TSFC').*CF.TSFC2SI.*CST.GravitySI.*parentesis./(a0);%.*sqrt(T./T0));
indexPar     = ~isnan(parametro);
[Wf_Wto_fit,gof] = fit(parametro(indexPar)',Wf_Wto(indexPar)','a*x^2+b*x+c','StartPoint',[0.1 0.1 0],'Lower',[0 -Inf,0],'Upper',[Inf Inf 0]); %order-->[a,b,c]
delta_EW         = loadFields(SP,'Weight.EW')-loadFields(SP,'Engine.TotalWeight')-0.2.*loadFields(SP,'Weight.MTOW')-Wfix;%[kg]
fus_parameter    = loadFields(SP,'Fuselage.fusLength').*(loadFields(SP,'Fuselage.fusWidth')+loadFields(SP,'Fuselage.fusHeight'))./2; %[m^2]
indexFus         = ~isnan(fus_parameter);
[delta_EW_fit,delta_EW_R2] = polyfitR2(log10(fus_parameter(indexFus)),log10(delta_EW(indexFus)),1);
fus_parameter_ac = DP.fusLength*(DP.fusWidth+DP.fusHeight)/2;

%Solve nonlinear equation to obtain MTOW before and after EW reduction by MTorres
options = optimoptions('fsolve',...
    'StepTolerance',1e-9,...
    'Display','none');
[MTOW.oldTorenbeek,~,exitflag,~] = fsolve(@(x)getTorenbeekMTOW(x, Wf_Wto_fit, delta_EW_fit, Wfix, DP, Parameters, CST, CF, 'EW_old'),30000,options);
[~,EW.oldTorenbeek, MFW.oldTorenbeek, parametro_ac_old, MFW_MTOW_ac_old] = getTorenbeekMTOW(MTOW.oldTorenbeek, Wf_Wto_fit, delta_EW_fit, Wfix, DP, Parameters, CST, CF, 'EW_old');
if ~isequal(exitflag,1)
    disp('El solver del MTOW por Torenbeek no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear exitflag
end
[MTOW.newTorenbeek,~,exitflag,~] = fsolve(@(x)getTorenbeekMTOW(x, Wf_Wto_fit, delta_EW_fit, Wfix, DP, Parameters, CST, CF, 'EW_new'),MTOW.oldTorenbeek,options);
[~,EW.newTorenbeek, MFW.newTorenbeek, parametro_ac_new, MFW_MTOW_ac_new] = getTorenbeekMTOW(MTOW.newTorenbeek, Wf_Wto_fit, delta_EW_fit, Wfix, DP, Parameters, CST, CF, 'EW_new');
if ~isequal(exitflag,1)
    disp('El solver del MTOW por Torenbeek no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear exitflag options
end


%Fracción de combustible vs parámetro random
if DP.ShowReportFigures
    figure()
    plot(parametro(indexPar),Wf_Wto(indexPar),'*','LineWidth',1,'Color',Parameters.Colors(1,:)); hold on; xl = xlim; yl = ylim;
    plot([0,max(parametro(indexPar))],polyval([Wf_Wto_fit.a, Wf_Wto_fit.b, Wf_Wto_fit.c],[0,max(parametro(indexPar))]),'--','LineWidth',1.25,'Color',Parameters.Colors(2,:))
    plot(parametro_ac_old, MFW_MTOW_ac_old,'+','LineWidth',1.25,'Color',Parameters.Colors(3,:))
    plot(parametro_ac_new, MFW_MTOW_ac_new,'+','LineWidth',1.25,'Color',Parameters.Colors(4,:))
    legend('Similar planes','Polynomial regression','Without weight reduction','With weight reduction','Location','southeast')
    legend('boxoff')
    xlabel('$$\frac{R}{a_{0}}\frac{TSFC_{0}}{\sqrt{\theta}}\left( \frac{1}{M\sqrt{A}}+0.068 P M \frac{l_f\left(b_f+h_f\right)}{2 MTOW}\right)$$','Interpreter','latex')
    ylabel('$$\frac{MFW}{MTOW}$$','Interpreter','latex')
    title('Estimation of fuel weight fraction for similar business jets','Interpreter','latex')
    xlim([xl(1)-0.25,xl(2)+0.1])
    ylim([yl(1)-0.05,yl(2)+0.025])
    %Show equation in the graph
    x = min(parametro(indexPar)) + 0.35*(max(parametro(indexPar))-min(parametro(indexPar)));
    y = polyval([Wf_Wto_fit.a, Wf_Wto_fit.b, Wf_Wto_fit.c],x);
    txt0 = '$$\ \ \leftarrow$$ $$MFW/MTOW=Ax^2+Bx$$';
    txt1 = strcat('$$\ \ \ \ \ \ A=',num2str(Wf_Wto_fit.a),'$$');
    txt2 = strcat('$$\ \ \ \ \ \ B=',num2str(Wf_Wto_fit.b),'$$');
    txt3 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(gof.rsquare),'$$');
    text(x,y,txt0,'Interpreter','latex','FontSize',11)
    text(x,y-0.0075,txt1,'Interpreter','latex','FontSize',11)
    text(x,y-0.0150,txt2,'Interpreter','latex','FontSize',11)
    text(x,y-0.0225,txt3,'Interpreter','latex','FontSize',11)
    clear xl yl x y 
    saveFigure(ME.FiguresFolder,'Torenbeek_FuelFraction')
end

%Incremento de empty weight vs parametro del fuselaje
if DP.ShowReportFigures
    figure()
    loglog(fus_parameter(indexFus),delta_EW(indexFus),'*','LineWidth',1.00,'Color',Parameters.Colors(1,:)); hold on;
    plot([min(fus_parameter(indexFus))-2.5,max(fus_parameter(indexFus))],10.^polyval(delta_EW_fit,log10([min(fus_parameter(indexFus))-2.5,max(fus_parameter(indexFus))])),'--','LineWidth',1.25,'Color',Parameters.Colors(2,:))
    plot(fus_parameter_ac,10^polyval(delta_EW_fit,log10(fus_parameter_ac)),'+','LineWidth',1.25,'Color',Parameters.Colors(4,:))
    legend('Similar planes','Logarithmic regression','Aircraft','Location','southeast')
    legend('boxoff')
    xlabel('$$\frac{l_f\left(b_f+h_f\right)}{2}\ [m^2]$$','Interpreter','latex')
    ylabel('$$\Delta EW\ [kg]$$','Interpreter','latex')
    title('Estimation of the fuselage dependent empty weight','interpreter','latex')
    xlim([min(fus_parameter(indexFus))-5,max(fus_parameter(indexFus))+5]);
    ylim([min(delta_EW(indexFus))-0.5e3,max(delta_EW(indexFus))+1e3]);
    %Show equation in the graph
    x = min(fus_parameter(indexFus)) + 0.15*(max(fus_parameter(indexFus))-min(fus_parameter(indexFus)));
    y = 10^polyval(delta_EW_fit,log10(x));
    txt0 = '$$\ \ \leftarrow$$ $$\log_{10}(\Delta EW)=A+B\log_{10}(\frac{l_f\left(b_f+h_f\right)}{2})$$';
    txt1 = strcat('$$\ \ \ \ \ \ A=',num2str(delta_EW_fit(2)),'$$');
    txt2 = strcat('$$\ \ \ \ \ \ B=',num2str(delta_EW_fit(1)),'$$');
    txt3 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(delta_EW_R2),'$$');
    text(x,y,txt0,'Interpreter','latex','FontSize',11)
    text(x,y-0.35e3,txt1,'Interpreter','latex','FontSize',11)
    text(x,y-0.65e3,txt2,'Interpreter','latex','FontSize',11)
    text(x,y-0.90e3,txt3,'Interpreter','latex','FontSize',11)
    saveFigure(ME.FiguresFolder,'Torenbeek_fuselage_EW')
end

%         txt3 = strcat('$$\ \ \ \ \ \ R^{2}=',num2str(rsquared),'$$');
%         text(x,y,txt,'Interpreter','latex','FontSize',11)

clear Wfix index a0 delta_EW delta_EW_fit delta_EW_R2 fus_parameter fus_parameter_ac parametro parentesis P T T0 Wf_Wto Wf_Wto_fit
clear gof indexFus indexPar MFW_MTOW_ac_old MFW_MTOW_ac_new parametro_ac_old parametro_ac_new reserveRange txt0 txt1 txt2 txt3


% Define more AC weights;
%     AC.Weight.MTOW = mean([MTOW.newRoskam,MTOW.newTorenbeek]);
    AC.Weight.MTOW = MTOW.newTorenbeek + 0.75 * (MTOW.newRoskam - MTOW.newTorenbeek);
    %AC.Weight.EW   = mean([EW.newRoskam,EW.newTorenbeek]);
    %AC.Weight.MFW  = mean([MFW.newRoskam,MFW.newTorenbeek]);
    [~, AC.Weight.MTOW, AC.Weight.EW, AC.Weight.MFW] = getWeights( AC.Weight.MTOW, ME, CST, CF, Parameters, 'EW_new');
    AC.Weight.MRW  = DP.MRW_MTOW * AC.Weight.MTOW;
    AC.Weight.MLW  = DP.MLW_MTOW * AC.Weight.MTOW;
    AC.Weight.TUL  = AC.Weight.MFW + ME.Payload;
    AC.Weight.OEW  = AC.Weight.EW + 0.005*AC.Weight.MTOW + ME.CrewWeight;
    AC.Weight.BOW  = AC.Weight.OEW;

%Define AC fuselage
    AC.Fuselage.fusHeight    = DP.fusHeight;
    AC.Fuselage.fusLength    = DP.fusLength;
    AC.Fuselage.fusWidth     = DP.fusWidth;
    AC.Fuselage.fuselage_AoA = DP.fuselage_AoA;
    AC.Fuselage.Volume       = DP.totalFusVolume;
    AC.Fuselage.cabLength    = DP.cabLength;
    AC.Fuselage.ln           = DP.ln;
    AC.Fuselage.la           = DP.la;
    AC.Fuselage.frontArea    = DP.frontArea;
    AC.Fuselage.Swet         = DP.Swet;
    AC.Fuselage.A_I          = DP.A_I;
    AC.Fuselage.A_II         = DP.A_II;
    AC.Fuselage.beta         = DP.tailConeAngle;



clear MTOW EW MFW

%% USEFUL FUNCTIONS DEFINITION:
function [ F, W_TO_guess, W_E, W_F, W_E_tent] = getWeights( x, ME, CST, CF, Parameters, W_E_Str )
%GETWEIGHTS: Gets the estimation of take-off gross weight (WTO), empty weight (WE) and mission
    % fuel weight (WF). All weights in kg

    % 1. Determine the mission payload weight (W_PL)
    W_PL = ME.Payload; %[kg]


    % 2. Guessing a likely value of W_TO_guess:
    %An initial guess is obtained by comparing the mission specification of the
    %airplane with the mission capabilities of similar airplanes.
    W_TO_guess = x; %[kg]


    % 3. Determination of mission fuel weight:
    %Eq 2.13
    M_ff = 1; %Mission Fuel Fraction
    for i=1:length(Parameters.fuelFraction(:))
        M_ff = M_ff*Parameters.fuelFraction(i).value;
    end

    W_F_res = 0; %NEEDS TO BE ESTABLISHED (FAR?), en mi caso ya están incluidas, así que avisa si las vas a añadir aquí
    W_F = (1 - M_ff)*W_TO_guess + W_F_res; %Eq 2.15 [kg]


    % Step 4. Calculate a tentative value for W_OE from:
    W_OE_tent = W_TO_guess - W_F - W_PL;  %Eq 2.4 [kg]


    % Step 5. Calculate a tentative value for W_E from:
    W_tfo = 0.005*W_TO_guess; % 0.5% of MTOW, taken from example pag.52. Note that W_tfo (trapped fuel-oil) is often neglected in this stage (page 7)
    W_E_tent = W_OE_tent - W_tfo - ME.CrewWeight;   %Eq 2.4.  [kg]

    % 4. Finding the allowable value for W_E
    if strcmp(W_E_Str,string('EW_old'))
        W_E = 10^((log10(W_TO_guess*CST.GravitySI*CF.N2lbf)-Parameters.Table_2_15.a)/Parameters.Table_2_15.b); %In lbf, remember that Roskam correlation is in lbf
        W_E = W_E*CF.lbf2N/CST.GravitySI; %W_E in kg
    elseif strcmp(W_E_Str,string('EW_new'))
        W_E = 10^((log10(W_TO_guess*CST.GravitySI*CF.N2lbf)-(Parameters.Table_2_15.a-Parameters.Table_2_15.b*log10(Parameters.EWnew_EWold)))/Parameters.Table_2_15.b); %In lbf, remember that Roskam correlation is in lbf
        W_E = W_E*CF.lbf2N/CST.GravitySI; %W_E in kg
    else
        error('Incorrect input string for Empty Weight, you must choose between "EW_old" and "EW_new".')
    end

    F = W_E_tent-W_E;

end


function [F, EW, MFW, parametro_ac, MFW_MTOW_ac] = getTorenbeekMTOW(x, Wf_Wto_fit, delta_EW_fit, Wfix, DP, Parameters, CST, CF, W_E_Str )
    MTOW = x; %MTOW in kg
    [T, a, P, ~] = atmosisa(DP.CruiseAltitude); %#ok<ASGLU>
    [T0,a0,~, ~] = atmosisa(0); %#ok<ASGLU>
    Mcruise      = DP.CruiseSpeed/a;
    parentesis_ac = (1./(Mcruise.*sqrt(DP.AspectRatio))) + (0.068.*P.*Mcruise.*DP.fusLength.*(DP.fusWidth+DP.fusHeight)./(2.*CST.GravitySI.*MTOW));
    Reserve_Range = 0.75*CF.hour2sec*DP.CruiseSpeed;
    parametro_ac  = (DP.Range.*1e3+Reserve_Range).*DP.CruiseTSFC.*CF.TSFC2SI.*CST.GravitySI.*parentesis_ac./(a0);%.*sqrt(T./T0));
    MFW_MTOW_ac   = polyval([Wf_Wto_fit.a, Wf_Wto_fit.b, Wf_Wto_fit.c],parametro_ac);
    fus_parameter_ac = DP.fusLength*(DP.fusWidth+DP.fusHeight)/2;
    delta_EW_ac      = 10^polyval(delta_EW_fit,log10(fus_parameter_ac));
    for i=1:length(Parameters.EngineOptions)
        usedEngine(i) = strcmp(Parameters.EngineOptions(i).Model,DP.EngineModel); %#ok<AGROW>
    end
    EnginesWeight =  DP.EngineNumber*Parameters.EngineOptions(find(usedEngine,1)).Weight;
    if strcmp(W_E_Str,string('EW_old'))
        EW = EnginesWeight + Wfix + delta_EW_ac + 0.2*MTOW;
    elseif strcmp(W_E_Str,string('EW_new'))
        EW = DP.EWnew_EWold*(EnginesWeight + Wfix + delta_EW_ac + 0.2*MTOW);
    else
        error('Incorrect input string for Empty Weight, you must choose between "EW_old" and "EW_new".')
    end
    MFW = MFW_MTOW_ac*MTOW;
    F   = EW + DP.Payload + (MFW_MTOW_ac-1)*MTOW;
end
