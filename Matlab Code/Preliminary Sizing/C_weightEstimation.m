% Estimating take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

% Get a first guess of W_TO from similar airplanes:
switch ME.MissionType
    case 5
        W_TO_guess = 30000; %[kg]
    case 11
        W_TO_guess = 40000; %[kg]
end


% Solve the iterative process without weight reduction due to application of composites on fuselage:
options = optimoptions('fsolve',...
                       'StepTolerance',1e-9,...
                       'Display','none');
[W_TO_old,~,exitflag,~] = fsolve(@(x)getWeights(x,ME,CST,CF,Parameters,'EW_old'),W_TO_guess,options);
 W_E_old = (10^((log10(W_TO_old*CST.GravitySI*CF.N2lbf)-Parameters.Table_2_15.a)/Parameters.Table_2_15.b))*CF.lbf2N/CST.GravitySI;
if ~isequal(exitflag,1)
    disp('El solver no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear W_TO_guess exitflag
end


% Solve the iterative process with weight reduction due to application of composites on fuselage:
[W_TO_new,~,exitflag,~] = fsolve(@(x)getWeights(x,ME,CST,CF,Parameters,'EW_new'),W_TO_old,options);
if ~isequal(exitflag,1)
    disp('El solver no ha logrado converger correctamente. Se debería revisar el resultado.')
    pause
else
    clear options exitflag
end


% Get the other weights:
[~, AC.Weight.MTOW, AC.Weight.EW, AC.Weight.MFW] = getWeights( W_TO_new, ME, CST, CF, Parameters, 'EW_new');


% Define more AC weights;
switch ME.MissionType
    case 5
        AC.Weight.MRW = DP.MRW_MTOW * AC.Weight.MTOW;
        AC.Weight.MLW = DP.MLW_MTOW * AC.Weight.MTOW;
        AC.Weight.TUL = AC.Weight.MFW + ME.Payload;
        AC.Weight.OEW = AC.Weight.EW + 0.005*AC.Weight.MTOW + ME.CrewWeight;
        AC.Weight.BOW = AC.Weight.OEW;
    case 11
end


%Create figure showing the convergence
switch ME.MissionType
    case 5
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
            y1(i) = W_E_tent; %Fuel fraction method
            y2(i) = W_E;      %New EW regresion, with the weight reduction corrections

        end
        plot(y2,x,'LineWidth',1.25,'Color',Parameters.Colors(3,:)); %<-- EW regresion with weight reduction corrections
        plot(y1,x,'LineWidth',1.25,'Color',Parameters.Colors(6,:)); %<-- EW calculated with the fuel fraction method
        %EW_old intersection point
        plot(W_E_old,W_TO_old,'+','LineWidth',2,'Color',Parameters.Colors(5,:));
        txt=['$$EW_{old}:\ $$',num2str(round(W_E_old)),'kg\ \ $$\rightarrow\ \ \ $$'];
        text(W_E_old,W_TO_old,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        txt=['$$MTOW_{old}:\ $$',num2str(round(W_TO_old)),'kg\ \ \ '];
        text(W_E_old,W_TO_old+1.5e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        %EW_new intersection point
        plot(AC.Weight.EW,AC.Weight.MTOW,'+','LineWidth',2,'Color',Parameters.Colors(4,:));
        txt=['EW:\ ',num2str(round(AC.Weight.EW)),'kg\ \ $$\rightarrow\ \ \ $$'];
        text(AC.Weight.EW,AC.Weight.MTOW,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        txt=['MTOW:\ ',num2str(round(AC.Weight.MTOW)),'kg\ \ \ '];
        text(AC.Weight.EW,AC.Weight.MTOW+1.5e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        %Legend
        h = findobj(gca,'Type','line');
        legend([h(5),h(4),h(3)],{'Similar Planes','Similar Planes (Weight correction)','Fuel-Fraction Method'},'Location','southeast')
        legend('boxoff')
        saveFigure(ME.FiguresFolder,'MTOW definition')
        clear h h1 h2 x y1 y2 W_E W_E_tent i objects txt
        end
        clear W_TO_old W_TO_new W_E_old
    case 11
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
            y1(i) = W_E_tent; %Fuel fraction method
            y2(i) = W_E;      %New EW regresion, with the weight reduction corrections

        end
        plot(y2,x,'LineWidth',1.25,'Color',Parameters.Colors(3,:)); %<-- EW regresion with weight reduction corrections
        plot(y1,x,'LineWidth',1.25,'Color',Parameters.Colors(6,:)); %<-- EW calculated with the fuel fraction method
        %EW_old intersection point
        plot(W_E_old,W_TO_old,'+','LineWidth',2,'Color',Parameters.Colors(5,:));
        txt=['$$EW_{old}:\ $$',num2str(round(W_E_old)),'kg\ \ $$\rightarrow\ \ \ $$'];
        text(W_E_old,W_TO_old,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        txt=['$$MTOW_{old}:\ $$',num2str(round(W_TO_old)),'kg\ \ \ '];
        text(W_E_old,W_TO_old+1.5e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        %EW_new intersection point
        plot(AC.Weight.EW,AC.Weight.MTOW,'+','LineWidth',2,'Color',Parameters.Colors(4,:));
        txt=['EW:\ ',num2str(round(AC.Weight.EW)),'kg\ \ $$\rightarrow\ \ \ $$'];
        text(AC.Weight.EW,AC.Weight.MTOW,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        txt=['MTOW:\ ',num2str(round(AC.Weight.MTOW)),'kg\ \ \ '];
        text(AC.Weight.EW,AC.Weight.MTOW+1.5e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        %Legend
        h = findobj(gca,'Type','line');
        legend([h(5),h(4),h(3)],{'Similar Planes','Similar Planes (Weight correction)','Fuel-Fraction Method'},'Location','southeast')
        legend('boxoff')
        saveFigure(ME.FiguresFolder,'MTOW definition')
        clear h h1 h2 x y1 y2 W_E W_E_tent i objects txt W_TO_old W_TO_new W_E_old
        end      
%         x1 = linspace(15000,50000, 30);
%         for i = 1:length(x1)
%             x=x1(i);
%             [ F, W_TO_guess, W_E, W_F, W_E_tent,m,M_ff ] = getWeights( x, ME, CST, CF, Parameters );
%             y1(i)= W_E_tent;
%             y2(i)=W_E;
%         end
%         
%         
%         figure()
%         plot(log10(y1),log10(x1),'b') ;
%         hold all
%         plot(log10(y2),log10(x1),'r')
%         ylabel('log10(MTOW)')
%         xlabel('log10(EW)')
%         legend('fraccion combustible','semejantes','Location','southeast')
end


%% USEFUL FUNCTIONS DEFINITION:
function [ F, W_TO_guess, W_E, W_F, W_E_tent] = getWeights( x, ME, CST, CF, Parameters, W_E_Str )
%GETWEIGHTS: Gets the estimation of take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

%% 1. Determine the mission payload weight (W_PL)
W_PL = ME.Payload; %[kg]


%% 2. Guessing a likely value of W_TO_guess:
%An initial guess is obtained by comparing the mission specification of the
%airplane with the mission capabilities of similar airplanes.
W_TO_guess = x; %[kg]


%% 3. Determination of mission fuel weight:
 %Eq 2.13
M_ff = 1; %Mission Fuel Fraction
for i=1:length(Parameters.fuelFraction(:))
    M_ff = M_ff*Parameters.fuelFraction(i).value;
end

W_F_res = 0; %NEEDS TO BE ESTABLISHED (FAR?), en mi caso ya están incluidas, así que avisa si las vas a añadir aquí
W_F = (1 - M_ff)*W_TO_guess + W_F_res; %Eq 2.15 [kg]


%% Step 4. Calculate a tentative value for W_OE from:
W_OE_tent = W_TO_guess - W_F - W_PL;  %Eq 2.4 [kg]


%% Step 5. Calculate a tentative value for W_E from:
W_tfo = 0.005*W_TO_guess; % 0.5% of MTOW, taken from example pag.52. Note that W_tfo (trapped fuel-oil) is often neglected in this stage (page 7)
W_E_tent = W_OE_tent - W_tfo - ME.CrewWeight;   %Eq 2.4.  [kg]

%% 4. Finding the allowable value for W_E
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