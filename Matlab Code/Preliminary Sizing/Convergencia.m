switch ME.MissionType
    case 5
        % EXPECTED TO BE RUN AFTER B_loadParameters, IF NOT, COMMENT THIS:
        h1=gcf;
        h2=figure();
        objects=allchild(h1); %#ok<NASGU>
        copyobj(get(h1,'children'),h2);
        
        %Define linspace values of MTOW
        x=linspace(18.6e3,45e3,30);
        for i = 1:length(x)
            [ ~, ~, W_E, ~, W_E_tent,~,~ ] = getWeights(x(i), ME, CST, CF, Parameters );
            y1(i)= W_E_tent; %fuel fraction method
            y2(i)= W_E; %regresion, already painted
        end
        plot(y1,x,'b','LineWidth',1.25,'Color',Parameters.Colors(3,:));
        plot(AC.Weight.EW,AC.Weight.MTOW,'+','LineWidth',2,'Color',Parameters.Colors(4,:));
        txt=['EW:\ ',num2str(round(AC.Weight.EW)),'kg\ \ $$\rightarrow\ \ \ $$'];
        text(AC.Weight.EW,AC.Weight.MTOW,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        txt=['MTOW:\ ',num2str(round(AC.Weight.MTOW)),'kg\ \ \ '];
        text(AC.Weight.EW,AC.Weight.MTOW+1.5e3,txt,'HorizontalAlignment','right','FontSize',10,'Interpreter','Latex')
        h = findobj(gca,'Type','line');
        legend([h(3),h(2)],{'Similar Planes','Fuel-Fraction Method'},'Location','southeast')
        legend('boxoff')
        saveFigure(ME.FiguresFolder,'MTOW definition')
        clear h h1 h2 x y1 y2 W_E W_E_tent i objects txt
        
    case 11
        x1 = linspace(15000,50000, 30);
        for i = 1:length(x1)
            x=x1(i);
            [ F, W_TO_guess, W_E, W_F, W_E_tent,m,M_ff ] = getWeights( x, ME, CST, CF, Parameters );
            y1(i)= W_E_tent;
            y2(i)=W_E;
        end
        
        
        figure()
        plot(log10(y1),log10(x1),'b') ;
        hold all
        plot(log10(y2),log10(x1),'r')
        ylabel('log10(MTOW)')
        xlabel('log10(EW)')
        legend('fraccion combustible','semejantes','Location','southeast')
end


function [ F, W_TO_guess, W_E, W_F, W_E_tent,m, M_ff ] = getWeights( x, ME, CST, CF, Parameters )
%GETWEIGHTS: Gets the estimation of take-off gross weight (WTO), empty weight (WE) and mission
% fuel weight (WF). All weights in kg

%% 1. Determine the mission payload weight (W_PL)
W_PL = ME.Payload; %[kg]


%% 2. Guessing a likely value of W_TO_guess:
%An initial guess is obtained by comparing the mission specification of the
%airplane with the mission capabilities of similar airplanes.
W_TO_guess = x;  %[kg]


%% 3. Determination of mission fuel weight:
 %Eq 2.13
M_ff = 1; %Mission Fuel Fraction
for i=1:length(Parameters.fuelFraction(:))
    M_ff = M_ff*Parameters.fuelFraction(i).value;
end

W_F_res = 0; %NEEDS TO BE ESTABLISHED (FAR?)
W_F = (1 - M_ff)*W_TO_guess + W_F_res; %Eq 2.15 %[kg]


%% Step 4. Calculate a tentative value for W_OE from:
W_OE_tent = W_TO_guess - W_F - W_PL;  %Eq 2.4 %[kg]


%% Step 5. Calculate a tentative value for W_E from:
W_tfo = 0.005*W_TO_guess; % 0.5% of MTOW, taken from example pag.52. Note that W_tfo (trapped fuel-oil) is often neglected in this stage (page 7)
W_E_tent = W_OE_tent - W_tfo - ME.CrewWeight;   %Eq 2.4.  %[kg]

%% 4. Finding the allowable value for W_E
W_E = 10^((log10(W_TO_guess*CST.GravitySI*CF.N2lbf)-Parameters.Table_2_15.a)/Parameters.Table_2_15.b); %In lbf, remember that Roskam correlation is in lbf
W_E = W_E*CF.lbf2N/CST.GravitySI; %W_E in kg

m= W_E/W_TO_guess;
% W_E = 10^((log10(W_TO_guess/CF.lbs_to_kg)-parameters.Table_2_15.a)/parameters.Table_2_15.b); %lb
% W_E = W_E*CF.lbs_to_kg;
% 
% disp(['WTO: ',num2str(W_TO_guess)])
% disp(['WF: ',num2str(W_F)])
% disp(['WE_tent: ',num2str(W_E_tent)])
% disp(['WE: ',num2str(W_E)])
F = W_E_tent-W_E;

end