

%Calcular alcance máximo en configuración de diseño
Winit = AC.Weight.MTOW;    %Cruise initial weight
Wend  = AC.Weight.OEW + ME.Payload; %Cruise final weight
[R, ~] = getRangeEndurance(Winit, Wend ,AC, DP, CST, CF, Polar);
disp(' ')
disp(['ALCANCE: ',num2str(R/1e3,'%10.2f'),' km.'])

%Payload-Range Diagram
if DP.ShowReportFigures
    getPayloadRangeDiagram(AC, CF, CST, DP, ME, Polar, ME.FiguresFolder, 'PayloadRangeDiagram');
%     getWeightRangeDiagram(AC, CF, CST, DP, ME, Polar, ME.FiguresFolder, 'WeightRangeDiagram')

end

clear Winit Wend



%% USEFUL FUNCTIONS
function [] = getPayloadRangeDiagram(AC, CF, CST, DP, ME, Polar, folderName, figureName)

%Get Points
Reservas = 360; %[km]
Ra = getRangeEndurance(AC.Weight.OEW + AC.Weight.TUL, AC.Weight.OEW + DP.maxPayload, AC, DP, CST, CF, Polar)/1e3 - Reservas;
Rb = getRangeEndurance(AC.Weight.OEW + AC.Weight.TUL, AC.Weight.OEW + DP.Payload, AC, DP, CST, CF, Polar)/1e3 - Reservas;
Rc = getRangeEndurance(AC.Weight.OEW + AC.Weight.MFW, AC.Weight.EW + ME.CrewWeight ,AC, DP, CST, CF, Polar)/1e3 - Reservas;
Rd = getRangeEndurance(AC.Weight.OEW + AC.Weight.MFW, AC.Weight.EW + ME.CrewWeight ,AC, DP, CST, CF, Polar)/1e3;


%Create figure
figure();
hold on
LegendStr = cell(0);
colors = get(gca,'colororder');

%Constant values
%     Rmax     = 10500; %[km]
%     
%     %MTOW
%     plot([0,Rmax],[AC.Weight.MTOW,AC.Weight.MTOW],'--','LineWidth',1.25,'Color',Parameters.Colors(1,:))
%     LegendStr{end+1} = 'MTOW';
%     
%     %MLW
%     plot([0,Rmax],[AC.Weight.MLW,AC.Weight.MLW],'--','LineWidth',1.25,'Color',Parameters.Colors(2,:))
%     LegendStr{end+1} = 'MLW';
%     
%     %MZFW
%     plot([0,Rmax],[AC.Weight.MZFW,AC.Weight.MZFW],'--','LineWidth',1.25,'Color',Parameters.Colors(3,:))
%     LegendStr{end+1} = 'MZFW';
%     
%     %OEW
%     plot([0,Rmax],[AC.Weight.OEW,AC.Weight.OEW],'--','LineWidth',1.25,'Color',Parameters.Colors(4,:))
%     LegendStr{end+1} = 'OEW';
    
    %Point A
    plot(Ra,DP.maxPayload,'*','LineWidth',1.00,'Color',colors(2,:))
    LegendStr{end+1} = 'Max. Payload';
    
    %Point B
    plot(Rb,DP.Payload,'*','LineWidth',1.00,'Color',colors(3,:))
    LegendStr{end+1} = 'Design Payload';
    
    %Point C
    plot(Rc,0,'*','LineWidth',1.00,'Color',colors(4,:))
    LegendStr{end+1} = 'Max. Range w/ Reserves';
    
    %Point D
    plot(Rd,0,'*','LineWidth',1.00,'Color',colors(5,:))
    LegendStr{end+1} = 'Max. Range w/o Reserves';
    
    %Lines
    plot([0, Ra],[DP.maxPayload, DP.maxPayload],'LineWidth',1.25,'Color',colors(1,:))
    plot([Ra,Rb],[DP.maxPayload, DP.Payload],'LineWidth',1.25,'Color',colors(1,:))
    plot([Rb,Rc],[DP.Payload, 0],'LineWidth',1.25,'Color',colors(1,:))
    plot([Rc,Rd],[0, 0],'LineWidth',1.25,'Color',colors(1,:))
    
    xlim([0,Rd+250])
    ylim([0,2e3])
    legend(LegendStr,'location','southwest','interpreter','latex')
    legend('boxoff')
    xlabel('Alcance (R) [km]','interpreter','latex')
    ylabel('Carga de pago (PL) [kg]','interpreter','latex')
    title('Diagrama de carga de pago - alcance','interpreter','latex')
    
    %Formating
     grafWidth   = 16;
     grafAR      = 0.45;
%      grid on
     set(gcf,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     %Saving
     folderName = strrep(folderName,'.','');
     folderPath = [pwd filesep folderName];
     figureName = strrep(figureName,'.','');
     format_Grafico = [folderPath filesep figureName];
     exists=exist(folderPath,'dir');
     if isequal(exists,0)
         mkdir(folderPath);
         saveas(gcf,format_Grafico,'epsc'); 
     elseif isequal(exists,7)
         saveas(gcf,format_Grafico,'epsc');
     else
         disp('Imposible to save the Figure, check what can be failing.')
     end
     
     
     %Detail
     xlim([9e3,Rd+150])
     ylim([0,2e3])
     %Formating
     grafWidth   = 16;
     grafAR      = 0.45;
%      grid on
     set(gcf,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     %Saving
     folderName = strrep(folderName,'.','');
     folderPath = [pwd filesep folderName];
     figureName = strrep(figureName,'.','');
     format_Grafico = [folderPath filesep figureName '_Detail'];
     exists=exist(folderPath,'dir');
     if isequal(exists,0)
         mkdir(folderPath);
         saveas(gcf,format_Grafico,'epsc'); 
     elseif isequal(exists,7)
         saveas(gcf,format_Grafico,'epsc');
     else
         disp('Imposible to save the Figure, check what can be failing.')
     end

end


function [] = getWeightRangeDiagram(AC, CF, CST, DP, ME, Polar, folderName, figureName)

%Get Points
Reservas = 360; %[km]
Ra = getRangeEndurance(AC.Weight.OEW + AC.Weight.TUL, AC.Weight.OEW + DP.maxPayload, AC, DP, CST, CF, Polar)/1e3 - Reservas;
Rb = getRangeEndurance(AC.Weight.OEW + AC.Weight.TUL, AC.Weight.OEW + DP.Payload, AC, DP, CST, CF, Polar)/1e3 - Reservas;
Rc = getRangeEndurance(AC.Weight.OEW + AC.Weight.MFW, AC.Weight.EW + ME.CrewWeight ,AC, DP, CST, CF, Polar)/1e3 - Reservas;
Rd = getRangeEndurance(AC.Weight.OEW + AC.Weight.MFW, AC.Weight.EW + ME.CrewWeight ,AC, DP, CST, CF, Polar)/1e3;


%Create figure
figure();
hold on
LegendStr = cell(0);
colors = get(gca,'colororder');

%Constant values
    Rmax     = 10750; %[km]
    
    %MTOW
    plot([0,Rmax],[AC.Weight.MTOW,AC.Weight.MTOW],'--','LineWidth',1.25,'Color',colors(1,:))
    LegendStr{end+1} = 'MTOW';
    
    %MLW
    plot([0,Rmax],[AC.Weight.MLW,AC.Weight.MLW],'--','LineWidth',1.25,'Color',colors(2,:))
    LegendStr{end+1} = 'MLW';
    
    %MZFW
    plot([0,Rmax],[AC.Weight.MZFW,AC.Weight.MZFW],'--','LineWidth',1.25,'Color',colors(3,:))
    LegendStr{end+1} = 'MZFW';
    
    %OEW
    plot([0,Rmax],[AC.Weight.OEW,AC.Weight.OEW],'--','LineWidth',1.25,'Color',colors(4,:))
    LegendStr{end+1} = 'OEW';
    
    %Point A
    plot(Ra,AC.Weight.OEW + DP.maxPayload,'*','LineWidth',1.00,'Color',colors(2,:))
    LegendStr{end+1} = 'Max. Payload';
    
    %Point B
    plot(Rb,AC.Weight.OEW + DP.Payload,'*','LineWidth',1.00,'Color',colors(3,:))
    LegendStr{end+1} = 'Design Payload';
    
    %Point C
    plot(Rc,AC.Weight.OEW,'*','LineWidth',1.00,'Color',colors(4,:))
    LegendStr{end+1} = 'Max. Range w/ Reserves';
    
    %Point D
    plot(Rd,AC.Weight.OEW,'*','LineWidth',1.00,'Color',colors(5,:))
    LegendStr{end+1} = 'Max. Range w/o Reserves';
    
    %Lines
    plot([0, Ra],[AC.Weight.OEW + DP.maxPayload, AC.Weight.OEW + DP.maxPayload],'LineWidth',1.25,'Color',colors(1,:))
    plot([Ra,Rb],[AC.Weight.OEW + DP.maxPayload, AC.Weight.OEW + DP.Payload],'LineWidth',1.25,'Color',colors(1,:))
    plot([Rb,Rc],[AC.Weight.OEW + DP.Payload, AC.Weight.OEW + 0],'LineWidth',1.25,'Color',colors(1,:))
    plot([Rc,Rd],[AC.Weight.OEW + 0, AC.Weight.OEW + 0],'LineWidth',1.25,'Color',colors(1,:))
    
    xlim([0,Rmax])
%     ylim([0,2e3])
    legend(LegendStr,'location','southwest','interpreter','latex')
    legend('boxoff')
    xlabel('Alcance (R) [km]','interpreter','latex')
    ylabel('Carga de pago (PL) [kg]','interpreter','latex')
    title('Diagrama de carga de pago - alcance','interpreter','latex')
    
    %Formating
     grafWidth   = 16;
     grafAR      = 0.45;
%      grid on
     set(gcf,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     %Saving
     folderName = strrep(folderName,'.','');
     folderPath = [pwd filesep folderName];
     figureName = strrep(figureName,'.','');
     format_Grafico = [folderPath filesep figureName];
     exists=exist(folderPath,'dir');
     if isequal(exists,0)
         mkdir(folderPath);
         saveas(gcf,format_Grafico,'epsc'); 
     elseif isequal(exists,7)
         saveas(gcf,format_Grafico,'epsc');
     else
         disp('Imposible to save the Figure, check what can be failing.')
     end
     
     
     %Detail
%      xlim([9e3,Rd+150])
%      ylim([0,2e3])
     %Formating
     grafWidth   = 16;
     grafAR      = 0.45;
%      grid on
     set(gcf,'DefaultLineLineWidth',1.5);
     set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
     set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
     %Saving
     folderName = strrep(folderName,'.','');
     folderPath = [pwd filesep folderName];
     figureName = strrep(figureName,'.','');
     format_Grafico = [folderPath filesep figureName '_Detail'];
     exists=exist(folderPath,'dir');
     if isequal(exists,0)
         mkdir(folderPath);
         saveas(gcf,format_Grafico,'epsc'); 
     elseif isequal(exists,7)
         saveas(gcf,format_Grafico,'epsc');
     else
         disp('Imposible to save the Figure, check what can be failing.')
     end

end


function fig = getGustDiagram()

    % Envolvente
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
    
    
    %Create figure
    fig = figure();
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

end


function [R, E] = getRangeEndurance(Winit, Wend ,AC, DP, CST, CF, Polar)

    %Atmospheric conditions
    [~, ~, ~, rho] = atmosisa(DP.CruiseAltitude);
    
    % Cruise range
    Weight_cruise = linspace(Winit*CST.GravitySI,Wend*CST.GravitySI,250);
    CL_cruise = Weight_cruise./(0.5*rho*DP.CruiseSpeed^2*AC.Wing.Sw);
    CD_cruise = interp1(Polar.CL,Polar.CD,CL_cruise);

    % Integration (Breguet? pls.)
    R = -trapz(Weight_cruise,DP.CruiseSpeed.*(CL_cruise./CD_cruise)./(Weight_cruise.*CST.GravitySI*DP.CruiseTSFC*CF.TSFC2SI));
    E = R/DP.CruiseSpeed;

end