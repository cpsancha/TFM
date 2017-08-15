%********************************************************************************************************
%*                                                                                                      *
%*                       E - CONFIGURATION DESIGN: PRELIMINARY DESIGN SEQUENCE I                        *
%*                                                                                                      *
%*    In this script...                                                                                 *
%*                                                                                                      *
%********************************************************************************************************


%% OVERALL CONFIGURATION
    % Land based Tandem Wing
    
    
%% FUSELAGE CONFIGURATION --> Chapter 4
    % Conventional
    
    
%% ENGINE TYPE, NUMBER AND DISPOSITION --> Chapter 5
    % Turbojet[x]/ Turbofan[]
    % 2[X] / 3[ ]
    % Pushers
    % In nacelles
    % On the wing [x] / Below the wing[] / Above the wing[ ] / On the fuselage[ ]

        offset=5; %Fila del primer motor
   for i=offset:length(turbopropDataBase)
      enginesPower(i) = ME.Powerplant.Number*turbopropDataBase(i).Power*CF.hp2watts;
      enginesSFC(i)   = turbopropDataBase(i).SFC;
      enginesModel{i} = turbopropDataBase(i).Model;
   end
        
        
index = find(enginesPower<1.1*AC.Engine.TotalPower & enginesPower>0.9*AC.Engine.TotalPower); %Find matching engines in a range
% index = find(enginesPower<1.1*AC.Engine.TotalPower);
indexSFC =  find(enginesSFC<0.5 & ~isnan(enginesSFC));
 
[~,indexmin] = min(enginesSFC(index)); % In that range, select the one whose SFC is minimum

AC.Engine.Model       = enginesModel{index(indexmin)};
AC.Engine.SFC         = enginesSFC(index(indexmin));
AC.Engine.TotalPower  = enginesPower(index(indexmin));
AC.Engine.Power       = turbopropDataBase(index(indexmin)).Power*CF.hp2watts;
AC.Engine.Weight      = turbopropDataBase(index(indexmin)).Weight*CF.lbm2kg;
AC.Engine.Length      = turbopropDataBase(index(indexmin)).Length*CF.in2m;
AC.Engine.Width       = turbopropDataBase(index(indexmin)).Width*CF.in2m;
AC.Engine.Swet        = pi*AC.Engine.Width*AC.Engine.Length; %Aproximacion cutre por un cilindro sin base xd

AC.Weight.Pto_MTOW = AC.Engine.TotalPower/ (AC.Weight.MTOW*CST.GravitySI);
%         [min,indexmin]= min(abs(enginesPower-AC.Engine.TotalPower)); % Closest to power selected
%%        
 figure(); hold on;
 title('Engine selection')
 clear LegendStr
   LegendStr=cell(1);
   grid on
 
%    Engines plot
%  for i=1:length(index) 
%  plot(Wto_S,ones(1,length(Wto_S)).*enginesPower(index(i))./(AC.Weight.MTOW.*CST.GravitySI),'--')
%  LegendStr{i} = enginesModel{index(i)};
%  end
 
%  for i=1:length(indexSFC ) 
%  plot(Wto_S,ones(1,length(Wto_S)).*enginesPower(indexSFC (i))./(AC.Weight.MTOW.*CST.GravitySI),'--')
%  LegendStr{i} = enginesModel{indexSFC(i)};
%  end
 
  for i=1:length(indexSFC) 
 plot(Wto_S,ones(1,length(Wto_S)).*enginesPower(indexSFC(i))./(AC.Weight.MTOW.*CST.GravitySI),'--')
 LegendStr{i} = enginesModel{indexSFC(i)};
  end
 

 %Design point plot
     plot(AC.Wing.WingLoading,AC.Weight.Pto_MTOW,'o')
 %Restrictions   
        plot(Wto_S,P_W.cr,'Color',Parameters.Colors(2,:));         
        plot(Wto_S,P_W.take_off,'Color',Parameters.Colors(3,:));
%         plot(Wto_S,P_W.take_off1,'Color',Parameters.Colors(4,:));
        plot(WingLoading.LandingRoskam.*ones(1,100),linspace(1,100,100),'LineWidth',1.25,'Color',Parameters.Colors(5,:));
        plot(Wto_S,P_W.cl.CS25121tr,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
        plot(Wto_S,P_W.cl.CS25111,'LineWidth',1.25,'Color',Parameters.Colors(7,:));
        plot(Wto_S,P_W.cl.CS25121nd,'LineWidth',1.25,'Color',Parameters.Colors(8,:));         
        plot(Wto_S,P_W.cl.CS25121er,'LineWidth',1.25,'Color',Parameters.Colors(9,:));
        plot(Wto_S,P_W.cl.CS25119,'LineWidth',1.25,'Color',Parameters.Colors(10,:));    
        plot(Wto_S,P_W.cl.CS25121ba,'LineWidth',1.25,'Color',Parameters.Colors(11,:));
        plot(WingLoading.Gust.*ones(1,100),linspace(0,100,100),'LineWidth',1.25,'Color',Parameters.Colors(12,:));
     %Formating
        xlim([250,max(Wto_S)])
        ylim([0,60])
        set(gcf,'Position',[450   200   700   525])
        xlabel('Wing Loading - MTOW/Sw [N/m^2]')
        ylabel('Power/Weight_{TO} [W/N]')
        grafWidth   = 16;
        grafAR      = 0.6;
        set(gcf,'DefaultLineLineWidth',1.5);
        set(gcf,'PaperUnits', 'centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
        set(gca,'FontSize',10,'FontName','Times new Roman','box','on')
        warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
%         [~,objs]=columnlegend(2,LegendStr,'Location','northwest','FontSize',8,'boxoff');
        drawnow; % make sure everything is finished rendering 
%         set(findall(objs, 'type', 'text'), 'fontsize', 8, 'interpreter', 'tex')
        warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');
        clear grafWidth grafAR showRoskamRequirements LegendStr objs
        
        
% Determination of proppeller diameter

Pto_DP2 = (250-140)/(1000)*sqrt(AC.Engine.Power/CF.hp2watts*ME.Cruise.Speed*3.6);
Dp = sqrt(AC.Engine.Power/CF.hp2watts /Pto_DP2 ); 
%  [~, a, ~, rho] = atmosisa(ME.Cruise.Altitude);
% D = a/pi/4*sqrt(0.8^2-ME.Cruise.Mach^2)

AC.Engine.propDiameter = Dp;
AC.Engine.Position= [NaN,(AC.Engine.Power/CF.hp2watts/100*1.65e-2+0.1)+Dp/2+AC.Fuselage.fusWidth/2, NaN];
        
%% WING CONFIGURATION --> Chapter 6 & 7
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
    
    
%% EMPENNAGE CONFIGURATION --> Chapter 8
    % HORIZONTAL TAIL ??
        % Fuselage mounted[ ] / Vertical tail mounted[ ] / V-tail[ ]
    % VERTICAL TAIL ??
        % Fuselage mounted[ ] / V-tail[ ]
    % STRAKERS ??

    
%% LANDING GEAR CONFIGURATION --> Chapter 9
    % Retractable
    % Tricicle[X] / Tandem[ ] / Outrigger[ ]
    % Fuselage mounted[ ] / Wing mounted[ ] / Nacelle mounted[ ]
    % Number of main gear struts
    % Number of tyres per strut
    % Retractions kinematics and available volumen to receive the gear

    
    
%% JET FUEL --> VOLUME
%The used fuel is known as Jet A-1 and it's density at 15 °C (59 °F) is 0.804 kg/l (6.71 lb/US gal)











%% OTHER IDEAS
% Your idea of moving the wings apart is sound. However, in order to make the plane statically stable, the rear wing needs to produce less lift per 
% area relative to the forward wing, making it less efficient. If you now run an optimizer which varies wing area and minimizes overall drag, you 
% will invariably end up with a conventional design in which the rear wing has 15% - 20% of the area of the main wing.

% What is the optimum distance between the wings? This depends entirely on your preferences. You need to know that stability increases linearly with 
% wing-tail distance while pitch damping increases with the square of it. If you want an agile aircraft, keep both close together. If you want to 
% reduce the tail surface and don't intend to fly aerobatics, move both apart. The key figure for stability is the tail volume, which is the product 
% of tail surface area and the distance between wing and tail. Look at airplanes you think are worth emulating and try to match the tail volume of your 
% design to theirs. But don't use a magic number.
















