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
    % Below the wing[X] / Above the wing[ ] / On the fuselage[ ]
    % Options:
        % [ ] 2*Rolls-Royce AE 3007A1E (39.7kN) -->  ~80kN ~kg ~lb/(lbf�h)--> Embraer Legacy family
        % [X] 2*Snecma Silvercrest (35 - 53 kN) --> ~100kN ~2*1040kg ~0.628lb/(lbf�h)--> Dassault Falcon 5X (51kN) - Cessna Citation Hemisphere (53kN)
        % [ ] 2*Rolls-Royce BR710A2-20 (65.6kN) --> ~130kN ~2*kg ~0.63lb/(lbf�h)--> Bombardier Global 5000
        % [ ] 2*P&W Canada PW800 (44 to 89 kN)  -->  ~90kN ~kg ~lb/(lbf�h)--> Gulfstream G500/G600 (67.36/69.75 kN)
        % [ ] 2*GE Passport (44 to 89 kN)       -->  ~90kN ~kg ~lb/(lbf�h)--> Bombardier Global 7000/8000 (73.4 kN)

        offset=5 %Fila del primer motor
   for i=offset:length(turbopropDataBase)
      enginesPower(i) = ME.Powerplant.Number*turbopropDataBase(i).Power*CF.hp2watts;
      enginesSFC(i)   = turbopropDataBase(i).SFC;
      enginesModel{i} = turbopropDataBase(i).Model;
   end
        
        
index = find(enginesPower<1.03*AC.Engine.TotalPower & enginesPower>0.97*AC.Engine.TotalPower); %Find matching engines in a range

[~,indexmin] = min(enginesSFC(index)); % In that range, select the one whose SFC is minimum

AC.Engine.Model       = enginesModel{index(indexmin)};
AC.Engine.SFC         = enginesSFC(index(indexmin));
AC.Engine.TotalPower  = enginesPower(index(indexmin));
AC.Engine.Power       = turbopropDataBase(index(indexmin)).Power*CF.hp2watts;
%         [min,indexmin]= min(abs(enginesPower-AC.Engine.TotalPower)); % Closest to power selected
        
 figure(4); hold on;
 title('Engine selection')
 clear LegendStr
   LegendStr=cell(0);
 
 for i=1:length(index) 
 plot(Wto_S,ones(1,length(Wto_S)).*enginesPower(index(i))./(AC.Weight.MTOW.*CST.GravitySI),'--')
 LegendStr{i} = enginesModel{index(i)};
 end
  
     plot(AC.Wing.WingLoading,AC.Weight.Pto_MTOW,'o')
     
        plot(Wto_S,P_W.cr,'Color',Parameters.Colors(2,:));         
        plot(Wto_S,P_W.take_off,'Color',Parameters.Colors(3,:));
        plot(Wto_S,P_W.take_off1,'Color',Parameters.Colors(4,:));
        plot(WingLoading.LandingRoskam.*ones(1,100),linspace(1,100,100),'LineWidth',1.25,'Color',Parameters.Colors(5,:));
        plot(Wto_S,P_W.cl.CS25121tr,'LineWidth',1.25,'Color',Parameters.Colors(6,:));
        plot(Wto_S,P_W.cl.CS25111,'LineWidth',1.25,'Color',Parameters.Colors(7,:));
        plot(Wto_S,P_W.cl.CS25121nd,'LineWidth',1.25,'Color',Parameters.Colors(8,:));         
        plot(Wto_S,P_W.cl.CS25121er,'LineWidth',1.25,'Color',Parameters.Colors(9,:));
        plot(Wto_S,P_W.cl.CS25119,'LineWidth',1.25,'Color',Parameters.Colors(10,:));    
        plot(Wto_S,P_W.cl.CS25121ba,'LineWidth',1.25,'Color',Parameters.Colors(11,:));
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
        [~,objs]=columnlegend(2,LegendStr,'Location','northwest','FontSize',8,'boxoff');
        drawnow; % make sure everything is finished rendering 
        set(findall(objs, 'type', 'text'), 'fontsize', 8, 'interpreter', 'tex')
        warning('on', 'MATLAB:handle_graphics:exceptions:SceneNode');
        clear grafWidth grafAR showRoskamRequirements LegendStr objs
        
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
%The used fuel is known as Jet A-1 and it's density at 15 �C (59 �F) is 0.804 kg/l (6.71 lb/US gal)











%% OTHER IDEAS
% Your idea of moving the wings apart is sound. However, in order to make the plane statically stable, the rear wing needs to produce less lift per 
% area relative to the forward wing, making it less efficient. If you now run an optimizer which varies wing area and minimizes overall drag, you 
% will invariably end up with a conventional design in which the rear wing has 15% - 20% of the area of the main wing.

% What is the optimum distance between the wings? This depends entirely on your preferences. You need to know that stability increases linearly with 
% wing-tail distance while pitch damping increases with the square of it. If you want an agile aircraft, keep both close together. If you want to 
% reduce the tail surface and don't intend to fly aerobatics, move both apart. The key figure for stability is the tail volume, which is the product 
% of tail surface area and the distance between wing and tail. Look at airplanes you think are worth emulating and try to match the tail volume of your 
% design to theirs. But don't use a magic number.















