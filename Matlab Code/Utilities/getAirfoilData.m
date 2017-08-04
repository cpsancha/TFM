function [Cl_alpha,alpha_zeroLift,Cm_ac,Cl,Cd,Cm] = getAirfoilData(desiredMach,desiredReynolds,Colors,FiguresFolder,plotFlag,varargin)

%CHECK IF EXIST desiredAlpha
if isempty(varargin)
    Cl = [];
    Cd = [];
    Cm = [];
    flag = false;
elseif length(varargin)==1
    desiredAlpha = varargin{1};
    flag = true;
else
    error('Unexpected number of inputs')
end


%LOAD AND JOIN STORED AIRFOIL DATA
    load('.\Temporary Stuff\Temp.mat')
    %alpha_Cn
        varnames = who('alpha_Cn*');
        values = cellfun(@eval, varnames, 'UniformOutput', false);
        Mach  = [];
        Re    = [];
        Cn    = [];
        Other = [];
        
        for i=1:length(varnames)
            Mach  = cat(1,Mach ,values{i}(:,3));
            Re    = cat(1,Re,   values{i}(:,4));
            Cn    = cat(1,Cn,   values{i}(:,2));
            Other = cat(1,Other,values{i}(:,1));
        end
        if plotFlag
            figure()
            hold on
            LegendStr=cell(0);
            for i=1:length(varnames)
                plot(values{i}(:,1),values{i}(:,2),'LineWidth',1.25,'Color',Colors(i,:));
                LegendStr{end+1} = ['$M_\infty=',num2str(values{i}(1,3)),'\ \ Re=',num2str(values{i}(1,4)/1e6),'e6$']; %#ok<AGROW>
            end
            xlabel('$\alpha\ [^o]$','interpreter','latex')
            ylabel('$C_l\ [-]$','interpreter','latex')
            xlim([-5,10])
            ylim([-0.25,1.2])
            legend(LegendStr,'Location','southeast','interpreter','latex')
            legend('boxoff')
%             set(gca,'FontSize',12,'FontName','Times new Roman','box','on')
%             warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
%             [~,objs]=columnlegend(2,LegendStr,'Location','southeast','FontSize',12,'padding', 0.5,'boxoff');
%             drawnow; % make sure everything is finished rendering 
%             set(findall(objs, 'type', 'text'),'FontName','Times new Roman','FontSize', 10, 'interpreter', 'latex')
            saveFigure(FiguresFolder,'SC(3)-0712(B) - Cl_alpha')
        end
        Cn_alpha  = cat(2,Mach,Re,Cn,Other);
    %Cd_Cn
        varnames = who('Cd_Cn*');
        values = cellfun(@eval, varnames, 'UniformOutput', false);
        Mach  = [];
        Re    = [];
        Cn    = [];
        Other = [];
        for i=1:length(varnames)
            Mach  = cat(1,Mach ,values{i}(:,3));
            Re    = cat(1,Re,   values{i}(:,4));
            Cn    = cat(1,Cn,   values{i}(:,2));
            Other = cat(1,Other,values{i}(:,1));
        end
        if plotFlag
            figure()
            hold on
            LegendStr=cell(0);
            for i=1:length(varnames)
                plot(values{i}(:,1),values{i}(:,2),'LineWidth',1.25,'Color',Colors(i,:));
                LegendStr{end+1} = ['$M_\infty=',num2str(values{i}(1,3)),'\ \ Re=',num2str(values{i}(1,4)/1e6),'e6$']; %#ok<AGROW>
            end
            xlabel('$C_d\ [-]$','interpreter','latex')
            ylabel('$C_l\ [-]$','interpreter','latex')
%             xlim([-5,10])
%             ylim([-0.25,1.2])
            legend(LegendStr,'Location','southeast','interpreter','latex')
            legend('boxoff')
%             set(gca,'FontSize',12,'FontName','Times new Roman','box','on')
%             warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
%             [~,objs]=columnlegend(2,LegendStr,'Location','southeast','FontSize',12,'padding', 0.5,'boxoff');
%             drawnow; % make sure everything is finished rendering 
%             set(findall(objs, 'type', 'text'),'FontName','Times new Roman','FontSize', 10, 'interpreter', 'latex')
            saveFigure(FiguresFolder,'SC(3)-0712(B) - Cl_Cd')
        end
        Cn_Cd  = cat(2,Mach,Re,Cn,Other);
    %Cm_Cn
        varnames = who('Cm_Cn*');
        values = cellfun(@eval, varnames, 'UniformOutput', false);
        Mach  = [];
        Re    = [];
        Cn    = [];
        Other = [];
        for i=1:length(varnames)
            Mach  = cat(1,Mach ,values{i}(:,3));
            Re    = cat(1,Re,   values{i}(:,4));
            Cn    = cat(1,Cn,   values{i}(:,2));
            Other = cat(1,Other,values{i}(:,1));
        end
        if plotFlag
            figure()
            hold on
            LegendStr=cell(0);
            for i=1:length(varnames)
                plot(values{i}(:,1),values{i}(:,2),'LineWidth',1.25,'Color',Colors(i,:));
                LegendStr{end+1} = ['$M_\infty=',num2str(values{i}(1,3)),'\ \ Re=',num2str(values{i}(1,4)/1e6),'e6$']; %#ok<AGROW>
            end
            xlabel('$C_m\ [-]$','interpreter','latex')
            ylabel('$C_l\ [-]$','interpreter','latex')
            xlim([-0.2,0.1])
%             ylim([-0.25,1.2])
            legend(LegendStr,'Location','southeast','interpreter','latex')
            legend('boxoff')
%             set(gca,'FontSize',12,'FontName','Times new Roman','box','on')
%             warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
%             [~,objs]=columnlegend(2,LegendStr,'Location','southeast','FontSize',12,'padding', 0.5,'boxoff');
%             drawnow; % make sure everything is finished rendering 
%             set(findall(objs, 'type', 'text'),'FontName','Times new Roman','FontSize', 10, 'interpreter', 'latex')
            saveFigure(FiguresFolder,'SC(3)-0712(B) - Cl_Cm')
        end
        Cn_Cm  = cat(2,Mach,Re,Cn,Other);


%LIFT CURVE SLOPE & ZERO LIFT ANGLE
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_alpha);
    lowMachLowReFit = polyfit(deg2rad(Cn_alpha(lowMachLowReIndex,4)),Cn_alpha(lowMachLowReIndex,3),1);
    lowMachUpReFit  = polyfit(deg2rad(Cn_alpha( lowMachUpReIndex,4)),Cn_alpha( lowMachUpReIndex,3),1);
    upMachLowReFit  = polyfit(deg2rad(Cn_alpha( upMachLowReIndex,4)),Cn_alpha( upMachLowReIndex,3),1);
    upMachUpReFit   = polyfit(deg2rad(Cn_alpha(  upMachUpReIndex,4)),Cn_alpha(  upMachUpReIndex,3),1);
    %Cl_alpha    
        lowMach_Cl_alpha = mean([lowMachLowReFit(1),lowMachUpReFit(1)]);
        upMach_Cl_alpha  = mean([ upMachLowReFit(1), upMachUpReFit(1)]);
        Cl_alpha = mean([lowMach_Cl_alpha,upMach_Cl_alpha]); 
    %alpha_zeroLift    
        lowMach_alpha_zeroLift = mean([rad2deg(-lowMachLowReFit(2)/lowMachLowReFit(1)),rad2deg(-lowMachUpReFit(2)/lowMachUpReFit(1))]);
        upMach_alpha_zeroLift  = mean([rad2deg(- upMachLowReFit(2)/ upMachLowReFit(1)),rad2deg(- upMachUpReFit(2)/ upMachUpReFit(1))]);
        alpha_zeroLift = mean([lowMach_alpha_zeroLift, upMach_alpha_zeroLift]); 


%AERODYNAMIC CENTER MOMENTUM COEFFICIENT
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_Cm);
    lowMachLowReCm_ac = mean(Cn_Cm(lowMachLowReIndex,4));
    lowMachUpReCm_ac  = mean(Cn_Cm( lowMachUpReIndex,4));
    upMachLowReCm_ac  = mean(Cn_Cm( upMachLowReIndex,4));
    upMachUpReCm_ac   = mean(Cn_Cm(  upMachUpReIndex,4));
    lowMach_Cm_ac = mean([lowMachLowReCm_ac,lowMachUpReCm_ac]);
    upMach_Cm_ac  = mean([ upMachLowReCm_ac, upMachUpReCm_ac]);
    Cm_ac = mean([lowMach_Cm_ac,upMach_Cm_ac]); 
        

if flag    
%NORMAL FORCE COEFFICIENT
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_alpha);
    lowMachLowReCl = interp1(Cn_alpha(lowMachLowReIndex,4),Cn_alpha(lowMachLowReIndex,3),desiredAlpha);
    lowMachUpReCl  = interp1(Cn_alpha( lowMachUpReIndex,4),Cn_alpha( lowMachUpReIndex,3),desiredAlpha);
    upMachLowReCl  = interp1(Cn_alpha( upMachLowReIndex,4),Cn_alpha( upMachLowReIndex,3),desiredAlpha);
    upMachUpReCl   = interp1(Cn_alpha(  upMachUpReIndex,4),Cn_alpha(  upMachUpReIndex,3),desiredAlpha);
    lowMach_Cl = mean([lowMachLowReCl,lowMachUpReCl]);
    upMach_Cl  = mean([ upMachLowReCl, upMachUpReCl]);
    Cl = mean([lowMach_Cl,upMach_Cl]);
    
    
%DRAG COEFFICIENT
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_Cd);
    lowMachLowReCd = interp1(Cn_Cd(lowMachLowReIndex,3),Cn_Cd(lowMachLowReIndex,4),Cl);
    lowMachUpReCd  = interp1(Cn_Cd( lowMachUpReIndex,3),Cn_Cd( lowMachUpReIndex,4),Cl);
    upMachLowReCd  = interp1(Cn_Cd( upMachLowReIndex,3),Cn_Cd( upMachLowReIndex,4),Cl);
    upMachUpReCd   = interp1(Cn_Cd(  upMachUpReIndex,3),Cn_Cd(  upMachUpReIndex,4),Cl);
    lowMach_Cd = mean([lowMachLowReCd,lowMachUpReCd]);
    upMach_Cd  = mean([ upMachLowReCd, upMachUpReCd]);
    Cd = mean([lowMach_Cd,upMach_Cd]);
    
    
%PITCHING MOMENT COEFFICIENT
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_Cm);
    lowMachLowReCm = interp1(Cn_Cm(lowMachLowReIndex,3),Cn_Cm(lowMachLowReIndex,4),Cl);
    lowMachUpReCm  = interp1(Cn_Cm( lowMachUpReIndex,3),Cn_Cm( lowMachUpReIndex,4),Cl);
    upMachLowReCm  = interp1(Cn_Cm( upMachLowReIndex,3),Cn_Cm( upMachLowReIndex,4),Cl);
    upMachUpReCm   = interp1(Cn_Cm(  upMachUpReIndex,3),Cn_Cm(  upMachUpReIndex,4),Cl);
    lowMach_Cm = mean([lowMachLowReCm,lowMachUpReCm]);
    upMach_Cm  = mean([ upMachLowReCm, upMachUpReCm]);
    Cm = mean([lowMach_Cm,upMach_Cm]);    
    
end    
    
    
end



function [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,dataArray)

    MachArray = unique(dataArray(:,1));
    [~,MachIndex] = min(abs(MachArray-desiredMach));
    if desiredMach>MachArray(MachIndex)
        if MachIndex == length(MachArray)
            lowMach = MachArray(MachIndex);
            upMach  = MachArray(MachIndex);
        else
            lowMach = MachArray(MachIndex);
            upMach  = MachArray(MachIndex+1);
        end
    elseif desiredMach<MachArray(MachIndex)
        if MachIndex == 1
            lowMach = MachArray(MachIndex);
            upMach  = MachArray(MachIndex);
        else
            lowMach = MachArray(MachIndex-1);
            upMach  = MachArray(MachIndex);
        end
    else
        lowMach = MachArray(MachIndex);
        upMach  = MachArray(MachIndex);
    end


    lowMachIndex = dataArray(:,1)==lowMach;
    upMachIndex  = dataArray(:,1)==upMach;

    lowMachReArray = unique(dataArray(lowMachIndex,2));
    upMachReArray  = unique(dataArray(upMachIndex,2));

    [~,lowMachReIndex] = min(abs(lowMachReArray-desiredReynolds));
    [~, upMachReIndex] = min(abs( upMachReArray-desiredReynolds));

    if desiredReynolds>lowMachReArray(lowMachReIndex)
        if lowMachReIndex == length(lowMachReArray)
            lowMachLowRe = lowMachReArray(lowMachReIndex);
            lowMachUpRe  = lowMachReArray(lowMachReIndex);
        else
            lowMachLowRe = lowMachReArray(lowMachReIndex);
            lowMachUpRe  = lowMachReArray(lowMachReIndex+1);
        end
    elseif desiredReynolds<lowMachReArray(lowMachReIndex)
        if lowMachReIndex == 1
            lowMachLowRe = lowMachReArray(lowMachReIndex);
            lowMachUpRe  = lowMachReArray(lowMachReIndex);
        else
            lowMachLowRe = lowMachReArray(lowMachReIndex-1);
            lowMachUpRe  = lowMachReArray(lowMachReIndex);
        end
    else
        lowMachLowRe = lowMachReArray(lowMachReIndex);
        lowMachUpRe  = lowMachReArray(lowMachReIndex);
    end

    if desiredReynolds>upMachReArray(upMachReIndex)
        if upMachReIndex == length(upMachReArray)
            upMachLowRe = upMachReArray(upMachReIndex);
            upMachUpRe  = upMachReArray(upMachReIndex);
        else
            upMachLowRe = upMachReArray(upMachReIndex);
            upMachUpRe  = upMachReArray(upMachReIndex+1);
        end
    elseif desiredReynolds<upMachReArray(upMachReIndex)
        if upMachReIndex == 1
            upMachLowRe = upMachReArray(upMachReIndex);
            upMachUpRe  = upMachReArray(upMachReIndex);
        else
            upMachLowRe = upMachReArray(upMachReIndex-1);
            upMachUpRe  = upMachReArray(upMachReIndex);
        end
    else
        upMachLowRe = upMachReArray(upMachReIndex);
        upMachUpRe  = upMachReArray(upMachReIndex);
    end


    lowMachLowReIndex = logical(lowMachIndex.*(dataArray(:,2)==lowMachLowRe));
    lowMachUpReIndex  = logical(lowMachIndex.*(dataArray(:,2)==lowMachUpRe));
    upMachLowReIndex  = logical(upMachIndex .*(dataArray(:,2)==upMachLowRe));
    upMachUpReIndex   = logical(upMachIndex .*(dataArray(:,2)==upMachUpRe));

end