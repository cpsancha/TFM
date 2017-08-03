function [Cl_alpha,alpha_zeroLift,Cm_ac] = getAirfoilData(desiredMach,desiredReynolds,varargin)

%CHECK IF EXIST desiredAlpha
if isempty(varargin)
elseif length(varargin)==1
    desiredAlpha = varargin(1);
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
        Cn_Cm  = cat(2,Mach,Re,Cn,Other);


%LIFT CURVE SLOPE & ZERO LIFT ANGLE
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_alpha);
    lowMachLowReFit = polyfit(deg2rad(Cn_alpha(lowMachLowReIndex,4)),Cn_alpha(lowMachLowReIndex,3),1);
    lowMachUpReFit  = polyfit(deg2rad(Cn_alpha( lowMachUpReIndex,4)),Cn_alpha( lowMachUpReIndex,3),1);
    upMachLowReFit  = polyfit(deg2rad(Cn_alpha( upMachLowReIndex,4)),Cn_alpha( upMachLowReIndex,3),1);
    upMachUpReFit   = polyfit(deg2rad(Cn_alpha(  upMachUpReIndex,4)),Cn_alpha(  upMachUpReIndex,3),1);
    %Cl_alpha    
        lowMach_Cl_alpha = mean(lowMachLowReFit(1),lowMachUpReFit(1));
        upMach_Cl_alpha  = mean( upMachLowReFit(1), upMachUpReFit(1));
        Cl_alpha = mean(lowMach_Cl_alpha,upMach_Cl_alpha); 
    %alpha_zeroLift    
        lowMach_alpha_zeroLift = mean(rad2deg(-lowMachLowReFit(2)/lowMachLowReFit(1)),rad2deg(-lowMachUpReFit(2)/lowMachUpReFit(1)));
        upMach_alpha_zeroLift  = mean(rad2deg(- upMachLowReFit(2)/ upMachLowReFit(1)),rad2deg(- upMachUpReFit(2)/ upMachUpReFit(1)));
        alpha_zeroLift = mean(lowMach_alpha_zeroLift, upMach_alpha_zeroLift); 


%AERODYNAMIC CENTER MOMENTUM COEFFICIENT
    [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,Cn_Cm);
    lowMachLowReCm_ac = mean(Cn_Cm(lowMachLowReIndex,4));
    lowMachUpReCm_ac  = mean(Cn_Cm( lowMachUpReIndex,4));
    upMachLowReCm_ac  = mean(Cn_Cm( upMachLowReIndex,4));
    upMachUpReCm_ac   = mean(Cn_Cm(  upMachUpReIndex,4));
    lowMach_Cm_ac = mean(lowMachLowReCm_ac,lowMachUpReCm_ac);
    upMach_Cm_ac  = mean( upMachLowReCm_ac, upMachUpReCm_ac);
    Cm_ac = mean(lowMach_Cm_ac,upMach_Cm_ac); 
        
%

end



function [lowMachLowReIndex,lowMachUpReIndex,upMachLowReIndex,upMachUpReIndex] = getMachReIndexes(desiredMach,desiredReynolds,dataArray)

    MachArray = unique(dataArray(:,1));
    [~,MachIndex] = min(abs(MachArray-desiredMach));
    if desiredMach>MachArray(MachIndex)
        lowMach = MachArray(MachIndex);
        upMach  = MachArray(MachIndex+1);
    elseif desiredMach<MachArray(MachIndex)
        lowMach = MachArray(MachIndex-1);
        upMach  = MachArray(MachIndex);
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
        lowMachLowRe = lowMachReArray(lowMachReIndex);
        lowMachUpRe  = lowMachReArray(lowMachReIndex+1);
    elseif desiredReynolds<lowMachReArray(lowMachReIndex)
        lowMachLowRe = lowMachReArray(lowMachReIndex-1);
        lowMachUpRe  = lowMachReArray(lowMachReIndex);
    else
        lowMachLowRe = lowMachReArray(lowMachReIndex);
        lowMachUpRe  = lowMachReArray(lowMachReIndex);
    end

    if desiredReynolds>upMachReArray(upMachReIndex)
        upMachLowRe = upMachReArray(upMachReIndex);
        upMachUpRe  = upMachReArray(upMachReIndex+1);
    elseif desiredReynolds<upMachReArray(upMachReIndex)
        upMachLowRe = upMachReArray(upMachReIndex-1);
        upMachUpRe  = upMachReArray(upMachReIndex);
    else
        upMachLowRe = upMachReArray(upMachReIndex);
        upMachUpRe  = upMachReArray(upMachReIndex);
    end


    lowMachLowReIndex = lowMachIndex.*(dataArray(:,2)==lowMachLowRe);
    lowMachUpReIndex  = lowMachIndex.*(dataArray(:,2)==lowMachUpRe);
    upMachLowReIndex  = upMachIndex .*(dataArray(:,2)==upMachLowRe);
    upMachUpReIndex   = upMachIndex .*(dataArray(:,2)==upMachUpRe);

end