function [SimilarPlanes, MV] = importSimilarPlanes( type, CST, CF )
%IMPORTSIMILARPLANES Summary of this function goes here
%   Detailed explanation goes here


%% Get Excel Path
sr=which(mfilename);
i=max(strfind(lower(sr),lower('MTORRES')))+6;
if i>0
  sr=sr(1:i);
else
  error('Cannot locate MTORRES directory. You must add the path to the MTorres directory.')
end

switch type
    case 5
        fileSP = fullfile(sr,'Aviones Semejantes',filesep,'Aviones_Semejantes.xlsx');
    case 11
        fileSP = fullfile(sr,'Aviones Semejantes',filesep,'Aviones_Semejantes_Seaplanes.xlsx');
        
end



%% Select Mission Specification
switch type
    case 5
        %Options
        numberSP = 11;
        sheetSP  = 'Aviones Semejantes Long-Range';
        initLetter = 'E';
        endLetter  = native2unicode(unicode2native(initLetter)+numberSP-1);
        
        %Load Plane Data Data from Excel
        excelPlane = importFile(fileSP, sheetSP,strcat(initLetter,'3:',endLetter,'4'));
        
        %Load General Data from Excel
        excelGeneralData = importFile(fileSP, sheetSP,strcat(initLetter,'6:',endLetter,'9'));
                
        %Load Engines from Excel
        excelEngines = importFile(fileSP, sheetSP,strcat(initLetter,'11:',endLetter,'23'));
                
        %Load Weights from Excel
        excelWeights = importFile(fileSP, sheetSP,strcat(initLetter,'25:',endLetter,'34'));
        
        %Load Payload from Excel
        excelPayload = importFile(fileSP, sheetSP,strcat(initLetter,'41:',endLetter,'49'));
        
        %Load Wing from Excel
        excelWing = importFile(fileSP, sheetSP,strcat(initLetter,'63:',endLetter,'83'));
        
        %Load Actuations from Excel
        excelActuations = importFile(fileSP, sheetSP,strcat(initLetter,'159:',endLetter,'175'));
        
    case 11
                %Options
        numberSP = 7;
        sheetSP  = 'Aviones Semejantes';
        initLetter = 'E';
        endLetter  = native2unicode(unicode2native(initLetter)+numberSP-1);
        
        %Load Plane Data Data from Excel
        excelPlane = importFile(fileSP, sheetSP,strcat(initLetter,'3:',endLetter,'4'));
        
        %Load General Data from Excel
        excelGeneralData = importFile(fileSP, sheetSP,strcat(initLetter,'6:',endLetter,'9'));
                
        %Load Engines from Excel
        excelEngines = importFile(fileSP, sheetSP,strcat(initLetter,'11:',endLetter,'24'));
                
        %Load Weights from Excel
        excelWeights = importFile(fileSP, sheetSP,strcat(initLetter,'26:',endLetter,'41'));
        
        %Load Payload from Excel
        excelPayload = importFile(fileSP, sheetSP,strcat(initLetter,'43:',endLetter,'54'));
        
        %Load Wing from Excel
        excelWing = importFile(fileSP, sheetSP,strcat(initLetter,'65:',endLetter,'85'));
        
        %Load Hull from Excel
        excelHull = importFile(fileSP, sheetSP,strcat(initLetter,'147:',endLetter,'149'));
        
        %Load Actuations from Excel
        excelActuations = importFile(fileSP, sheetSP,strcat(initLetter,'165:',endLetter,'187'));
    otherwise
        error('There are only two cases, businessJet (5) or Amphibious (11). Choose one.')
end

%% Create Similar Planes Structure
SimilarPlanes = cell(numberSP,1);
switch type
    case 5
        for i=1:numberSP
            %Create aircraft
            SimilarPlanes{i} = aircraft();

            %Plane Name
            SimilarPlanes{i}.Model        = string(excelPlane{1,i});
            SimilarPlanes{i}.Manufacturer = string(excelPlane{2,i});

            %General Data
            SimilarPlanes{i}.FirstFlight = string(excelGeneralData{1,i});
            SimilarPlanes{i}.Height      = double(excelGeneralData{2,i});
            SimilarPlanes{i}.Length      = double(excelGeneralData{3,i});
            SimilarPlanes{i}.Wingspan    = double(excelGeneralData{4,i});

            %Weights
            SimilarPlanes{i}.Weight.MTOW = double(excelWeights{1,i});
            SimilarPlanes{i}.Weight.MRW  = double(excelWeights{2,i});
            SimilarPlanes{i}.Weight.EW   = double(excelWeights{3,i});
            SimilarPlanes{i}.Weight.OEW  = double(excelWeights{4,i});
            SimilarPlanes{i}.Weight.BOW  = double(excelWeights{5,i});
            SimilarPlanes{i}.Weight.MPL  = double(excelWeights{6,i});
            SimilarPlanes{i}.Weight.MFW  = double(excelWeights{7,i});
            SimilarPlanes{i}.Weight.TUL  = double(excelWeights{8,i});
            SimilarPlanes{i}.Weight.MZFW = double(excelWeights{9,i});
            SimilarPlanes{i}.Weight.MLW  = double(excelWeights{10,i});

            %Engines
            SimilarPlanes{i}.Engine.Number       = double(excelEngines{1,i});
            SimilarPlanes{i}.Engine.PositionStr  = string(excelEngines{2,i});
            SimilarPlanes{i}.Engine.Type         = string(excelEngines{3,i});
            SimilarPlanes{i}.Engine.Manufacturer = string(excelEngines{4,i});
            SimilarPlanes{i}.Engine.Model        = string(excelEngines{5,i});
            SimilarPlanes{i}.Engine.Weight       = double(excelEngines{6,i});
            SimilarPlanes{i}.Engine.Thrust       = double(excelEngines{7,i})*1e3; %kN --> N
            SimilarPlanes{i}.Engine.TSFC         = double(excelEngines{8,i});
            SimilarPlanes{i}.Engine.SFC          = double(excelEngines{9,i});
            SimilarPlanes{i}.Engine.etaPropeller = double(excelEngines{10,i});
            SimilarPlanes{i}.Engine.Diameter     = double(excelEngines{11,i});
            SimilarPlanes{i}.Engine.Length       = double(excelEngines{12,i});
            SimilarPlanes{i}.Engine.TotalThrust  = double(excelEngines{13,i})*1e3; %kN --> N
            SimilarPlanes{i}.Engine.TotalWeight  = double(excelEngines{1,i})*...
                                                   double(excelEngines{6,i});
            
            %Payload
            SimilarPlanes{i}.Payload.crew      = double(excelPayload{1,i});
            SimilarPlanes{i}.Payload.paxMin    = double(excelPayload{2,i});
            SimilarPlanes{i}.Payload.paxMax    = double(excelPayload{3,i});
            SimilarPlanes{i}.Payload.beds      = double(excelPayload{4,i});
            SimilarPlanes{i}.Payload.cabLength = double(excelPayload{5,i});
            SimilarPlanes{i}.Payload.cabWide   = double(excelPayload{6,i});
            SimilarPlanes{i}.Payload.cabHeight = double(excelPayload{7,i});
            SimilarPlanes{i}.Payload.cabVolume = double(excelPayload{8,i});
            SimilarPlanes{i}.Payload.bagVolume = double(excelPayload{9,i});

            %Wing
            SimilarPlanes{i}.Wing.Sw           = double(excelWing{4,i});
            SimilarPlanes{i}.Wing.WingSpan     = double(excelWing{5,i});
            SimilarPlanes{i}.Wing.RealSemiSpan = double(excelWing{6,i});
            SimilarPlanes{i}.Wing.TipChord     = double(excelWing{7,i});
            SimilarPlanes{i}.Wing.RootChord    = double(excelWing{8,i});
            SimilarPlanes{i}.Wing.CMG          = double(excelWing{9,i});
            SimilarPlanes{i}.Wing.CMA          = double(excelWing{10,i});
            SimilarPlanes{i}.Wing.AspectRatio  = double(excelWing{15,i});
            SimilarPlanes{i}.Wing.TaperRatio   = double(excelWing{16,i});
            SimilarPlanes{i}.Wing.Sweep        = double(excelWing{17,i});
            SimilarPlanes{i}.Wing.Dihedral     = double(excelWing{18,i});
            SimilarPlanes{i}.Wing.Airfoil      = string(excelWing{19,i});
            SimilarPlanes{i}.Wing.WingLoading  = double(excelWing{21,i});
            SimilarPlanes{i}.Wing.LongPos      = double(excelWing{2,i});
            SimilarPlanes{i}.Wing.Root_LE      = double(excelWing{3,i});
            SimilarPlanes{i}.Wing.CMA_LE       = double(excelWing{11,i});
            SimilarPlanes{i}.Wing.CMA_14       = double(excelWing{12,i});
            SimilarPlanes{i}.Wing.CMA_b        = double(excelWing{13,i});
            SimilarPlanes{i}.Wing.TipSweep     = double(excelWing{14,i});
        %     SimilarPlanes{i}.Wing.CLmax_Cr     = double(excelWing{17,i});
        %     SimilarPlanes{i}.Wing.CLmaxTO      = double(excelWing{17,i});
        %     SimilarPlanes{i}.Wing.CLmaxL       = double(excelWing{17,i});

            %Actuations
            SimilarPlanes{i}.Actuations.Vmax      = double(excelActuations{1,i});
            SimilarPlanes{i}.Actuations.MMO       = double(excelActuations{2,i});
            SimilarPlanes{i}.Actuations.Vcruise   = double(excelActuations{3,i});
            SimilarPlanes{i}.Actuations.Mcruise   = double(excelActuations{4,i});
            SimilarPlanes{i}.Actuations.Vstall    = double(excelActuations{5,i});
            SimilarPlanes{i}.Actuations.Vstall_TO = double(excelActuations{6,i});
            SimilarPlanes{i}.Actuations.Vstall_L  = double(excelActuations{7,i});
            SimilarPlanes{i}.Actuations.Vto       = double(excelActuations{8,i});
            SimilarPlanes{i}.Actuations.Vapprox   = double(excelActuations{9,i});
            SimilarPlanes{i}.Actuations.Vasc      = double(excelActuations{10,i});
            SimilarPlanes{i}.Actuations.Range     = double(excelActuations{11,i});
            SimilarPlanes{i}.Actuations.Endurance = double(excelActuations{12,i});
            SimilarPlanes{i}.Actuations.Hmax      = double(excelActuations{13,i});
            SimilarPlanes{i}.Actuations.Hcruise   = double(excelActuations{14,i});
            SimilarPlanes{i}.Actuations.Sto       = double(excelActuations{16,i});
            SimilarPlanes{i}.Actuations.Sl        = double(excelActuations{17,i});

        end

    case 11
        for i=1:numberSP
            %Create aircraft
            SimilarPlanes{i} = aircraft();

            %Plane Name
            SimilarPlanes{i}.Model        = string(excelPlane{1,i});
            SimilarPlanes{i}.Manufacturer = string(excelPlane{2,i});

            %General Data
            SimilarPlanes{i}.FirstFlight = string(excelGeneralData{1,i});
            SimilarPlanes{i}.Height      = double(excelGeneralData{2,i});
            SimilarPlanes{i}.Length      = double(excelGeneralData{3,i});
            SimilarPlanes{i}.Wingspan    = double(excelGeneralData{4,i});

            %Weights
            SimilarPlanes{i}.Weight.MTOW = double(excelWeights{1,i});
            SimilarPlanes{i}.Weight.MRW  = double(excelWeights{2,i});
            SimilarPlanes{i}.Weight.EW   = double(excelWeights{3,i});
            SimilarPlanes{i}.Weight.OEW  = double(excelWeights{4,i});
            SimilarPlanes{i}.Weight.BOW  = double(excelWeights{5,i});
            SimilarPlanes{i}.Weight.MPL  = double(excelWeights{6,i});
            SimilarPlanes{i}.Weight.MFW  = double(excelWeights{7,i});
            SimilarPlanes{i}.Weight.TUL  = double(excelWeights{8,i});
            SimilarPlanes{i}.Weight.MZFW = double(excelWeights{9,i});
            SimilarPlanes{i}.Weight.MLW  = double(excelWeights{10,i});

            %Engines
            SimilarPlanes{i}.Engine.Number       = double(excelEngines{1,i});
            SimilarPlanes{i}.Engine.PositionStr  = string(excelEngines{2,i});
            SimilarPlanes{i}.Engine.Type         = string(excelEngines{3,i});
            SimilarPlanes{i}.Engine.Manufacturer = string(excelEngines{4,i});
            SimilarPlanes{i}.Engine.Model        = string(excelEngines{5,i});
            SimilarPlanes{i}.Engine.Weight       = double(excelEngines{6,i});
            SimilarPlanes{i}.Engine.Thrust       = double(excelEngines{7,i});
            SimilarPlanes{i}.Engine.TSFC         = double(excelEngines{8,i});
            SimilarPlanes{i}.Engine.SFC          = double(excelEngines{9,i});
            SimilarPlanes{i}.Engine.etaPropeller = double(excelEngines{10,i});
            SimilarPlanes{i}.Engine.Diameter     = double(excelEngines{11,i});
            SimilarPlanes{i}.Engine.Length       = double(excelEngines{12,i});

            %Payload
            SimilarPlanes{i}.Payload.crew      = double(excelPayload{1,i});
            SimilarPlanes{i}.Payload.paxMin    = double(excelPayload{2,i});
            SimilarPlanes{i}.Payload.paxMax    = double(excelPayload{3,i});
            SimilarPlanes{i}.Payload.beds      = double(excelPayload{4,i});
            SimilarPlanes{i}.Payload.cabLength = double(excelPayload{5,i});
            SimilarPlanes{i}.Payload.cabWide   = double(excelPayload{6,i});
            SimilarPlanes{i}.Payload.cabHeight = double(excelPayload{7,i});
            SimilarPlanes{i}.Payload.cabVolume = double(excelPayload{8,i});
            SimilarPlanes{i}.Payload.bagVolume = double(excelPayload{9,i});

            %Wing
            SimilarPlanes{i}.Wing.Sw           = double(excelWing{4,i});
            SimilarPlanes{i}.Wing.WingSpan     = double(excelWing{5,i});
%             SimilarPlanes{i}.Wing.RealSemiSpan = double(excelWing{6,i});
            SimilarPlanes{i}.Wing.TipChord     = double(excelWing{7,i});
            SimilarPlanes{i}.Wing.RootChord    = double(excelWing{8,i});
            SimilarPlanes{i}.Wing.CMG          = double(excelWing{9,i});
            SimilarPlanes{i}.Wing.CMA          = double(excelWing{10,i});
            SimilarPlanes{i}.Wing.AspectRatio  = double(excelWing{15,i});
            SimilarPlanes{i}.Wing.TaperRatio   = double(excelWing{16,i});
%             SimilarPlanes{i}.Wing.Sweep        = double(excelWing{17,i});
%             SimilarPlanes{i}.Wing.Dihedral     = double(excelWing{18,i});
%             SimilarPlanes{i}.Wing.Airfoil      = string(excelWing{19,i});
%             SimilarPlanes{i}.Wing.WingLoading  = double(excelWing{21,i});
%             SimilarPlanes{i}.Wing.LongPos      = double(excelWing{2,i});
%             SimilarPlanes{i}.Wing.Root_LE      = double(excelWing{3,i});
%             SimilarPlanes{i}.Wing.CMA_LE       = double(excelWing{11,i});
%             SimilarPlanes{i}.Wing.CMA_14       = double(excelWing{12,i});
%             SimilarPlanes{i}.Wing.CMA_b        = double(excelWing{13,i});
%             SimilarPlanes{i}.Wing.TipSweep     = double(excelWing{14,i});
            SimilarPlanes{i}.Wing.WingLoading  = double(excelWing{21,i});
   
            
        %     SimilarPlanes{i}.Wing.CLmax     = double(excelWing{17,i});
        %     SimilarPlanes{i}.Wing.CLmax_TO      = double(excelWing{17,i});
        %     SimilarPlanes{i}.Wing.CLmax_L       = double(excelWing{17,i});

            %Actuations
            SimilarPlanes{i}.Actuations.Vmax      = double(excelActuations{1,i});
            SimilarPlanes{i}.Actuations.MMO       = double(excelActuations{2,i});
            SimilarPlanes{i}.Actuations.Vcruise   = double(excelActuations{3,i});
            SimilarPlanes{i}.Actuations.Mcruise   = double(excelActuations{4,i});
            SimilarPlanes{i}.Actuations.Vstall    = double(excelActuations{5,i});
            SimilarPlanes{i}.Actuations.Vstall_TO = double(excelActuations{6,i});
            SimilarPlanes{i}.Actuations.Vstall_TOw= double(excelActuations{7,i});
            SimilarPlanes{i}.Actuations.Vstall_L  = double(excelActuations{8,i});
            SimilarPlanes{i}.Actuations.Vstall_Lw = double(excelActuations{9,i});
            SimilarPlanes{i}.Actuations.Vto       = double(excelActuations{10,i});
            SimilarPlanes{i}.Actuations.Vtow      = double(excelActuations{11,i});            
            SimilarPlanes{i}.Actuations.Vl       = double(excelActuations{12,i});           
            SimilarPlanes{i}.Actuations.Vlw       = double(excelActuations{13,i});
            SimilarPlanes{i}.Actuations.Vapprox   = double(excelActuations{14,i});
            SimilarPlanes{i}.Actuations.Vasc      = double(excelActuations{15,i});
            SimilarPlanes{i}.Actuations.Range     = double(excelActuations{16,i});
            SimilarPlanes{i}.Actuations.Endurance = double(excelActuations{17,i});
            SimilarPlanes{i}.Actuations.Hmax      = double(excelActuations{18,i});
            SimilarPlanes{i}.Actuations.Hcruise   = double(excelActuations{19,i});
            SimilarPlanes{i}.Actuations.Sto       = double(excelActuations{20,i});
            SimilarPlanes{i}.Actuations.Sl        = double(excelActuations{21,i});
            SimilarPlanes{i}.Actuations.Stow       = double(excelActuations{22,i});            
            SimilarPlanes{i}.Actuations.Slw        = double(excelActuations{23,i});

        end

end



%% Additional Corrections
switch type
    case 5
        for i=1:numberSP
            %Calculate efficiency
            if (~isempty(SimilarPlanes{i}.Actuations.Range)   && ...
                ~isempty(SimilarPlanes{i}.Engine.TSFC)        && ...
                ~isempty(SimilarPlanes{i}.Actuations.Vcruise) && ...
                ~isempty(SimilarPlanes{i}.Weight.MTOW)        && ...
                ~isempty(SimilarPlanes{i}.Weight.MFW))    
            SimilarPlanes{i}.Actuations.L_D = (SimilarPlanes{i}.Actuations.Range*1e3*...
                SimilarPlanes{i}.Engine.TSFC*CF.TSFC2SI*CST.GravitySI)/(SimilarPlanes{i}...
                .Actuations.Vcruise*log(SimilarPlanes{i}.Weight.MTOW/(SimilarPlanes{i}...
                .Weight.MTOW-SimilarPlanes{i}.Weight.MFW)));
            end
            
            %Calculate Wing CLmax from WingLoading, Hcruise and Vcruise
            if (~isempty(SimilarPlanes{i}.Wing.WingLoading)  && ...
                ~isempty(SimilarPlanes{i}.Actuations.Hcruise) && ...
                ~isempty(SimilarPlanes{i}.Actuations.Vcruise))
                    [~,~,~,rho] = atmosisa(SimilarPlanes{i}.Actuations.Hcruise);
                    SimilarPlanes{i}.Wing.CLmax = 2*SimilarPlanes{i}.Wing.WingLoading*CST.GravitySI/...
                                                 (rho*SimilarPlanes{i}.Actuations.Vcruise^2);
            end
            
            %Calculate Wing CLmax_TO from Vstall_TO and MTOW (WingLoading)
            if (~isempty(SimilarPlanes{i}.Wing.WingLoading)  && ...
                ~isempty(SimilarPlanes{i}.Actuations.Vstall_TO))
                    [~,~,~,rho] = atmosisa(0);
                    SimilarPlanes{i}.Wing.CLmax_TO = 2*SimilarPlanes{i}.Wing.WingLoading*CST.GravitySI/...
                                                 (rho*SimilarPlanes{i}.Actuations.Vstall_TO^2);
            end

            %Calculate Wing CLmax_L from Vstall_L and MLW
            if (~isempty(SimilarPlanes{i}.Wing.WingLoading)  && ...
                ~isempty(SimilarPlanes{i}.Actuations.Vstall_L))
                    [~,~,~,rho] = atmosisa(0);
                    SimilarPlanes{i}.Wing.CLmax_L = 2*(SimilarPlanes{i}.Weight.MLW/SimilarPlanes{i}.Wing.Sw)*CST.GravitySI/...
                                                 (rho*SimilarPlanes{i}.Actuations.Vstall_L^2);
            end

        end
        
    case 11
           for i=1:numberSP

                %Calculate Wing CLmax from WingLoading, Hcruise and Vcruise
                if (~isempty(SimilarPlanes{i}.Wing.WingLoading)  && ...
                    ~isempty(SimilarPlanes{i}.Actuations.Hcruise) && ...
                    ~isempty(SimilarPlanes{i}.Actuations.Vcruise))
                        [~,~,~,rho] = atmosisa(SimilarPlanes{i}.Actuations.Hcruise);
                        SimilarPlanes{i}.Wing.CLmax = 2*SimilarPlanes{i}.Wing.WingLoading/...
                                                     (rho*SimilarPlanes{i}.Actuations.Vcruise^2);
                end

                %Calculate Wing CLmax_TO from Vstall_TO and MTOW (WingLoading)
                if (~isempty(SimilarPlanes{i}.Wing.WingLoading)  && ...
                    ~isempty(SimilarPlanes{i}.Actuations.Vstall_TOw))
                        [~,~,~,rho] = atmosisa(0);
                        SimilarPlanes{i}.Wing.CLmax_TO = 2*SimilarPlanes{i}.Wing.WingLoading/...
                                                     (rho*SimilarPlanes{i}.Actuations.Vstall_TOw^2);
                end


                %Calculate Wing CLmax_L from Vstall_L and MLW
                if (~isempty(SimilarPlanes{i}.Wing.WingLoading)  && ...
                    ~isempty(SimilarPlanes{i}.Actuations.Vstall_Lw))
                        [~,~,~,rho] = atmosisa(0);
                        SimilarPlanes{i}.Wing.CLmax_L = 2*(SimilarPlanes{i}.Weight.MLW/SimilarPlanes{i}.Wing.Sw)*CST.GravitySI/...
                                                     (rho*SimilarPlanes{i}.Actuations.Vstall_Lw^2);
                end

           end

        
       
        
end


%% Get mean values from similar planes
switch type
    case 5
    
    case 11  
    SP=SimilarPlanes;
        for i=1:length(SP)
                %CLmax
                if ~isempty(SP{i}.Wing.CLmax)
                    index(i)   = true; %#ok<SAGROW>
                    CLmax(i) = SP{i}.Wing.CLmax; %#ok<SAGROW>
                else
                    index(i)   = false; %#ok<SAGROW>
                    CLmax(i) = NaN; %#ok<SAGROW>
                end  
                %CLmax_TO
                if ~isempty(SP{i}.Wing.CLmax_TO)
                    index(i)   = true; %#ok<SAGROW>
                    CLmax_TO(i) = SP{i}.Wing.CLmax_TO; %#ok<SAGROW>
                else
                    index(i)   = false; %#ok<SAGROW>
                    CLmax_TO(i) = NaN; %#ok<SAGROW>
                end   
                %CLmax_L
                if ~isempty(SP{i}.Wing.CLmax_L)
                    index(i)   = true; %#ok<SAGROW>
                    CLmax_L(i) = SP{i}.Wing.CLmax_L; %#ok<SAGROW>
                else
                    index(i)   = false; %#ok<SAGROW>
                    CLmax_L(i) = NaN; %#ok<SAGROW>
                end 
        end
        
            MV.CL_max = mean(CLmax(index));
            MV.CL_max_TO = mean(CLmax_TO(index));
            MV.CL_max_L = mean(CLmax_L(index));

            clear index i CLmax CLmax_TO CLmax_L



end

