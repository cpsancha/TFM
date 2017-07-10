function CST = importConstants()
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
fileSP = fullfile(sr,'Aviones Semejantes',filesep,'Aviones_Semejantes.xlsx');



%% Load data from Excel
    sheetCST  = 'Constants';
    ConstNum  = length(xlsread(fileSP,sheetCST));
    Constants = importFile(fileSP, sheetCST,strcat('A1:C',num2str(ConstNum)));


%% Create Constants (CST) Structure
CST = struct();
for i=1:ConstNum
    CST.(char(Constants{i,1})) = double(Constants{i,3});   
end


%% Additional Fields




end

