function setPath()
%SETPATH Summary of this function goes here
%   Detailed explanation goes here

%Get current folder
sr=which(mfilename);
i=max(strfind(lower(sr),lower(mfilename)))-2;
if i>0
  sr=sr(1:i);
else
  error('Cannot locate setPath.m. You must add the path to the MTorres directory.')
end

%Define existing folders
%Use the following if subfolder --> fullfile(sr,'Control',filesep,'DynamicSimulation',...
sd  =   {...
        fullfile(sr,'Classes Definition'),...
        fullfile(sr,'Configuration Design'),...
        fullfile(sr,'Digitalized Data'),...
        fullfile(sr,'Preliminary Sizing'),...
        fullfile(sr,'Similar Planes'),...
        fullfile(sr,'Utilities')        
        };
        
%Add defined folders to path        
for i=1:length(sd)
    addpath(sd{i});
end









end

