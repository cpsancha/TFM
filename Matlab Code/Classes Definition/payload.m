classdef payload < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        beds      double        %[double]  Nº of seats that can be converted into sleeping beds [-]
        crew      double        %[double]  Nº of necessary crew members [-]
        paxMax    double        %[double]  Nº of passengers in high-density seating [-]
        paxMin    double        %[double]  Nº of passengers in executive seating [-]
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = payload(varargin)
        end
                     
    end
    
    
end

