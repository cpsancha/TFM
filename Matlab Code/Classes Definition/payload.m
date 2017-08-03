classdef payload < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        beds      double        %[double]  N� of seats that can be converted into sleeping beds [-]
        crew      double        %[double]  N� of necessary crew members [-]
        paxMax    double        %[double]  N� of passengers in high-density seating [-]
        paxMin    double        %[double]  N� of passengers in executive seating [-]
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = payload(varargin)
        end
                     
    end
    
    
end

