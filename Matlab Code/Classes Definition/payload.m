classdef payload < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        crew      double        %[double]  Nº of necessary crew members [-]
        paxMin    double        %[double]  Nº of passengers in executive seating [-]
        paxMax    double        %[double]  Nº of passengers in high-density seating [-]
        beds      double        %[double]  Nº of seats that can be converted into sleeping beds [-]
        cabLength double        %[double]  Net length of the cabin excluding cockpit and baggage [m]
        cabWide   double        %[double]  Net wide of the cabin [m]
        cabHeight double        %[double]  Height of the cabin [m]
        cabVolume double        %[double]  Total volume of the cabin [m^3]
        bagVolume double        %[double]  Total volume of the baggage [m^3]
 
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = payload(varargin)
        end
                     
    end
    
    
end

