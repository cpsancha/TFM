classdef fuselage < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        fusLength      double    %[double]  Total length of the fuselage [m]
        fusWidth       double    %[double]  Maximum width of the fuselage [m]
        fusHeight      double    %[double]  Maximum height of the fuselage [m]
        cabLength      double    %[double]  Net length of the cabin excluding cockpit and baggage [m]
        cabWidth       double    %[double]  Net wide of the cabin [m]
        cabHeight      double    %[double]  Height of the cabin [m]
        cabVolume      double    %[double]  Total volume of the cabin [m^3]
        bagVolume      double    %[double]  Total volume of the baggage [m^3]
        frontArea      double    %[double]  Maximum front area [m^2]
        minHeight      double    %[double]  Height of the lower point of the fuselage till the ground [m]
        finenessRatio  double    %[double]  The ratio of the length of a body to its maximum width (Esbeltez) [-]
        fusHeightWidth double    %[double]  The ratio between the fuselage height and its width [-]        
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = fuselage(varargin)
        end
                     
    end
    
    
end

