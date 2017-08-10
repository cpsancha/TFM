classdef fuselage < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        bagVolume      double    %[double]  Total volume of the baggage [m^3]
        cabHeight      double    %[double]  Height of the cabin [m]
        cabLength      double    %[double]  Net length of the cabin excluding cockpit and baggage [m]
        cabVolume      double    %[double]  Total volume of the cabin [m^3]
        cabWidth       double    %[double]  Net wide of the cabin [m]
        finenessRatio  double    %[double]  The ratio of the length of a body to its maximum width (Esbeltez) [-]
        frontArea      double    %[double]  Maximum front area [m^2]
        fuselage_AoA   double    %[double]  Angle of Attack of the fuselage [º]
        fusHeight      double    %[double]  Maximum height of the fuselage [m]
        fusHeightWidth double    %[double]  The ratio between the fuselage height and its width [-]    
        fusLength      double    %[double]  Total length of the fuselage [m]
        fusWidth       double    %[double]  Maximum width of the fuselage [m]
        minHeight      double    %[double]  Height of the lower point of the fuselage till the ground [m]  
        Volume         double    %[double]  Total exterior volume of the fuselage [m^3]
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = fuselage(varargin)
        end
                     
    end
    
    
end

