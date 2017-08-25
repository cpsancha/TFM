classdef fuselage < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        A_I            double    % Planform area of section I [m^2] Figure F.13 Torenbeek
        A_II           double    % Planform area of section II [m^2]
        beta           double    % Upsweep angle of the afterbody [º]
        bagVolume      double    %[double]  Total volume of the baggage [m^3]
        cabHeight      double    %[double]  Height of the cabin [m]
        cabLength      double    %[double]  Net length of the cabin excluding cockpit and baggage (includes gallery) [m]
        cabVolume      double    %[double]  Total volume of the cabin [m^3]
        cabWidth       double    %[double]  Net wide of the cabin [m]
        cabPos         double    %longitudinal position of start of passenger cabin [m]
        finenessRatio  double    %[double]  The ratio of the length of a body to its maximum width (Esbeltez) [-]
        frontArea      double    %[double]  Maximum front area [m^2]
        fuselage_AoA   double    %[double]  Angle of Attack of the fuselage [º]
        fusHeight      double    %[double]  Maximum height of the fuselage [m]
        fusHeightWidth double    %[double]  The ratio between the fuselage height and its width [-]    
        fusLength      double    %[double]  Total length of the fuselage [m]
        fusWidth       double    %[double]  Maximum width of the fuselage [m]
        Swet           double    % exterior wet surface [m^2]
        la             double    % length of the afterbody [m]
        ln             double    % length of the nose [m]
        minHeight      double    %[double]  Height of the lower point of the fuselage till the ground [m]  
        Volume         double    %[double]  Total exterior volume of the fuselage [m^3]
        cabinFrac      double    % Cockpit position in fraction of fuselage lenght( for xcg position) [-]
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = fuselage(varargin)
        end
                     
    end
    
    
end

