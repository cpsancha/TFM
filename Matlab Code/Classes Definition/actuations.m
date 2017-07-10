classdef actuations < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        Range       double        %[double]  Range [km]
        Endurance   double        %[double]  Endurance [hr]
        MMO         double        %[double]  Max. Mach Operating [-]
        Vmax        double        %[double]  Max speed [m/s]
        Mcruise     double        %[double]  Cruise Mach [-]
        Vcruise     double        %[double]  Cruise speed [m/s]
        L_D         double        %[double]  Efficiency [-]
        Vto         double        %[double]  Take Off speed [m/s]
        Vstall      double        %[double]  Stall Speed in clean configuration [m/s]
        Vstall_TO   double        %[double]  Stall Speed in Take-Off configuration [m/s]
        Vstall_L    double        %[double]  Stall Speed in Landing configuration[m/s]
        Vapprox     double        %[double]  Approximation Speed (1.3*Vstall_L) [m/s]
        Vasc        double        %[double]  Ascensional Speed [m/s]
        Hmax        double        %[double]  Max. Operating Altitude [m]
        Hcruise     double        %[double]  Typical cruise Altitude [m]
        Sto         double        %[double]  Take off distance [m]
        Sl          double        %[double]  Landing distance [m]
        
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = actuations(varargin)
        end
                     
    end
    
    
end

