classdef actuations < handle
    %ACTUATIONS A MATLAB class to contain the elementary actuations of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable)
        Endurance   double        %[double]  Endurance [hr]
        Hcruise     double        %[double]  Typical cruise Altitude [m]
        Hmax        double        %[double]  Max. Operating Altitude [m]
        L_D         double        %[double]  Efficiency [-]
        Mcruise     double        %[double]  Cruise Mach [-]
        MMO         double        %[double]  Max. Mach Operating [-]
        Range       double        %[double]  Range [km]
        Sl          double        %[double]  Landing distance [m]
        Slw         double        %[double]  Landing distance water [m]
        Sto         double        %[double]  Take off distance [m]
        Stow        double        %[double]  Take off distance water [m]
        Vapprox     double        %[double]  Approximation Speed (1.3*Vstall_L) [m/s]
        Vasc        double        %[double]  Ascensional Speed [m/s]
        Vcruise     double        %[double]  Cruise speed [m/s]
        Vl          double
        Vlw         double
        Vmax        double        %[double]  Max speed [m/s]
        Vstall      double        %[double]  Stall Speed in clean configuration [m/s]
        Vstall_L    double        %[double]  Stall Speed in Landing configuration[m/s]
        Vstall_Lw   double       %[double]  Stall Speed in Landing configuration[m/s]
        Vstall_TO   double        %[double]  Stall Speed in Take-Off configuration [m/s]   
        Vstall_TOw  double       %[double]  Stall Speed in Take-Off water [m/s]
        Vto         double        %[double]  Take Off speed [m/s] también llamada Lift-Off speed o V2 (1.2*Vstall_TO)
        Vtow        double
        Wf_Wto      double        % Fuel used to get the range in actuations divided by Wto
    end
    
    
    properties (Dependent, GetObservable) 
    end
   
    
    methods
        % Class constructor
        function obj = actuations(varargin)
        end
                     
    end
    
    
end

