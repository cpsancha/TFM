classdef engines < handle
    %ENGINES A MATLAB class to contain the elementary data of an aircraft
    %engines.
    %   Detailed explanation goes here
    
    properties (SetObservable)
        Diameter     double     %[double] Maximum diameter of the engine [m]
        propDiameter double     % Propeller diameter [m]
        etaPropeller double     %[double] Propeller efficiency [-]
        Length       double     %[double] Maximum length of the engine,inculding propeller [m]
        Manufacturer string     %[string] Engine's manufacturer.
        Model        string     %[string] Engine's model.
        Number       double     %[double] Number of engines in the plane.
        Position     double     %[double] Position of the CoG (x y z) of the engine [m].
        PositionStr  string     %[string] Verbose Position of the engine/engines
        Power        double     %[double] Power at take-off [W]   
        Pylon_t_c    double     %[double] Relative width of the pylon airfoil [-]
        Pylon_Swet   double     %[double] Pylon total wet surface [m^2]
        Pylon_Sweep  double     %[double] Pylon sweep angle [º]
        SFC          double     %[double] Specific fuel consumption [lb/(shp·h)]
        Swet         double     % Wet surface of the engine (prop not included) [m^2]
        Thrust       double     %[double] Thrust of the engine [N].
        TotalPower   double     %[double] Total take-off power [kW]
        TotalThrust  double     %[double] Total thrust of the engines [N].
        TotalWeight  double     %[double] Total mass of the engines [kg].
        TSFC         double     %[double] Thrust specific fuel consumption [lb/(lbf·h)]
        TSFC_TO      double     %[double] Thrust specific fuel consumption at Take-Off [lb/(lbf·h)]
        Type         string     %[string] Type of engine (jet, propeller, turbofan,...)
        Weight       double     %[double] Mass of the engine [kg]. ONE ENGINE
        Width        double     % Maximum width of the nacelle [m]
    end
    
    
    properties (Dependent, GetObservable)

    end
    
    
    methods
        % Class constructor
        function obj = engines(varargin)
            if (nargin > 0) && (isa(varargin{1},'aircraft'))
                addlistener(obj,'TotalThrust','PostSet',@(src,evnt)aircraft.TtoMTOWModification(src,evnt,varargin{1}));
            end
        end
    end
        
    
end

