classdef engines < handle
    %ENGINES A MATLAB class to contain the elementary data of an aircraft
    %engines.
    %   Detailed explanation goes here
    
    properties (SetObservable)
        Manufacturer string     %[string] Engine's manufacturer.
        Model        string     %[string] Engine's model.
        Type         string     %[string] Type of engine (jet, propeller, turbofan,...)
        Number       double     %[double] Number of engines in the plane.
        Position     double     %[double] Position of the CoG (x y z) of the engine [m].
        PositionStr  string     %[string] Verbose Position of the engine/engines
        Weight       double     %[double] Mass of the engine [kg].
        Thrust       double     %[double] Thrust of the engine [kN].
        TSFC         double     %[double] Thrust specific fuel consumption [lb/(lbf·h)]
        SFC          double     %[double] Specific fuel consumption [lb/(shp·h)]
        etaPropeller double     %[double] Propeller efficiency [-]
        Diameter     double     %[double] Maximum diameter of the engine [m]
        Length       double     %[double] Maximum length of the engine [m]
    end
    
    
    properties (Dependent, GetObservable)
        TotalThrust double      %[double] Total thrust of the engine [kN].
        TotalWeight double      %[double] Total mass of the engine [kg].
    end
    
    
    methods
        % Class constructor
        function obj = engines(varargin)
            if (nargin > 0) && (isa(varargin{1},'aircraft'))
                addlistener(obj,'Number','PostSet',@(src,evnt)aircraft.TtoMTOWModification(src,evnt,varargin{1}));
                addlistener(obj,'Thrust','PostSet',@(src,evnt)aircraft.TtoMTOWModification(src,evnt,varargin{1}));
            end
        end
        
        % TotalThrust getter function
        function TotalThrust = get.TotalThrust(obj)
            TotalThrust = obj.Number * obj.Thrust;
        end
        
        %  TotalWeight getter function
        function TotalWeight = get.TotalWeight(obj)
            TotalWeight = obj.Number * obj.Weight;
        end
    end
        
    
end

