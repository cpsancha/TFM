classdef weights < handle
    %WEIGHTS A MATLAB class to contain the elementary weights of an aircraft
    %   Detailed explanation goes here
    
    properties (SetObservable) 
        MTOW     double    %[double]  Maximum Take Off Weight [kg]
        MRW      double    %[double]  Maximum Ramp Weight [kg]
        EW       double    %[double]  Empty Weight, full equiped [kg]
        OEW      double    %[double]  Operational Empty Weight [kg]
        BOW      double    %[double]  Basic Operational Weight [kg]
        MPL      double    %[double]  Maximum Payload Weight [kg]
        MFW      double    %[double]  Maximum Fuel Weight [kg]
        TUL      double    %[double]  Total Useful Load [kg]
        MZFW     double    %[double]  Maximum Zero Fuel Weight [kg]
        MLW      double    %[double]  Maximum Landing Weight [kg]
        Tto_MTOW double    %[double]  Ratio of total take-off thrust and MTOW [N/kg]
    end
    
    
    properties (Dependent, GetObservable) 
        OEW_MTOW double    %[double]  Ratio of the OEW and the MTOW [-]
        MPL_MTOW double    %[double]  Ratio of the MPL and the MTOW [-]
        MFW_MTOW double    %[double]  Ratio of the MFW and the MTOW [-]
    end
   
    
    methods
        % Class constructor
        function obj = weights(varargin)
            if (nargin > 0) && (isa(varargin{1},'aircraft'))
                addlistener(obj,'MTOW','PostSet',@(src,evnt)aircraft.TtoMTOWModification(src,evnt,varargin{1}));
            end
        end
                    
        % OEW_MTOW getter function
        function OEW_MTOW = get.OEW_MTOW(obj)
            OEW_MTOW = obj.OEW / obj.MTOW;
        end
        
        % MPL_MTOW getter function
        function MPL_MTOW = get.MPL_MTOW(obj)
            MPL_MTOW = obj.MPL / obj.MTOW;
        end
        
        % MFW_MTOW getter function
        function MFW_MTOW = get.MFW_MTOW(obj)
            MFW_MTOW = obj.MFW / obj.MTOW;
        end  
    end
    
    
end

