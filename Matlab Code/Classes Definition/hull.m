classdef hull < handle
    %WING Summary of this class goes here
    %   Detailed explanation goes here
    % 
    % Mean Aerodynamic Chord (MAC - CMA) extra information:
    %   Rectangular wing: MAC = chordlength
    %   Deltawing: MAC = 2/3 inner cord length
    %   Trapezoid wing: MAC = 2/3 * innerchord * ((1+lambda+lambda^2)/(1+lambda)) --> (lambda = outerchord / innerchord)
    %   Elliptical wing: MAC = 8/3 * pi * innerchord
    %   MAC is used to have a dimensionless value for a position (like ca, cl...) for example: COG = 29% of MAC (this works for all wings, if MAC is calculated correctly)
    %
    
    properties
        Beam              double  %Beam [m]
        Beta              double  %Deadrise angle [º]
        Length            double  %Hull length [m]
        Length_Beam       double  %Hull to beam ratio
        Lf                double  %Forebody length [m]
        La                double  %Afterbody length [m]
        Lf_Beam           double  %forebody ratio
        Cv0               double  %Static load coeficcient
        Height            double  %
    end
    
    methods
         % Class constructor
        function obj = hull(varargin)
        end
        
    end
    
end

