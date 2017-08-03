classdef wing < handle
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
        Airfoil         struct  %[struct] Data of the used Airfoil [-]
        alpha_zeroLift  double  %[double] Angle of attack of the full wing for zero lift [º]
        AspectRatio     double  %[double] Aspect Ratio (b/CMA) (b^2/Sw) [-]
        CL_alpha        double  %[double] Pendiente de la curva de sustentacion del perfil
        CLdesign        double  %[double] Cl of design of the airfoil
        CLmax           double  %[double] Maximum Cl of the wing at cruise
        CLmax_L         double  %[double] Maximum Cl of the wing at Landing
        CLmax_TO        double  %[double] Maximum Cl of the wing at Take-Off
        CMA             double  %[double] Mean Aerodynamic Chord [m]
        CMA_14          double  %[double] Longitudinal Position of the point 1/4 of the Mean Aerodynamic Chord
        CMA_b           double  %[double] Distancia del fuselaje a la CMA de forma horizontal (Distancia entre Croot-->CMA)
        CMA_LE          double  %[double] Longitudinal Position of the Mean Aerodynamic Chord Leading Edge
        CMG             double  %[double] Mean Geometric Chord [m]
        Dihedral        double  %[double] Dihedral of the wing [º]
        LongPos         double  %[double] Definida como el cociente entre la distancia longitudinal del punto un cuarto de la cuerda media aerodinámica al morro del avión y la longitud del fuselaje.
        MachDiv         double  %[double] Mach of divergence, for which drag increases a lot [-]
        RealSemiSpan    double  %[double] Distance from the fuselage (root chord) to the tip chord
        Root_LE         double  %[double] Longitudinal Position of the Root Leading Edge 
        RootChord       double  %[double] Chord at the root [m]
        RootWidth       double  %[double] t at the root [m]
        Snet            double  %[double] Total wing surface subtracting the portion contained in the fuselage [m^2]
        Sw              double  %[double] Total Wing Surface [m^2]
        Sweep_12        double  %[double] Sweep of the wing at the point 1/2 of CMA [º]
        Sweep_14        double  %[double] Sweep of the wing at the point 1/4 of CMA [º]
        Sweep_LE        double  %[double] Sweep of the wing at the leading edge [º]
        t_c             double  %[double] Thickness ratio. Maximum width/chord length [-]
        TaperRatio      double  %[double] Wing taper ratio (tipChord/rootChord) [-]        
        TaperRatio_i    double  %[double] Inborad Wing taper ratio (tipChord/rootChord) [-]           
        TaperRatio_o    double  %[double] Outboard Wing taper ratio (tipChord/rootChord) [-]
        TipChord        double  %[double] Chord ath the tip [m]
        TipSweep        double  %[double] Distancia longitudinal del borde de ataque de la raiz, al borde de ataque en la punta debido a la flecha
        TipTwist        double  %[double] Geometrical twist in the section located at the tip of the wing [º]
        Torsion         double  %[double] Ángulo de torsión de la sección en punta del ala, [RADIANES]
        WingLoading     double  %[double] MTOW/Sw [kg/m^2]
        WingSpan        double  %[double] Span of the wing [m]
        Incidence       double  % iw , Angle of incidence mesaured from fuselage floor [º]
    end
    
    methods
         % Class constructor
        function obj = wing(varargin)
            obj.Airfoil = struct('Name', '');
        end
        
    end
    
end

