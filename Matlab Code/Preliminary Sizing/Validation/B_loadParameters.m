%% Constants
g = 9.80665;

%% Conversion factors

CF.nm_to_m = 1852;
CF.mph_to_ms = 2.236936;
CF.sm_to_m   =  1.6093e+03;
CF.kts_to_ms = 0.514444;
CF.lbs_to_kg = convmass(1,'lbm','kg');
CF.hp_to_watts = 745.7 ;
CF.c_p_to_SI = CF.lbs_to_kg*(1/CF.hp_to_watts)*(1/3600)*g;

%% Table 2.2: Suggested Values For L/D, c_j, n_p and c_p for several mission phases:
switch ME.Type
    case 5
        
% 5.Business Jets
parameters.Cruise.L_D = [10, 12];
parameters.Cruise.c_j = [0.5, 0.9]/3600;  %lbs/lbs/hr to 1/s
parameters.Cruise.c_p = NaN;              %Bussiness JETS can't be propeller driven
parameters.Cruise.n_p = NaN;              %Bussiness JETS can't be propeller driven

parameters.Loiter.L_D = [12, 14];
parameters.Loiter.c_j = [0.4, 0.6]/3600;   %lbs/lbs/hr to 1/s
parameters.Loiter.c_p = NaN;               %Bussiness JETS can't be propeller driven
parameters.Loiter.n_p = NaN;               %Bussiness JETS can't be propeller driven

    case 11
%11. Flying boats, amphibious, float airplanes
parameters.Cruise.L_D = [10, 12];
parameters.Cruise.c_j = [0.5, 0.9]/3600;          %lbs/lbs/hr to 1/s
parameters.Cruise.c_p = [0.5, 0.7]*CF.c_p_to_SI;  %lbs/hp/hr to N/Watts/s
parameters.Cruise.n_p = 0.82;

parameters.Loiter.L_D = [13, 15];
parameters.Loiter.c_j = [0.4, 0.6]/3600;          %lbs/lbs/hr to 1/s
parameters.Loiter.c_p = [0.5, 0.7]*CF.c_p_to_SI;  %lbs/hp/hr to N/Watts/s
parameters.Loiter.n_p = 0.77;

    otherwise
        warning('Unexpected mission type.')
end

%% Table 2.1: Suggested Fuel-Fractions for Several Mission Phases

switch ME.Type
    case 5
% 5.Business Jets
parameters.fuelFraction(1).phase = 'W1/W_TO: Engine start, warm-up';
parameters.fuelFraction(1).value = 0.990;

parameters.fuelFraction(2).phase = 'W2/W1: Taxi';
parameters.fuelFraction(2).value = 0.995;

parameters.fuelFraction(3).phase = 'W3/W2: Take-off';
parameters.fuelFraction(3).value = 0.995;

parameters.fuelFraction(4).phase = 'W4/W3: Climb';
parameters.fuelFraction(4).value = 0.980;  %Can be calculated from Breguet's equation

parameters.fuelFraction(5).phase = 'W5/W4: Cruise';
parameters.fuelFraction(5).value = NaN;  %From Breguet's equation

parameters.fuelFraction(6).phase = 'W6/W5: Loiter';
parameters.fuelFraction(6).value = NaN;  %From Breguet's equation

parameters.fuelFraction(7).phase = 'W7/W6: Descent';
parameters.fuelFraction(7).value = 0.990; 

parameters.fuelFraction(8).phase = 'W8/W7: Landing, taxi and shut-down';
parameters.fuelFraction(8).value = 0.992; 

    case 11
%11. Flying boats, amphibious, float airplanes
parameters.fuelFraction(1).phase = 'W1/W_TO: Engine start, warm-up';
parameters.fuelFraction(1).value = 0.992;

parameters.fuelFraction(2).phase = 'W2/W1: Taxi';
parameters.fuelFraction(2).value = 0.996;

parameters.fuelFraction(3).phase = 'W3/W2: Take-off';
parameters.fuelFraction(3).value = 0.996;

parameters.fuelFraction(4).phase = 'W4/W3: Climb';
parameters.fuelFraction(4).value = 0.990;  %Can be calculated from Breguet's equation

parameters.fuelFraction(5).phase = 'W5/W4: Cruise';
parameters.fuelFraction(5).value = 0.863;  %From Breguet's equation

parameters.fuelFraction(6).phase = 'W6/W5: Loiter';
parameters.fuelFraction(6).value = 1;  %From Breguet's equation

parameters.fuelFraction(7).phase = 'W7/W6: Descent';
parameters.fuelFraction(7).value = 0.992; 

parameters.fuelFraction(8).phase = 'W8/W7: Landing, taxi and shut-down';
parameters.fuelFraction(8).value = 0.992; 

    otherwise
        warning('Unexpected mission type.')
end


%% Breguet´s equation to get fuel fraction during cruise and loiter:

% %Updating estructure 'parameters' with computed fuel fractions from Breguet's equation 
% parameters = getFuelFraction('W5/W4: Cruise',ME,parameters);
% parameters = getFuelFraction('W6/W5: Loiter',ME,parameters);

%% Table 2.15: Regresion constants A and B of equation 2.16:
switch ME.Type
    case 5
% 5.Business Jets
   parameters.A = 0.2678;
   parameters.B = 0.9979; 
    case 11
%11. Flying boats, amphibious, float airplanes
   parameters.A = 0.0966;
   parameters.B = 1.0298; 
end 








%% Auxiliar function definition

function [parameters] = getFuelFraction( phase, ME, parameters )
%DOC: Computes fuel fraction using Breguet's equation for range/endurance
%according to the mission phase , mission especification and propulsion paramaters.
switch phase
    case 'W5/W4: Cruise'
        %Cruise
        switch ME.Powerplant.Type
            case 'jet'
        cj  = parameters.Cruise.c_j(end);  %Worst case: maximo consumo
        L_D = parameters.Cruise.L_D(1);    %Worst case: minima eficiencia

        parameters.fuelFraction(5).value=1/exp(ME.Cruise.Range*cj/(ME.Cruise.Speed*L_D));
            case 'propeller'

        R   = ME.Cruise.Range;
        cp  = parameters.Cruise.c_p(end);  %Worst case: maximo consumo
        np  = parameters.Cruise.n_p(1);    %n_p is not an array, just in case
        L_D = parameters.Cruise.L_D(1);    %Worst case: minima eficiencia

        parameters.fuelFraction(5).value=1/exp(R*(cp/np)/L_D);
        end
        
    case 'W6/W5: Loiter'
        %Loiter
        switch ME.Powerplant.Type
            case 'jet'

        cj  = parameters.Loiter.c_j(end);  %Worst case: maximo consumo
        L_D = parameters.Loiter.L_D(1);    %Worst case: minima eficiencia

        parameters.fuelFraction(6).value=1/exp(ME.Loiter.E_ltr*cj/L_D);
            case 'propeller'

        cp  = parameters.Loiter.c_p(end);  %Worst case: maximo consumo
        np  = parameters.Loiter.n_p(1) ;   %n_p is not an array, just in case
        L_D = parameters.Loiter.L_D(1);    %Worst case: minima eficiencia

        parameters.fuelFraction(6).value=1/exp(ME.Loiter.E_ltr*ME.Loiter.V_ltr*(cp/np)/L_D);
        end
    otherwise
        warning('Wrong mission phase specifiec')
end

end



