%% Mission especifications (ME)

%Choose one:
ME.Type = 5 ; % Business Jets
ME.Type = 11; % Flying boats, amphibious, float airplanes

%% Payload
Passengers = 6;

passengerWeight = convmass(175,'lbm','kg');
passengerBaggageWeight   = convmass(40,'lbm','kg');

%Mission payload weight in kg:
ME.Payload =  Passengers*(passengerWeight+passengerBaggageWeight);   

ME.Payload = 1250;

%% Crew
Crew = 2; %Check FAR 91.215 for minimun crew members

crewWeight = convmass(175,'lbm','kg');
crewBaggageWeight   = convmass(30,'lbm','kg');

%Mission crew weight in kg:
ME.Crew = Crew*(crewWeight+crewBaggageWeight);

%% Take-off
ME.TakeOff = NaN;

%% Climb
ME.Climb.E_cl = 1;    %Time to climb   in seconds
ME.Climb.V_cl = NaN;  %Climb speed in m/s (horizontal)

%% Cruise

ME.Cruise.Range = 2659472;     % in m

ME.Cruise.Altitude = nan;      % in m

ME.Cruise.Speed = 243.3320;    % m/s

%% Loiter
ME.Loiter.E_ltr = 60*30;             %Loiter time in seconds
ME.Loiter.V_ltr = ME.Cruise.Speed;  %Loiter speed in m/s 



%% Landing
ME.Landing = NaN;

%% Powerplant

%Choose one:
ME.Powerplant.Type = 'propeller';
ME.Powerplant.Type = 'jet';      

%%

ME.Pressurization = NaN;

ME.Mission_Profile = NaN;



