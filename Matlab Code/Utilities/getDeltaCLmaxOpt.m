function deltaCLmax = getDeltaCLmaxOpt(x, Swf_S , cle_c, cl_alpha)
%This function obtains the total increment of CLmax in the wing using the
%following parameters:
% Swf_S: flap area versus wing area (ver figura)
% cf_c : chord flap versus chord
%delta_f : flap deflection [º]
% cle_c: chord leading edge devices
% cl_alpha: airfoil lift slope
cf_c= x(1);
delta_f=x(2);

load('alphaDeltaf.mat')
load('alphaDeltaf15.mat')
a15 = interp1(alphaDeltaf15(:,1),alphaDeltaf15(:,2),delta_f,'linear','extrap');
a40 = interp1(alphaDeltaf(:,1),alphaDeltaf(:,2),delta_f,'linear','extrap');
alpha_delta_f = interp1([0.15 , 0.4] , [a15 , a40], cf_c,'linear','extrap')

load('K_Clmax_Cl.mat')
K = interp1( K_Clmax_Cl(:,1),K_Clmax_Cl(:,2), cf_c,'linear','extrap')

deltaCLmax =( cle_c * cl_alpha * (1 + cf_c)* alpha_delta_f * (delta_f*pi/180) * K * Swf_S)^-1;


%Optimization input
% @(x)getDeltaCLmax(x, 0.5 , 1, 2*pi)
% % Result
% x = [0.3297 44.2597]

end