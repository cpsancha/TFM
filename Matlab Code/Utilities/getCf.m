function [Cf] = getCf( Reynolds, xt_c)
x=[1e6,    1e7];
y=[0.0045, 0.003];
CF0=10.^polyval(polyfit(log10(x),log10(y),1),log10(Reynolds));

x=[1e6,    1e7];
y = [ 0.0013, 0.00041];
CF1=10.^polyval(polyfit(log10(x),log10(y),1),log10(Reynolds));

Cf = ones(1,length(Reynolds));
for ii=1:length(Reynolds)
    Cf(ii)=interp1([0,1],[CF0(ii),CF1(ii)],xt_c);
end
% loglog(linspace(1e6,1e7,length(Reynolds)),Cf)
end