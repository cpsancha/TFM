function [ fit, RSquared ] = polyfitR2( x, y, degree )
%POLYFITR2 Summary of this function goes here
%   Detailed explanation goes here

    fit = polyfit(x,y,degree);
    polydata = polyval(fit,x);
    sstot = sum((y - mean(y)).^2);
    ssres = sum((y - polydata).^2);
    RSquared = 1 - (ssres / sstot);
    
end

