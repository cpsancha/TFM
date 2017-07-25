function [ R2 ] = getR2( x, y, fit )
%GETR2 Summary of this function goes here
%   Detailed explanation goes here
    
    polydata = polyval(fit,x);
    sstot = sum((y - mean(y)).^2);
    ssres = sum((y - polydata).^2);
    R2 = 1 - (ssres / sstot);

end

