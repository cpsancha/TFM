function [ A,B ] = getWTORegression( )
%GETWTOGUESS Summary of this function goes here
%   Detailed explanation goes here
    type = 5;
    SimilarPlanes = importSimilarPlanes(type);
    for i=1:length(SimilarPlanes)
        if ~isempty(SimilarPlanes{i}.Weight.OEW)
            MTOW(i) = SimilarPlanes{i}.Weight.MTOW;
            OEW(i)  = SimilarPlanes{i}.Weight.OEW;
        else
            MTOW(i) = NaN;
            OEW(i)  = NaN;
        end
    end
    
end

