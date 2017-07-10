function [ fieldSP ] = loadFields( SP, fieldString )
%LOADFIELDS Summary of this function goes here
%   Detailed explanation goes here

for i=1:length(SP)
    str = strcat('SP{',num2str(i),'}.',fieldString);
    if ~isempty(eval(str))
        field(i) = eval(str); %#ok<AGROW>
    else
        switch class(eval(str))
            case 'double'
                field(i) = NaN; %#ok<AGROW>
            case {'char','string'}
                field(i) = ''; %#ok<AGROW>
            otherwise
                field(i) = NaN; %#ok<AGROW>
        end
    end
end
fieldSP = field;

end

