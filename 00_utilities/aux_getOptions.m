function [options] = aux_getOptions(options, inputArgs)
% aux_sFunGen_casadi( functionName, Options )
%
%   Author1     : Biafra Ahanonu (http://bahanonu.com/syscarut/articles/133/ 2014.02.17 [22:21:49])
%   Author2     : Eugen Nuss (e.nuss@irt.rwth-aachen.de)
% 
%   Date        : Winter 2018
%
%   Description : Gets default options for a function, replaces with 
%                 inputArgs inputs if they are present
%
%   Parameters  : Options -> Structure with options as fieldnames
% 
%                 inputArgs -> Varargin containing name-value pairs passed 
%                              from parent function
% 
%   Return      : -
% 
%-------------------------------------------------------------------------%
% List of valid options to accept, simple way to deal with illegal user input
validOptions = fieldnames(options);

% Loop over each input name-value pair, check whether name is valid and overwrite fieldname in options structure.
for ii = 1:2:length(inputArgs)
    val = inputArgs{ii};
    if ischar(val)
        % allow input of an options structure that overwrites existing fieldnames with its own, for increased flexibility
        if strcmp('Options',val)
            inputOptions = inputArgs{ii+1};
            [options] = mirrorRightStruct(inputOptions,options);
        elseif ~isempty(strcmp(val,validOptions))
            options.(val) = inputArgs{ii+1};
        end
    else
        continue;
    end
end

function [pullStruct] = mirrorRightStruct(pushStruct,pullStruct)
    % overwrites fields in pullStruct with those in pushStruct, other pullStruct fields rename intact
    % more generally, copies fields in pushStruct into pullStruct, if there is an overlap in field names, pushStruct overwrites.
    pushNames = fieldnames(pushStruct);
    for name = 1:length(pushNames)
        iName = pushNames{name};
        pullStruct.(iName) = pushStruct.(iName);
    end