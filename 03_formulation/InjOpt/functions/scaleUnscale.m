
function [parsedVars,scaleVec,offsVec] = scaleUnscale(unParsedVars,...
                                     VarNames,optOpts,scaleChar,isgrad)
    if isfield(optOpts,'ScaleOffsVars') && ...
            ~isempty(optOpts.ScaleOffsVars)
        % Get Opts
        UnScaleOpts.Vars    = optOpts.ScaleOffsVars;
        UnScaleOpts.Offs    = optOpts.Offs;
        UnScaleOpts.Scale   = optOpts.Scale;
        % Scale Variables
        scaleVec = ones(size(unParsedVars));
        offsVec  = zeros(size(unParsedVars));
        for i = 1:numel(VarNames)
            indUnScale = ismember(UnScaleOpts.Vars,VarNames{i});
            if any(indUnScale)
                scaleVec(i) = UnScaleOpts.Scale(indUnScale);
                offsVec(i)  = UnScaleOpts.Offs(indUnScale);
            end
        end
        if strcmp(scaleChar,'unscale')
            if isgrad
                parsedVars = unParsedVars.*scaleVec;
                offsVec    = zeros(size(unParsedVars));
            else
                parsedVars = unParsedVars.*scaleVec + offsVec;
            end
        else
            if isgrad
                parsedVars = unParsedVars./scaleVec;
                offsVec    = zeros(size(unParsedVars));
            else
                parsedVars = (unParsedVars - offsVec)./scaleVec;
            end
        end
    else
        parsedVars = unParsedVars;
        scaleVec   = ones(size(unParsedVars));
        offsVec    = zeros(size(unParsedVars));
    end
end