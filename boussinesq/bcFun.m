function res = bcFun(r, allVars, Ndom, whichEnd)
    % allVars is a cell array: {A_{-N}...A_{N}, B_{-N}...B_{N}, C_{-N}...C_{N}}
    % whichEnd: 'left' or 'right'
    nModes = 2*Ndom + 1;
    Ac = allVars(1:nModes);
    Bc = allVars(nModes+1:2*nModes);
    Cc = allVars(2*nModes+1:3*nModes);

    res = [];
    for k = -Ndom:Ndom
        res = [res; Ac{k+Ndom+1}]; %#ok<AGROW>
    end
    for k = -Ndom:Ndom
        target = 0;
        if strcmp(whichEnd,'left') && k == 0
            target = 1;
        end
        res = [res; Bc{k+Ndom+1} - target]; %#ok<AGROW>
    end
    for k = -Ndom:Ndom
        res = [res; Cc{k+Ndom+1}]; %#ok<AGROW>
    end
end