
%% ===================== Setup =====================
Ndom   = 3;
nModes = 2*Ndom + 1;
dom    = [1, 5];
lambda = 0.5;
nTotal = 3*nModes;

varNames = cell(1, nTotal);
for k = -Ndom:Ndom, varNames{k+Ndom+1}            = sprintf('A%d', k+Ndom); end
for k = -Ndom:Ndom, varNames{nModes + k+Ndom+1}   = sprintf('B%d', k+Ndom); end
for k = -Ndom:Ndom, varNames{2*nModes + k+Ndom+1} = sprintf('C%d', k+Ndom); end
argList = strjoin(varNames, ',');

%% ===================== Write a real operator file to disk =====================
opFile = fullfile(pwd, 'genOp.m');
fid = fopen(opFile, 'w');
fprintf(fid, 'function out = genOp(r,%s)\n', argList);
fprintf(fid, '    out = buildSystem(r, {%s}, %d, %.15g, [%.15g, %.15g]);\n', ...
        argList, Ndom, lambda, dom(1), dom(2));
fprintf(fid, 'end\n');
fclose(fid);
rehash path

opFcn = @genOp;
N_op  = chebop(opFcn, dom);


N_op.lbc = @genLBC;
N_op.rbc = @genRBC;

%% ===================== Solve =====================
sol = N_op \ 0;

