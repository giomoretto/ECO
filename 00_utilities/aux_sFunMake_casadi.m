function aux_sFunMake_casadi( functionList, Options )
% aux_sFunMake_casadi( functionList, Options )
%
%   Author1     : Eugen Nuss (e.nuss@irt.rwth-aachen.de)
% 
%   Last change : Summer 2019
%
%   Description : Generates and compiles CasADi functions
%
%   Parameters  : functionList -> Cell array containing CasADi function 
%                                 objects or single CasADi function object
% 
%                 Options -> Struct containing options (optional)
% 
%   Return      : -
% 
%-------------------------------------------------------------------------%


%-------------------------------------------------------------------------%
%                           Set options                                   %
%-------------------------------------------------------------------------%
% Set standard options
OptionsStd.common_sourcefile = '';
OptionsStd.output_is_sparse = false;
OptionsStd.output_is_matrix = false;
OptionsStd.build_directory = 'build';
OptionsStd.generate_rtwmakecfg = 0;

% Get user defined options
if nargin == 2
    Options = aux_getOptions( OptionsStd, {'Options', Options} );
else
    Options = OptionsStd; 
end

if ~isempty( Options.common_sourcefile )
    % Set include header to defined common sourcefile
    Options.include_header = {Options.common_sourcefile};
    sourceFile             = Options.common_sourcefile;
    
elseif iscell(functionList) && isempty( Options.common_sourcefile)
    % If no common filename passed but several functions are to be created
    sourceFile             = 'common_sourcefile';
    Options.include_header = {sourceFile};
    
else
    sourceFile             = functionList.name;
    
end

% Check if functionList contains cell array of function objects
if ~iscell(functionList)
    buff = functionList;
    functionList = cell(1);
    functionList{1} = buff;
end

% Set options for CasADi c-code generation of CasADi function
OptionsCg.casadi_real      = 'real_T';
OptionsCg.casadi_int       = 'int_T';
OptionsCg.with_header      = true;

%-------------------------------------------------------------------------%




% Switch to temporary folder for code generation
if exist( [pwd '/' Options.build_directory], 'dir' )~=7
    mkdir( Options.build_directory )
end
currentDir = cd;
cd(Options.build_directory)




%-------------------------------------------------------------------------%
%                           Generate c-code                               %
%-------------------------------------------------------------------------%
% Initialize CasADi codegenerator object
import casadi.CodeGenerator

cg = CodeGenerator( sourceFile, OptionsCg);

% Add CasADi-functions to be generated
for ii = 1:length(functionList)
    cg.add( functionList{ii} );
    aux_sFunGen_casadi( functionList{ii}.name, Options); % functionname in c
    fprintf('Codegen complete: created sfun_%s.c \n', functionList{ii}.name);
end

% Include matlab header for use simulink utilities
cg.add_include('simstruc.h');

% Generate c-code
cg.generate();
fprintf('Codegen complete: created %s.c \n', sourceFile);
fprintf('Codegen complete: created %s.h \n\n', sourceFile);

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
%                          Compile s-functions                            %
%-------------------------------------------------------------------------%
% Generate matlab s-function for specified CasADi-functions
for ii = 1:length(functionList)
    cmd = ['mex sfun_' functionList{ii}.name '.c ' sourceFile '.c'];
    
    % Compile s-function
    eval(cmd);
    fprintf('Compilation complete: created sfun_%s.mexw64 \n', functionList{ii}.name);
end
disp(' ');

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
%                          Generate rtwmakecfg                            %
%-------------------------------------------------------------------------%
% Codegenerate customized rtwmakecfg
if Options.generate_rtwmakecfg
	fid = fopen('rtwmakecfg.m', 'wt');
	fprintf(fid, 'function MakeInfo = rtwmakecfg()\n');
	fprintf(fid, '%% MakeInfo = rtwmakecfg()\n');
	fprintf(fid, '%%\n');
	fprintf(fid, '%%   Author1     : Eugen Nuss (e.nuss@irt.rwth-aachen.de)\n');
	fprintf(fid, '%%\n');
	fprintf(fid, '%%   Last change : Spring 2019\n');
	fprintf(fid, '%%\n');
	fprintf(fid, '%%   Description : Loads user defined make options for simulink code\n');
	fprintf(fid, '%%                 generation. Needs to be located at the same \n');
	fprintf(fid, '%%                 place as s-function .mex-File.\n');
	fprintf(fid, '%%\n');
	fprintf(fid, '%%   Parameters  : -\n');
	fprintf(fid, '%%\n');
	fprintf(fid, '%%   Return      : MakeInfo -> Struct containing configure options\n');
	fprintf(fid, '%%\n');
	fprintf(fid, '%%-------------------------------------------------------------------------%%\n');
	fprintf(fid, '%% Define additional include directories\n');
	fprintf(fid, 'MakeInfo.includePath = { fullfile(pwd, ''/'')};\n');
	fprintf(fid, '\n');
	fprintf(fid, '%% Define additional source files\n');
	fprintf(fid, ['MakeInfo.sources     = { ''' Options.common_sourcefile '.c''};\n']);
	fprintf(fid, '\n');
	fprintf(fid, '%% Issue a message\n');
	fprintf(fid, 'seperatorLine = char(ones(1,120) * ''~'');\n');
	fprintf(fid, '\n');
	fprintf(fid, 'fprintf(''%%s\\n\\n'', seperatorLine);\n');
	fprintf(fid, 'fprintf(''Running rtwmakecfg from folder: %%s\\n\\n'',pwd);\n');
	fprintf(fid, ['fprintf(''Adding legacy code to make process for ' sourceFile '\\n\\n'')\n']);
	fprintf(fid, 'fprintf('' - additional include directories:\\n'');\n');
	fprintf(fid, 'fprintf(''       %%s\\n\\n'', MakeInfo.includePath{:});\n');
	fprintf(fid, 'fprintf('' - additional source files:\\n'');\n');
	fprintf(fid, 'fprintf(''       %%s\\n\\n'', MakeInfo.sources{:});\n');
	fprintf(fid, 'fprintf(''%%s\\n'', seperatorLine);\n');
	fclose(fid);
end

%-------------------------------------------------------------------------%


% Go back to project root folder
cd(currentDir);