function aux_sFunGen_casadi( functionName, Options )
% aux_sFunGen_casadi( functionName, Options )
%
%   Author1     : Eugen Nuss (e.nuss@irt.rwth-aachen.de)
% 
%   Date        : Spring 2019
%
%   Description : Generates a matlab s-function specific for the CasADi
%                 codegeneration API. Suitable for both sparse and dense
%                 outputs as well as matrix outputs.
%
%   Parameters  : functionName -> String containing the C-CODE NAME of the
%                                 CasADi function to be generated
%                                 ( casadiFun = Function('C-CODE NAME',...
%                                   or casadiFun.name )
% 
%                 Options -> Struct containing options (optional)
% 
%   Return      : -
% 
%-------------------------------------------------------------------------%
% Set standard options
OptionsStd.output_is_sparse     = false;
OptionsStd.output_is_matrix     = true;
OptionsStd.input_is_matrix      = false;
OptionsStd.include_header       = {functionName};
OptionsStd.enable_debug_dspace  = false;
OptionsStd.enable_tictoc_dspace = false;
OptionsStd.options_s_function   = {'SS_OPTION_USE_TLC_WITH_ACCELERATOR', ...
                                  'SS_OPTION_CAN_BE_CALLED_CONDITIONALLY', ...
                                  'SS_OPTION_EXCEPTION_FREE_CODE', ...
                                  'SS_OPTION_WORKS_WITH_CODE_REUSE', ...
                                  'SS_OPTION_SFUNCTION_INLINED_FOR_RTW', ...
                                  'SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME'};

% Get user defined options
if nargin == 2
    Options = aux_getOptions(OptionsStd, {'Options', Options} );
else
    Options = OptionsStd;
end




fid = fopen(['sfun_' functionName '.c'], 'wt');


fprintf(fid, '// Must specify the S_FUNCTION_NAME as the name of the S-function\n');
fprintf(fid, ['#define S_FUNCTION_NAME  sfun_' functionName '\n']);
fprintf(fid, '#define S_FUNCTION_LEVEL 2\n\n');

fprintf(fid, '// Need to include simstruc.h for the definition of the SimStruct and\n');
fprintf(fid, '// its associated macro definitions\n');
fprintf(fid, '#ifndef __SIMSTRUC__\n');
fprintf(fid, '#include "simstruc.h"\n');
fprintf(fid, '#endif\n\n');

fprintf(fid, '// Specific header file(s) required by the legacy code function\n');
for ii = 1:length(Options.include_header)
    fprintf(fid, ['#include "' Options.include_header{ii} '.h"\n\n']);
end
if Options.enable_debug_dspace || Options.enable_tictoc_dspace
    fprintf(fid, '// Include header for dSPACE msg_info_printf, RTLIB_TIC_START and RTLIB_TIC_READ command on MABX\n');
    fprintf(fid, '#ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '	#include <brtenv.h> \n');
    fprintf(fid, '#endif\n\n\n\n\n');
end



%-------------------------------------------------------------------------%
%   Initialize function                                                   %
%-------------------------------------------------------------------------%
fprintf(fid, '/* Function: mdlInitializeSizes ===========================================\n');
fprintf(fid, ' * Abstract:\n');
fprintf(fid, ' *   The sizes information is used by Simulink to determine the S-function\n');
fprintf(fid, ' *   blocks characteristics (number of inputs, outputs, states, etc.).\n');
fprintf(fid, ' */\n');
fprintf(fid, 'static void mdlInitializeSizes(SimStruct *S)\n');
fprintf(fid, '{\n');
fprintf(fid, '	    //\n');
fprintf(fid, '	  ////\n');
fprintf(fid, '    //  // \n');
fprintf(fid, '	    // \n');
fprintf(fid, '	    // \n');
fprintf(fid, '	    //\n');
fprintf(fid, '	////////////////////////////////////////////////////////\n');
fprintf(fid, '    /* Set number of simulink s-function block parameters */\n');
fprintf(fid, '	/* Double click on simulink block 					  */\n');
fprintf(fid, '	////////////////////////////////////////////////////////\n');
fprintf(fid, '    ssSetNumSFcnParams(S, 0);\n');
fprintf(fid, '    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {\n');
fprintf(fid, '        return; // Parameter mismatch will be reported by Simulink\n');
fprintf(fid, '    }\n');
fprintf(fid, '	//--------------------------//\n\n\n');

	
	
fprintf(fid, '	 /////\n');
fprintf(fid, '	//   //\n');
fprintf(fid, '	   //\n');
fprintf(fid, '	 //\n');
fprintf(fid, '	///////\n');
fprintf(fid, '	////////////////////////////////////////\n');
fprintf(fid, '    /* Read in CasADi function dimensions */\n');
fprintf(fid, '	////////////////////////////////////////\n');
fprintf(fid, ['    int_T n_in  = ' functionName '_n_in();\n']);
fprintf(fid, ['    int_T n_out = ' functionName '_n_out();\n']);
fprintf(fid, '    int_T sz_arg, sz_res, sz_iw, sz_w;\n');
fprintf(fid, ['    ' functionName '_work(&sz_arg, &sz_res, &sz_iw, &sz_w);\n']);
fprintf(fid, '	//--------------------------//\n\n\n');
	
	
	
fprintf(fid, '	  //\n');
fprintf(fid, '	//  //\n');
fprintf(fid, '	   //\n');
fprintf(fid, '	    //\n');
fprintf(fid, '	// // \n');
fprintf(fid, '	//////////////////////////////////////////////\n');
fprintf(fid, '    /* Set up simulink work vectors             */\n');
fprintf(fid, '	//////////////////////////////////////////////\n');
fprintf(fid, '    ssSetNumRWork(S, sz_w);\n');
fprintf(fid, '    ssSetNumIWork(S, sz_iw);\n');
if Options.output_is_sparse
    fprintf(fid, '    ssSetNumPWork(S, sz_arg+sz_res+n_out);\n');
    fprintf(fid, '	// Set the number of work vectors\n');
    fprintf(fid, '    if (!ssSetNumDWork(S, n_out)) return;\n');
else
    fprintf(fid, '    ssSetNumPWork(S, sz_arg+sz_res);\n');
end
fprintf(fid, '    ssSetNumNonsampledZCs(S, 0);\n');
fprintf(fid, '	//--------------------------//\n\n\n');
		
	
	
fprintf(fid, '      // //\n');
fprintf(fid, '	 //  //\n');
fprintf(fid, '	//////\n');
fprintf(fid, '		//\n');
fprintf(fid, '		//\n');
fprintf(fid, ' 	////////////////////////////////////////\n');
fprintf(fid, '    /* Set up simulink input/output ports */\n');
fprintf(fid, '	////////////////////////////////////////\n');
fprintf(fid, '    int_T ii;\n');
fprintf(fid, '    const int_T* sp;\n');
fprintf(fid, '	// Inputs\n');
fprintf(fid, '    if (!ssSetNumInputPorts(S, n_in)) return;\n');
fprintf(fid, '    for (ii=0; ii<n_in; ++ii) {\n');
fprintf(fid, ['      sp = ' functionName '_sparsity_in(ii);\n']);
if Options.input_is_matrix
    fprintf(fid, '      // Dense matrix inputs assumed here\n');
    fprintf(fid, '      ssSetInputPortMatrixDimensions(S, ii, sp[0], sp[1]);\n');
else
    fprintf(fid, '      // Dense vector inputs assumed here\n');
    fprintf(fid, '      ssSetInputPortWidth(S, ii, sp[0]);\n');
end
fprintf(fid, '      ssSetInputPortDirectFeedThrough(S, ii, 1);\n');
fprintf(fid, '    }\n\n');
	
fprintf(fid, '	// Outputs	\n');
if Options.enable_tictoc_dspace
    fprintf(fid, '    if (!ssSetNumOutputPorts(S, n_out+1)) return;\n');
else
    fprintf(fid, '    if (!ssSetNumOutputPorts(S, n_out)) return;\n');
end
fprintf(fid, '    for (ii=0; ii<n_out; ++ii) {\n');
fprintf(fid, ['      sp = ' functionName '_sparsity_out(ii);\n']);
if Options.output_is_matrix
    fprintf(fid, '      ssSetOutputPortMatrixDimensions(S, ii, sp[0], sp[1]);\n');
else
    fprintf(fid, '      ssSetOutputPortWidth(S, ii, sp[0]);\n');
end
if Options.output_is_sparse
    fprintf(fid, '	  ssSetDWorkWidth(S, ii, sp[2+sp[1]]);\n');
end
fprintf(fid, '    }\n');
if Options.enable_tictoc_dspace
    fprintf(fid, '    // Add output for turnaround time\n');
    fprintf(fid, '    ssSetOutputPortWidth(S, n_out, 1);\n\n');
end
fprintf(fid, '	ssSetOutputPortOutputExprInRTW(S, 0, 0);\n');
fprintf(fid, '	//--------------------------//\n\n\n');

	
    
fprintf(fid, '	//////\n');
fprintf(fid, '	//\n');
fprintf(fid, '	/////\n');
fprintf(fid, '		//\n');
fprintf(fid, '	/////\n');
fprintf(fid, '	//////////////////////////////////////////////\n');
fprintf(fid, '    /* Not yet characterized 					*/\n');
fprintf(fid, '	//////////////////////////////////////////////\n');
fprintf(fid, '    /* This S-function can be used in referenced model simulating in normal mode */\n');
fprintf(fid, '    ssSetModelReferenceNormalModeSupport(S, MDL_START_AND_MDL_PROCESS_PARAMS_OK);\n\n');

fprintf(fid, '    /* Set the number of sample time */\n');
fprintf(fid, '    ssSetNumSampleTimes(S, 1);\n\n');

fprintf(fid, '    /* Set the compliance with the SimState feature */\n');
fprintf(fid, '    ssSetSimStateCompliance(S, USE_DEFAULT_SIM_STATE);\n\n');


fprintf(fid, '    /**\n');
fprintf(fid, '     * All options have the form SS_OPTION_<name> and are documented in\n');
fprintf(fid, '     * matlabroot/simulink/include/simstruc.h. The options should be\n');
fprintf(fid, '     * bitwise ord together as in\n');
fprintf(fid, '     *    ssSetOptions(S, (SS_OPTION_name1 | SS_OPTION_name2))\n');
fprintf(fid, '     */\n');
fprintf(fid, '    ssSetOptions(S,\n');
for ii = 1:length(Options.options_s_function)-1
    fprintf(fid, ['        ' Options.options_s_function{ii} ' |\n']);
end
fprintf(fid, ['        ' Options.options_s_function{end} '\n']);
fprintf(fid, '    );\n');
fprintf(fid, '	//--------------------------//\n\n\n');

	
if OptionsStd.enable_debug_dspace
    fprintf(fid, '	   //\n');
    fprintf(fid, '	 //  //\n');
    fprintf(fid, '	// \n');
    fprintf(fid, '	//////\n');
    fprintf(fid, '	//   //\n');
    fprintf(fid, '	//////\n');
    fprintf(fid, '	//////////////////////////////////////////////\n');
    fprintf(fid, '    /* Debugging stuff							*/\n');
    fprintf(fid, '	//////////////////////////////////////////////\n');
    fprintf(fid, '	/*\n');
    fprintf(fid, '	#ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '		msg_info_set(1, 1, "Test");\n');
    fprintf(fid, '		msg_info_printf(1, 2, "sz_arg %i", sz_arg);\n');
    fprintf(fid, '		msg_info_printf(1, 3, "sz_res %i", sz_res);\n');
    fprintf(fid, '		msg_info_printf(1, 4, "iw %i", sz_iw);\n');
    fprintf(fid, '		msg_info_printf(1, 5, "w %i", sz_w);\n');
    fprintf(fid, '	#endif\n');
    fprintf(fid, '	*/\n');
    fprintf(fid, '	//--------------------------//\n');
end
fprintf(fid, '}\n\n\n');

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
%   Initialize sample times function                                      %
%-------------------------------------------------------------------------%
fprintf(fid, '/* Function: mdlInitializeSampleTimes =====================================\n');
fprintf(fid, ' * Abstract:\n');
fprintf(fid, ' *   This function is used to specify the sample time(s) for your\n');
fprintf(fid, ' *   S-function. You must register the same number of sample times as\n');
fprintf(fid, ' *   specified in ssSetNumSampleTimes.\n');
fprintf(fid, ' */\n');
fprintf(fid, 'static void mdlInitializeSampleTimes(SimStruct *S)\n');
fprintf(fid, '{\n');
fprintf(fid, '    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);\n');
fprintf(fid, '    ssSetOffsetTime(S, 0, FIXED_IN_MINOR_STEP_OFFSET);\n');

fprintf(fid, '    #if defined(ssSetModelReferenceSampleTimeDefaultInheritance)\n');
fprintf(fid, '    ssSetModelReferenceSampleTimeDefaultInheritance(S);\n');
fprintf(fid, '    #endif\n');
fprintf(fid, '}\n\n\n');

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
%   Model output function                                                 %
%-------------------------------------------------------------------------%
fprintf(fid, '/* Function: mdlOutputs ===================================================\n');
fprintf(fid, ' * Abstract:\n');
fprintf(fid, ' *   In this function, you compute the outputs of your S-function\n');
fprintf(fid, ' *   block. Generally outputs are placed in the output vector(s),\n');
fprintf(fid, ' *   ssGetOutputPortSignal.\n');
fprintf(fid, ' */\n');
fprintf(fid, 'static void mdlOutputs(SimStruct *S, int_T tid)\n');
fprintf(fid, '{\n');

if Options.enable_tictoc_dspace
    fprintf(fid, '    // Start measuring turnaround time\n');
    fprintf(fid, '    #ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '        RTLIB_TIC_START();\n');
    fprintf(fid, '    #endif\n\n\n');
end	

fprintf(fid, '		//\n');
fprintf(fid, '	  ////\n');
fprintf(fid, '    //  // \n');
fprintf(fid, '	    // \n');
fprintf(fid, '	    // \n');
fprintf(fid, '	    //\n');
fprintf(fid, '	/////////////////////////////////\n');
fprintf(fid, '    /* Read in function dimensions */\n');
fprintf(fid, '	/////////////////////////////////\n');
fprintf(fid, ['    int_T n_in  = ' functionName '_n_in();\n']);
fprintf(fid, ['    int_T n_out = ' functionName '_n_out();\n']);
fprintf(fid, '    int_T sz_arg, sz_res, sz_iw, sz_w;\n');
fprintf(fid, ['    ' functionName '_work(&sz_arg, &sz_res, &sz_iw, &sz_w);\n']);
fprintf(fid, '	//--------------------------//\n\n\n');



	
fprintf(fid, '	 /////\n');
fprintf(fid, '	//   //\n');
fprintf(fid, '	   //\n');
fprintf(fid, '	 //\n');
fprintf(fid, '	///////\n');
fprintf(fid, '	//////////////////////////////\n');
fprintf(fid, '	/* Set up work vectors      */\n');
fprintf(fid, '	//////////////////////////////\n');
fprintf(fid, '	void** p = ssGetPWork(S);\n');
fprintf(fid, '    const real_T** arg = (const real_T**) p;\n');
fprintf(fid, '	p += sz_arg;\n');
if Options.output_is_sparse
    fprintf(fid, '    real_T** res = (real_T**) p;\n');
    fprintf(fid, '	p += sz_res;\n');
    fprintf(fid, '	real_T** y = (real_T**) p;\n');
else
    fprintf(fid, '    real_T** y = (real_T**) p;\n');
end
fprintf(fid, '    real_T* w = ssGetRWork(S);\n');
fprintf(fid, '    int_T* iw = ssGetIWork(S);\n');
fprintf(fid, '	//--------------------------//\n\n\n');
	
	
	
fprintf(fid, '	  //\n');
fprintf(fid, '	//  //\n');
fprintf(fid, '	   //\n');
fprintf(fid, '	    //\n');
fprintf(fid, '	// // \n');
fprintf(fid, '	////////////////////////////////////////\n');
fprintf(fid, '	/* Point to input and output buffers  */\n');
fprintf(fid, '	////////////////////////////////////////\n');
fprintf(fid, '    int_T ii;\n');
	
fprintf(fid, '	// Inputs\n');
fprintf(fid, '    for (ii=0; ii<n_in;++ii) {\n');
fprintf(fid, '      arg[ii] = *ssGetInputPortRealSignalPtrs( S, ii );\n\n');


if OptionsStd.enable_debug_dspace 
    fprintf(fid, '	  // Debugging\n');
    fprintf(fid, '	  /*\n');
    fprintf(fid, '	  #ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '		msg_info_printf(1, 1, "Address of start u%i array arg[%i] %p", ii, ii, arg[ii]);\n');
    fprintf(fid, '		msg_info_printf(1, 1, "Address of start u%i array arg[%i] %p", ii, ii, *(arg+ii)); \n');
    fprintf(fid, '		msg_info_printf(1, 1, "Content of first element u%i array *arg[%i] %f", ii, ii, *arg[ii]);\n');
    fprintf(fid, '		msg_info_printf(1, 1, "Address of second element of u%i array arg[%i]+1 %p", ii, ii, arg[ii]+1);\n');
    fprintf(fid, '		msg_info_printf(1, 1, "Content of second element of u%i array*(arg[%i]+1) %f", ii, ii, *(arg[ii]+1));\n');
    fprintf(fid, '		msg_info_printf(1, 1, "Content of second element of u%i array arg[%i][1] %f", ii, ii, arg[ii][1]);\n');
    fprintf(fid, '	  #endif\n');
    fprintf(fid, '	  */\n');
end
fprintf(fid, '    }\n\n');
	
fprintf(fid, '	// Outputs\n');
fprintf(fid, '    for (ii=0; ii<n_out;++ii) {\n');
fprintf(fid, '		// Simulink block output\n');
fprintf(fid, '		y[ii] = ssGetOutputPortRealSignal(S,ii);\n');
if Options.output_is_sparse
    fprintf(fid, '		// Output buffer for CasADi\n');
    fprintf(fid, '		res[ii] = (real_T*) ssGetDWork(S,ii); // So res[ii] wont point to nowhere\n');
end
fprintf(fid, '    }\n');
if Options.enable_tictoc_dspace
    fprintf(fid, '    #ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '        real_T* clk_tiTurnAround = ssGetOutputPortRealSignal(S,n_out);\n');
    fprintf(fid, '    #endif\n');
end	
fprintf(fid, '	//--------------------------//\n\n\n');
	
	
	
fprintf(fid, '      // //\n');
fprintf(fid, '	 //  //\n');
fprintf(fid, '	//////\n');
fprintf(fid, '		//\n');
fprintf(fid, '		//\n');
fprintf(fid, ' 	//////////////////////////////	\n');
fprintf(fid, '	/* Run the CasADi function  */\n');
fprintf(fid, '	//////////////////////////////\n');
if Options.output_is_sparse
    fprintf(fid, ['    ' functionName '( arg, res, iw, w, 0 );\n']);
else
    fprintf(fid, ['    ' functionName '( arg, y, iw, w, 0 );\n']);
end
fprintf(fid, '	//--------------------------//\n\n\n');
	

if Options.output_is_sparse && Options.output_is_matrix
    fprintf(fid, '	//////\n');
    fprintf(fid, '	//\n');
    fprintf(fid, '	/////\n');
    fprintf(fid, '		//\n');
    fprintf(fid, '	/////\n');
    fprintf(fid, '	//////////////////////////////\n');
    fprintf(fid, '	/* Decompress CCS format    */\n');
    fprintf(fid, '	//////////////////////////////\n');
    fprintf(fid, '	int_T col, row, jj, n_row, n_col, nnz, valuesInCol, kk;\n');
    fprintf(fid, '	const int_T* sp;\n');
    fprintf(fid, '	for (jj=0; jj<n_out; ++jj){\n');
    fprintf(fid, '		// Get sparsity information of casadi function output\n');
    fprintf(fid, ['		sp = ' functionName '_sparsity_out(jj);\n']);
    fprintf(fid, '		n_row = sp[0]; 				// Number of rows of output\n');
    fprintf(fid, '		n_col = sp[1]; 				// Number of columns of output\n');
    fprintf(fid, '		nnz   = sp[2+sp[1]];		// Number of nonzero elements\n');
    fprintf(fid, '      ii = 0;\n');
	fprintf(fid, '      col = 0;\n');
    fprintf(fid, '      valuesInCol = sp[2+col+1]-sp[2+col];\n');
	fprintf(fid, '      while(ii<nnz){\n');
    fprintf(fid, '        while(valuesInCol == 0){\n');
    fprintf(fid, '            col++;\n');
    fprintf(fid, '            valuesInCol = sp[2+col+1]-sp[2+col];\n');
    fprintf(fid, '        }\n');
    fprintf(fid, '        for(kk=0; kk<valuesInCol; kk++){\n');
    fprintf(fid, '            row = sp[2+n_col+1+ii];\n');
    fprintf(fid, '            y[jj][n_row*col+row] = res[jj][ii]; \n');
    fprintf(fid, '            ii++;\n');
    fprintf(fid, '        }\n');
    fprintf(fid, '        valuesInCol = 0;\n');
	fprintf(fid, '      }\n');
    fprintf(fid, '	}\n');
    fprintf(fid, '	//--------------------------//\n\n\n');
    
elseif Options.output_is_sparse
    fprintf(fid, '	//////\n');
    fprintf(fid, '	//\n');
    fprintf(fid, '	/////\n');
    fprintf(fid, '		//\n');
    fprintf(fid, '	/////\n');
    fprintf(fid, '	//////////////////////////////\n');
    fprintf(fid, '	/* Decompress CCS format    */\n');
    fprintf(fid, '	//////////////////////////////\n');
    fprintf(fid, '	int_T row, jj, nnz;\n');
    fprintf(fid, '	const int_T* sp;\n');
    fprintf(fid, '	for (jj=0; jj<n_out; ++jj){\n');
    fprintf(fid, '		// Get sparsity information of casadi function output\n');
    fprintf(fid, ['		sp = ' functionName '_sparsity_out(jj);\n']);
    fprintf(fid, '		nnz = sp[3];		// Number of nonzero elements if vector\n');
    fprintf(fid, '		for(ii=0; ii<nnz; ii++){\n');
    fprintf(fid, '			row = sp[4+ii]; // First row information if vector \n');
    fprintf(fid, '			y[jj][row] = res[jj][ii];\n');
    fprintf(fid, '		}\n');
    fprintf(fid, '	}	\n');
    fprintf(fid, '	//--------------------------//\n\n\n');    
end
	
if OptionsStd.enable_debug_dspace
    fprintf(fid, '	/*\n');
    fprintf(fid, '		#ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '		for(ii=0; ii<28; ii++){\n');
    fprintf(fid, '			msg_info_printf(1, 1, "Content of result: %f", res[0][ii]);\n');
    fprintf(fid, '		}\n');
    fprintf(fid, '	#endif\n');
    fprintf(fid, '	*/\n');
end

if Options.enable_tictoc_dspace
    fprintf(fid, '\n');
    fprintf(fid, '	// Store turnaround time in variable\n');
    fprintf(fid, '	#ifndef  MATLAB_MEX_FILE \n');
    fprintf(fid, '	  clk_tiTurnAround[0] = RTLIB_TIC_READ();\n');
    fprintf(fid, '	#endif\n\n\n');
end

fprintf(fid, '}\n\n\n');

%-------------------------------------------------------------------------%





%-------------------------------------------------------------------------%
%   Model terminate function                                              %
%-------------------------------------------------------------------------%
fprintf(fid, '/* Function: mdlTerminate =================================================\n');
fprintf(fid, ' * Abstract:\n');
fprintf(fid, ' *   In this function, you should perform any actions that are necessary\n');
fprintf(fid, ' *   at the termination of a simulation.\n');
fprintf(fid, ' */\n');
fprintf(fid, 'static void mdlTerminate(SimStruct *S)\n');
fprintf(fid, '{\n');
fprintf(fid, '}\n\n\n');

%-------------------------------------------------------------------------%





fprintf(fid, '/* Required S-function trailer */\n');
fprintf(fid, '#ifdef    MATLAB_MEX_FILE\n');
fprintf(fid, '# include "simulink.c"\n');
fprintf(fid, '#else\n');
fprintf(fid, '# include "cg_sfun.h"\n');
fprintf(fid, '#endif\n');




fclose(fid);