function storeUnconstrainedSolution(ocp)
   
   % Save current path
     currPath = pwd;
     
   % Build temporary folder
     templateFolder = what('05_example');
     filePath = fullfile(templateFolder.path,'MATLAB','temp');
     if exist(filePath,'dir')
         rmpath(filePath);
         rmdir(filePath,'s');
     end
     mkdir(filePath);
     addpath(filePath);
     
   % Save unconstrained ocp
     fileName = 'ocp_unconstrained';
     cd(filePath);
     ocp.store_iterate(fileName);
     
   % Return to current path
     cd(currPath);

   end