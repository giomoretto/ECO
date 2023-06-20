   function loadUnconstrainedSolution(ocp)
   
   % Save current path
     currPath = pwd;
     
   % Load unconcstrained ocp
     templateFolder = what('05_example');
     filePath = fullfile(templateFolder.path,'MATLAB','temp');
     fileName = 'ocp_unconstrained';
     cd(filePath);
     ocp.load_iterate(fileName);
     
   % Return to current path
     cd(currPath);
     
   end