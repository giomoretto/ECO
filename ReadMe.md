%% Info
Open source repository for ECO

- ECO stands for economic combustion optimization
- ECO solves for the optimal injector inputs of a direct-injection 
    compression-ignition engine
- The continuous-time OCP is formulated with the software packages acados 
    and casadi and discretized by an explicit fourth-order Runge-kutta 
    integration scheme and multiple shooting. The resulting NLP is solved 
    by sequential quadratic programming (SQP) using HPIPM with 
    a condensed horizon of 5. Furthermore, the Hessian is regularized by a 
    Levenberg-Marquardt Regularization.


%% Additional Software Packages
Please install additional software in the folder '04_external'
- casadi v.3.4.5 (https://web.casadi.org/)
- acados v0.1.9 or later (https://github.com/acados/acados)


%% Exemplary scripts
In order to run the exemplary scripts, please do the following

1) add all folders and subfolders to the working path 
    command window: run path_set.m 
2) create s-functions necessary to run ECO in simulink
    command window: run main_InjOpt.m
3) run exemplary file for use of ECO in matlab
    command window: run mainMATLAB.m
4) run exemplary file for use of ECO in Simulink
    command window: run mainSimulink.m
    

The exemplary file for matlab highlights how ECO can be used to obtain 
optimal injector inputs under varying constraints which are set by the 
user. Make sure to initialize the solver with initial values which are 
close to the solution previously obtained. Otherwise, robust convergence to
a feasible solution is not guaranteed.

The exemplary file for simulink highlights how the solutions of ECO can be
used in a real-time controller as feedforward control inputs. The bounds, 
which are not met, can be adapted by a simple PI-controller. The simulink 
model presented can be used for compilation on embedded systems


%% Cite the code
[![DOI](https://zenodo.org/badge/639391121.svg)](https://zenodo.org/badge/latestdoi/639391121)
