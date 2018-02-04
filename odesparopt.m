function odesparopt(errorfcn, xinit, modelfun, hsfun, bestx_file, cmaes_var, populationsize, randomSeed, MaxFunEvals)
%% starts parameter optimization for a given system of odes
% errorfcn - error function / target function that should be minimized
% xinit - initial parameter guess
% modelfun - system of ODEs for model simulation
% hsfun - additional function required for the model simulation
% bestx_file - result log file
% cmaes_var - initial variance of parameter mutations by CMA-ES
% populationsize - number of populations searching for a minimized solution
  
%%% CMAES
% modelfun = @necroptosis_odesystem
% hsfun = @hs_dotfun
%---------------------------------------------------------

%% Covariance Matrix Adapation - Evolutionary Strategy: instable. stops without reason
opts = {};
opts.LBounds = 0; opts.UBounds = 100; % bounded search space
opts.Seed = randomSeed; % set a random number generator seed
opts.PopSize = populationsize; % population size, hint: start with small populations (~5), increase if results not satisfactory
opts.EvalParallel = 'yes';
opts.MaxFunEvals = MaxFunEvals;
opts.EvalInitialX = 'no';
opts.Restarts = 3;  % each restart doubles the population size 
opts.LogFilenamePrefix = bestx_file;

% run CMAES
xmin_cmaes = cmaes(errorfcn, xinit, cmaes_var, opts, modelfun, hsfun);

% print optimized parameters
disp(xmin_cmaes')
end