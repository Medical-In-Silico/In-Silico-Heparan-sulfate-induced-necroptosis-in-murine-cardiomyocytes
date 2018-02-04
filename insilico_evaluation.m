function [t, y, broken] = insilico_evaluation(odes, hs_fun, y0, x, required_timesteps)
%insilico_evaluation utilizes the MATLAB ode-solver for evaluation of the current in silico model.
%% args
% required_timesteps: time steps that will be evaluated

y_start = y0(1); % initial HS value, required for hs_fun

% ODE system setup to provide the insilico model
odes_handler = @(t,y) odes(t,y,x, @(t) hs_fun(t,y_start));

% ode solver options
opts = odeset('NonNegative', [1:length(y0)]', 'RelTol', 1e-6, 'AbsTol', 1e-6);

% call a suitable ODE solver
[t,y] = ode23(odes_handler, required_timesteps, y0, opts);

if length(t) ~= length(required_timesteps) % numerical simulation failed
    broken = true;
    y = 100*ones(length(required_timesteps),length(y0));
    t = required_timesteps;
else
    broken = false;
end

end