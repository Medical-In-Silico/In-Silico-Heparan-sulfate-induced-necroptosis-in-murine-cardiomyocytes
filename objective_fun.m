function [ sse_error ] = objective_fun(x, odes, hs_fun)
back2base_weight = 0; % if greater than 0, this is attracting the baseline after a long period of time
eval_timesteps = [0,4,8,16,24,48];

% place / component ...
%  1     2      3        4         5      6      7      8         9        10     ll     l2      l3
%'Heparan sulfate','pERK 1/2','Cytochrome C','cleaved Caspase 3','cleaved PARP','ERK 1/2','Caspase 3','PARP','Apoptosis','pRIP3','TNF alpha','RIP3','Necroptosis' % places

% training selection
% Use clPARP for training (1), for validation (0)
% l,2,3-row: HS 5, 10, 20; 1,2,3,4-col: 4,8,16,24h
select_training_samples_clPARP = [[1 1 0 0];[0 1 1 0];[0 0 1 1]];

%% initialize SSE array
sse_add = 100*ones(1,13);

% negative penalty, terminate directly as ODE solvers may be confused from
% imaginary parts in ode solutions (which might occur when negative
% parameters are applied)
negative_params = x(x<0);
negative_params_norm1 = norm(negative_params, 1);

if negative_params_norm1 > 0
    sse_error = 10^16*(1+negative_params_norm1);
    return
end

% penalty for to simulated values that exceeds a specified manifold of its intial value within the simulation. 
% '0' means no penalty.
upper_limit_penalty = 1;
upperBoundMultiplier = [...
    10 0 0 0 0 0 0 0 0 0 0 0 0;  ...HS
    0 4 0 0 0 0 0 0 0 0 0 0 0;   ...pERK
    0 0 4 0 0 0 0 0 0 0 0 0 0;   ...CytoC
    0 0 0 4 0 0 0 0 0 0 0 0 0;   ...cleavedCasp3
    0 0 0 0 4 0 0 0 0 0 0 0 0;   ...cleavedPARP
    0 0 0 0 0 4 0 0 0 0 0 0 0;   ...ERK
    0 0 0 0 0 0 4 0 0 0 0 0 0;   ...Casp3
    0 0 0 0 0 0 0 4 0 0 0 0 0;   ...PARP
    0 0 0 0 0 0 0 0 0 0 0 0 0;    ...Apoptosis
    0 0 0 0 0 0 0 0 0 4 0 0 0;   ...pRIP3
    0 0 0 0 0 0 0 0 0 0 4 0 0;   ...TNFa
    0 0 0 0 0 0 0 0 0 0 0 4 0;   ...RIP3
    0 0 0 0 0 0 0 0 0 0 0 0 0;    ...Necroptosis
];

lower_limit_penalty = 3;
lowerBound = [0  ...HS
    0.85   ...pERK
    0.85 ...Cytochrome C
    0.85...cleaved Caspase 3
    0.85...cleaved PARP
    0.85...ERK
    0.85...Caspase 3
    0.85...PARP
    0.85 ...Apoptosis
    0.85...pRIP3
    0.85...TNFa
    0.85...RIP3
    0.85 ]; %Necroptosis


% data - this not the original data.
% required to form the objective
clPARP_HS_5 = [1 1 1 1];
clPARP_HS_10 = [1 1 1 1];
clPARP_HS_20 = [1 1 1 1;

%% HS=0, h=16
y0 = [0 1 1 1 1 1 1 1 1 1 1 1 1]; % steady without HS induction

y_target = [0 1 1 1 1 1 1 1 1 1 1 1 1];
[t, y, broken] = insilico_evaluation(odes,hs_fun,y0, x, eval_timesteps);
if broken
    sse_error = nan;
    return
end
sse_add(1) = sum((y(3,[2:8,10:12])-y_target(1,[2:8,10:12])).^2);

if back2base_weight ~= 0
    y_target = [0 1 1 1 1 1 1 1 1 1 1 1 1];
    sse_add(2) = back2base_weight * sum((y(6,[2:5 10:11])-y_target([2:5 10:11])).^2);
else
    sse_add(2) = 0;
end        

if upper_limit_penalty ~= 0
    exceed_residual = (y - repmat((upperBoundMultiplier*y0')',size(y,1),1));
    exceed_residual(exceed_residual < 0) = 0;
    sse_add(3) = upper_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
else
    sse_add(3) = 0;
end

if lower_limit_penalty ~= 0
    exceed_residual = (y - repmat(lowerBound,size(y,1),1));
    exceed_residual(exceed_residual > 0) = 0;
    sse_add(3) = sse_add(3) + lower_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
end

%% HS=5
y0 = [5 1 1 1 1 1 1 1 1 1 1 1 1];
[t, y, broken] = insilico_evaluation(odes,hs_fun,y0, x, eval_timesteps);
if broken
    sse_error = nan;
    return
end

% for t->infinity, all back on baseline
if back2base_weight ~= 0
    y_target = [0 1 1 1 1 1 1 1 1 1 1 1 1];
    sse_add(4) = back2base_weight * sum((y(6,[2:5 10:11])-y_target([2:5 10:11])).^2);
else
    sse_add(4) = 0;
end        

if upper_limit_penalty ~= 0
    exceed_residual = (y - repmat((upperBoundMultiplier*y0')',size(y,1),1));
    exceed_residual(exceed_residual < 0) = 0;
    sse_add(5) = upper_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
else
    sse_add(5) = 0;
end 

if lower_limit_penalty ~= 0
    exceed_residual = (y - repmat(lowerBound,size(y,1),1));
    exceed_residual(exceed_residual > 0) = 0;
    sse_add(5) = sse_add(5) + lower_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
end

sse_add(6) = sum((y((ones(1,6) == [0 select_training_samples_clPARP(1,:) 0]), 5) - clPARP_HS_5(ones(1,4) == select_training_samples_clPARP(1,:))').^2);

%% HS=10
y0 = [10 1 1 1 1 1 1 1 1 1 1 1 1];
[t, y, broken] = insilico_evaluation(odes,hs_fun,y0, x, eval_timesteps);
if broken
    sse_error = nan;
    return
end
%% HS=10, h=16
y_target = [1 1 1 1 clPARP_HS_10(3) 1 1 1 1 1 1 1 1];
sse_add(7) = sum((y(4,2:4)-y_target(2:4)).^2) + sum((y(4,10:11)-y_target(10:11)).^2);

% for t->infinity, all back on baseline
if back2base_weight ~= 0
    y_target = [0 1 1 1 1 1 1 1 1 1 1 1 1];
    sse_add(8) = back2base_weight * sum((y(6,[2:5 10:11])-y_target([2:5 10:11])).^2);
else
    sse_add(8) = 0;
end        

if upper_limit_penalty ~= 0
    exceed_residual = (y - repmat((upperBoundMultiplier*y0')',size(y,1),1));
    exceed_residual(exceed_residual < 0) = 0;
    sse_add(9) = upper_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
else
    sse_add(9) = 0;
end 
if lower_limit_penalty ~= 0
    exceed_residual = (y - repmat(lowerBound,size(y,1),1));
    exceed_residual(exceed_residual > 0) = 0;
    sse_add(9) = sse_add(9) + lower_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
end 

sse_add(10) = sum((y(ones(1,6) == [0 select_training_samples_clPARP(2,:) 0], 5) - clPARP_HS_10(ones(1,4) == select_training_samples_clPARP(2,:))').^2);

%% HS=20
y0 = [20 1 1 1 1 1 1 1 1 1 1 1 1];
[t, y, broken] = insilico_evaluation(odes,hs_fun,y0, x, eval_timesteps);

if broken
    sse_error = nan;
    return
end

% h = 36, t->infinity, all back on baseline
if back2base_weight ~= 0
    y_target = [0 1 1 1 1 1 1 1 1 1 1 1 1];
    sse_add(11) = back2base_weight * sum((y(6,[2:5 10:11])-y_target([2:5 10:11])).^2);
else
    sse_add(11) = 0;
end

if upper_limit_penalty ~= 0
    exceed_residual = (y - repmat((upperBoundMultiplier*y0')',size(y,1),1));
    exceed_residual(exceed_residual < 0) = 0;
    sse_add(12) = upper_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
else
    sse_add(12) = 0;
end

if lower_limit_penalty ~= 0
    exceed_residual = (y - repmat(lowerBound,size(y,1),1));
    exceed_residual(exceed_residual > 0) = 0;
    sse_add(12) = sse_add(12) + upper_limit_penalty * sum(exceed_residual([2:8 10,11,12]).^2);
end

sse_add(13) = sum((y(ones(1,6) == [0 select_training_samples_clPARP(3,:) 0], 5) - clPARP_HS_20(ones(1,4) == select_training_samples_clPARP(3,:))').^2);

sse_error = sum(sse_add); % sum up all square errors
end