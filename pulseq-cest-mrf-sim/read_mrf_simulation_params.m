%% read yaml parameter file
%
% Input:  yaml_fn: filename of the .yaml parameter file
%
% Output: PMEX:    input struct for mex function
%         dict:    dictionary variable, for now filled with parameters
%                  combinations
% Modified by N.Vladimirov 2023 to generate multipool case

function [PMEX, dict] = read_mrf_simulation_params(yaml_fn)

% init output
PMEX = [];
dict = [];

%% check for file
if ~exist(yaml_fn, 'file')
    error('yaml parameter file does not exist!')
end

%% read struct
params = yaml.ReadYaml(yaml_fn);

%% get water pool
if ~isfield(params, 'water_pool')
    error('Water pool must be defined in "water_pool"');
end
wp = params.water_pool;
if ~isfield(wp, 'f') || ~isfield(wp, 't1') || ~isfield(wp, 't2')
    error('"water_pool" must contain "f", "t1" and "t2"');
end

% fisrt fill possible dict variables
dict.variables.t1w = get_variable_param(wp.t1);
dict.variables.t2w = get_variable_param(wp.t2);

% add to pmex
PMEX.WaterPool.f  = wp.f;
PMEX.WaterPool.R1 = 1.0/dict.variables.t1w(1);
PMEX.WaterPool.R2 = 1.0/dict.variables.t2w(1);


%% CEST pools
num_pools = 0;
if isfield(params, 'cest_pool')
    cp = params.cest_pool;
    pool_names = fieldnames(cp);
    num_pools = numel(pool_names);
    PMEX.CESTPool(num_pools) = struct('R1',[],'R2',[],'f',[],'k',[],'dw',[]);

    % only 1 cest pool so far
    for p = 1:num_pools
        cpool = cp.(pool_names{p});
        if ~isfield(cpool, 'f') || ...
                ~isfield(cpool, 't1') || ~isfield(cpool, 't2') || ...
                ~isfield(cpool, 'k') || ~isfield(cpool, 'dw')
            error([pool_names{p} ' must contain "f/c", "t1" , "t2", "k" and "dw"']);
        end
        
        % Dynamic field names based on pool number
        t1_field_name = sprintf('t1s_%d', p);
        t2_field_name = sprintf('t2s_%d', p);
        f_field_name = sprintf('fs_%d', p);
        k_field_name = sprintf('ksw_%d', p);
        
        % fill the dictionary struct
        dict.variables.(t1_field_name) = get_variable_param(cpool.t1);
        dict.variables.(t2_field_name) = get_variable_param(cpool.t2);
        dict.variables.(f_field_name)  = get_variable_param(cpool.f);
        dict.variables.(k_field_name)  = get_variable_param(cpool.k);
        
        % fill the pmex struct
        PMEX.CESTPool(p).R1 = 1.0/dict.variables.(t1_field_name)(1);
        PMEX.CESTPool(p).R2 = 1.0/dict.variables.(t2_field_name)(1);
        PMEX.CESTPool(p).f  = dict.variables.(f_field_name)(1);
        PMEX.CESTPool(p).k  = dict.variables.(k_field_name)(1);
        PMEX.CESTPool(p).dw = cpool.dw;

    end
else
    warning('No CEST pools found in param files! specify with "cest_pool"');
end

%% MT pool
if isfield(params, 'mt_pool')
    mt = params.mt_pool;
    if ~isfield(mt, 'f')  || ...
            ~isfield(mt, 't1') || ~isfield(mt, 't2') || ...
            ~isfield(mt, 'k') || ~isfield(mt, 'dw') || ~isfield(mt, 'lineshape')
        error('"mt_pool" must contain "f", "t1", "t2", "k", "dw" and "lineshape"');
    end
    if ~strcmp(mt.lineshape, 'SuperLorentzian') && ~strcmp(mt.lineshape, 'Lorentzian') && ...
            ~strcmp(mt.lineshape, 'None')
        error([mt.lineshape ' is invalid. Please use "None", "Lorentzian" or "SuperLorentzian"']);
    end
    
    % fill the dictionary struct
    dict.variables.t1ss = get_variable_param(mt.t1);
    dict.variables.t2ss = get_variable_param(mt.t2);
    dict.variables.fss = get_variable_param(mt.f);
    dict.variables.kssw = get_variable_param(mt.k);
      
    % fill the pmex strut
    PMEX.MTPool.R1 = 1.0/dict.variables.t1ss(1);
    PMEX.MTPool.R2 = 1.0/dict.variables.t2ss(1);
    PMEX.MTPool.f  = dict.variables.fss(1);
    PMEX.MTPool.k  = dict.variables.kssw(1);
    PMEX.MTPool.dw = mt.dw;
    PMEX.MTPool.Lineshape = mt.lineshape;
    
else
    warning('No MT pool found in param files! specify with "mt_pool"');
end


%% Put together an initial Magnetization vector (fully relaxed)
% [MxA, MxB, MxD, MyA, MyB, MyD, MzA, MzB, MzD, MzC]
% -> A: Water Pool, B: 1st CEST Pool, D: 2nd CEST Pool, C: MT Pool
% Cest pools would continue in the same way with E, F, G ...
nTotalPools = num_pools+1; % cest + water
PMEX.M = zeros(nTotalPools*3,1);
PMEX.M(nTotalPools*2+1,1)= PMEX.WaterPool.f;
for ii = 2:nTotalPools
    PMEX.M(nTotalPools*2+ii,1)= PMEX.CESTPool(ii-1).f;
end
if isfield(PMEX, 'MTPool') && size(PMEX.M,1) == nTotalPools*3 % add MT pool
    PMEX.M = [PMEX.M; PMEX.MTPool.f];
end

%% scale init vector
if isfield(params, 'scale')
    PMEX.M = PMEX.M * params.scale;
end

%% scanner parameters
if ~isfield(params, 'b0') || ~isfield(params, 'gamma')
    error('Parameter file must contain "b0" and "gamma"');
end
PMEX.Scanner.B0    = params.b0;    % field strength [T]
PMEX.Scanner.Gamma = params.gamma; % gyromagnetic ratio [rad/uT]

if isfield(params, 'b0_inhom')
    PMEX.Scanner.B0Inhomogeneity = params.b0_inhom;
end

if isfield(params, 'rel_b1')
    PMEX.Scanner.relB1 = params.rel_b1;
end


%% more optinal paramters
if isfield(params, 'verbose')
    PMEX.Verbose = logical(params.verbose);
end
if isfield(params, 'reset_init_mag')
    PMEX.ResetInitMag = logical(params.reset_init_mag);
end
if isfield(params, 'max_pulse_samples')
    PMEX.MaxPulseSamples = params.max_pulse_samples;
end

% for some reason yamlmatlab can read arrays into cells...
    function p = get_variable_param(p)
        if (iscell(p))
            p = cell2mat(p);
        end
    end
end

