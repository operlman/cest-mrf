% This function runs the dictionary generation for a specific .seq and
% .yaml pair and saves the dictionary in a mat file
%
% Input:  seq_fn:      filename of the .seq-file
%         param_fn:    filename of the .yaml parameter file
%         dict_fn:     (optional) filename of resulting dictionary .mat file 
%         num_workers: (optional) number of workers for parpool
%

function dict = generate_mrf_cest_dictionary(seq_fn, param_fn, dict_fn, num_workers)
% if no input, let the user use the ui
if nargin < 2
    [seq_fn, seq_fp] = uigetfile({'*.seq','All .seq Files'},'Choose .seq-file for simulation');
    seq_fn = fullfile(seq_fp, seq_fn);
    [param_fn, param_fp] = uigetfile({'*.yaml; *.yml','All .yaml Files'},'Choose .yaml-file for simulation');
    param_fn = fullfile(param_fp, param_fn);
end

% dict file same as yaml file if not specified
if nargin < 3
    [p_fp, p_fn] = fileparts(param_fn);
    dict_fn = fullfile(p_fp, [p_fn '.mat']);
end

%% check for files
if ~exist(seq_fn, 'file')
    error('.seq file does not exist!')
end

%% read .yaml file
[PMEX, dict] = read_mrf_simulation_params(param_fn);

%% all variable parameters are now in dict.variables
% for now, lets use the old naming convention, but we can make it flexible
% later
% calculate entire number of combinations
if ~isfield(dict, 'variables')
    error('No parameter variation in yaml file...');
end

[dict, numComb] = prepare_dictionary(dict);

%% check adc events in seq file for allocation
seq = mr.Sequence();
seq.read(seq_fn);
numADC = 0; % loop through sequence and check for ADC blocks -> number of images
for cBlock= 1:numel(seq.blockEvents)
    if ~isempty(seq.getBlock(cBlock).adc)
        numADC = numADC+1;
    end
end
disp(['Found ' num2str(numADC) ' ADC events in seq file.']);

%% check number of pools for M indexing
nTotalPools = 1;
if isfield(PMEX, 'CESTPool')
    nTotalPools = nTotalPools + numel(PMEX.CESTPool);
end

%% Set up parallel pool

if isempty(gcp('nocreate'))
    % no running parallel pool found, start one
    if nargin == 4
        pp = parcluster;
        pp.NumWorkers = num_workers;
    else
        pp= parpool; % standard profile
    end
else
    pp = gcp; % get current running pool
end

%% distribute combinations across workers
numWorkers = pp.NumWorkers;
combIds = 1:numComb;
workerIds = mat2cell(combIds(:), diff(fix(linspace(0, numComb, numWorkers+1))), 1);

%% start parallel dictionary generation

parfor w = 1:pp.NumWorkers
    PMEX_local = PMEX;
    dict_local = dict;
    pulseq_cest_mex('init', PMEX_local, seq_fn);
    waterSignalPar{w} = zeros(numADC,numel(workerIds{w}));
    idx = 1;
    for c = workerIds{w}(:)'
        PMEX_local.WaterPool.R1 = 1.0/dict_local.t1w(c);
        PMEX_local.WaterPool.R2 = 1.0/dict_local.t2w(c);
        if isfield(PMEX_local, 'CESTPool')
            PMEX_local.CESTPool.R1 = 1.0/dict_local.t1s(c);
            PMEX_local.CESTPool.R2 = 1.0/dict_local.t2s(c);
            PMEX_local.CESTPool.f = dict_local.fs(c);
            PMEX_local.CESTPool.k = dict_local.ksw(c);
        end
        if isfield(PMEX_local, 'MTPool')
            PMEX_local.MTPool.R1 = 1.0/dict_local.t1ss(c);
            PMEX_local.MTPool.R2 = 1.0/dict_local.t2ss(c);
            PMEX_local.MTPool.f = dict_local.fss(c);
            PMEX_local.MTPool.k = dict_local.kssw(c);
        end
        pulseq_cest_mex('update', PMEX_local);
        Mout = pulseq_cest_mex('run');
        waterSignalPar{w}(:,idx) = sqrt(Mout(1,:).^2+Mout(nTotalPools+1,:).^2);
        idx = idx+1;
    end
    pulseq_cest_mex('close');
end

%% Save signal in dict
dict.sig = cell2mat(waterSignalPar);

% Removing unnecessary fields
dict = rmfield(dict,'variables');

%% save dict
save(dict_fn, 'dict');

