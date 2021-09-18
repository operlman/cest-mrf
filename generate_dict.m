%% A batch file that generates a yaml and seq file for a phantom dictionary and runs the dictionary generation

% The code is divided in 3 different sections:
% 1. Write .yaml file with info about the probe and what parameter to vary.
% 2. Write a CW CEST-MRF with single-shot readout .seq file.
% 3. Run the dictionary generation.

% Get filepath
if strcmp(mfilename, 'LiveEditorEvaluationHelperESectionEval')
    script_fp = fileparts(matlab.desktop.editor.getActiveFilename);
else
    script_fp = fileparts(which(mfilename));
end
yaml_fn = fullfile(script_fp, 'scenario.yaml');
seq_fn  = fullfile(script_fp, 'acq_protocol.seq');

% dict fn is optional. If none is specified the file has just the same name
% as th yaml_fn, just with a mat ending
dict_fn = fullfile(script_fp, 'dict.mat');

%% Write the yaml-file
% We vary T1w, T2w, Ksw and M0s in one file

% Small description
yaml_struct.description = 'A parameter file for Iohexol phantom CEST-MRF experiment';

% Water_pool
yaml_struct.water_pool.t1 = (2500:50:3300) ./ 1000; %vary t1
yaml_struct.water_pool.t2 = (600:50:1200) ./ 1000; % vary t2
yaml_struct.water_pool.f = 1;

% Solute pool
yaml_struct.cest_pool.Amine.t1 = 2800 ./ 1000;  % fix solute t1
yaml_struct.cest_pool.Amine.t2 = 40 ./ 1000; % fix solute t2
yaml_struct.cest_pool.Amine.k  = 100:10:1400;  % vary solute excjhabge rate
yaml_struct.cest_pool.Amine.f  = (10:5:120) .* 3 ./ 110000; % solute concentration * protons / water concentration
yaml_struct.cest_pool.Amine.dw = 3.0;  % fixed solute exchange rate at 3 ppm

% Fill initial magnetization info
% this is important now for the mrf simulation! For the regular pulseq-cest
% simulation, we usually assume athat the magnetization reached a steady
% state after the readout, which means we can set the magnetization vector
% to a specific scale, e.g. 0.5. This is because we do not simulate the
% readout there. For mrf we include the readout in in the simulation, which
% means we need to carry the same magnetization vector through the entire
% sequence. To avoid that the magnetization vector gets set to the initial
% value after each readout, we need to set reset_init_mag to false
yaml_struct.scale          = 1;
yaml_struct.reset_init_mag = 0;

% Fill scanner info
yaml_struct.b0       = 9.4;         % [T]
yaml_struct.gamma    = 267.5153;  % [rad / uT]
yaml_struct.b0_inhom = 0;
yaml_struct.rel_b1   = 1;

% Fill additional info
yaml_struct.verbose = 0;             % no unneccessary console log wished
yaml_struct.max_pulse_samples = 100; % block pulses are only 1 sample anyways

% Write the file
yaml.WriteYaml(yaml_fn, yaml_struct);

%% Write the seq file for a 2d experiment
% for more info about the seq file, check out the pulseq-cest repository
seq_defs.n_pulses      = 1       ; % number of pulses
seq_defs.tp            = 3   ; % pulse duration [s]
seq_defs.td            = 0   ; % interpulse delay [s]
seq_defs.Trec          = 1      ; % delay before readout [s]
seq_defs.Trec_M0       =  nan      ; % delay before m0 readout [s]
seq_defs.M0_offset     = nan   ; % dummy m0 offset [ppm]
seq_defs.DCsat         = seq_defs.tp/(seq_defs.tp+seq_defs.td); % duty cycle
seq_defs.offsets_ppm   = ones(30,1) .* 3.0;   % offset vector [ppm]
seq_defs.num_meas      = numel(seq_defs.offsets_ppm)   ;        % number of repetition
seq_defs.Tsat          = seq_defs.n_pulses.*(seq_defs.tp+ seq_defs.td)-seq_defs.td;
seq_defs.B0            = 9.4        ; % B0 [T]
[~, seqid] = fileparts(seq_fn);
seq_defs.seq_id_string = seqid    ; % unique seq id
% we vary B1 for the dictiobnary generation
seq_defs.B1pa          = [5, 5, 3, 3.75, 2.5, 1.75, 5.5, 6, 3.75, ...
    5.75, 0.25, 3, 6, 4.5, 3.75, 3.5, 3.5, 0, 3.75, 6, 3.75, 4.75, 4.5, ...
    4.25, 3.25, 5.25, 5.25, 0.25, 4.5, 5.25];

% >>> Gradients and scanner limits - see pulseq doc for more info
% Mostly relevant for clinical scanners
% lims =  mr.opts('MaxGrad',30,'GradUnit','mT/m',...
%     'MaxSlew',100,'SlewUnit','T/m/s', ...
%     'rfRingdownTime', 50e-6, 'rfDeadTime', 200e-6, 'rfRasterTime',1e-6);
% <<<

% gamma
gyroRatio_hz  = 42.5764;            % for H [Hz/uT]
gyroRatio_rad = gyroRatio_hz*2*pi;  % [rad/uT]


% This is the info for the 2d readout sequence. As gradients etc ar
% simulated as delay, we can just add a delay afetr the imaging pulse for
% simulation which has the same duration as the actual sequence
% the flip angle of the readout sequence:
imagingPulse = mr.makeBlockPulse(deg2rad(60), 'Duration', 2.1e-3);
% =='system', lims== should be added for clinical scanners

% the duration of the readout sequence:
te = 20e-3;
imagingDelay = mr.makeDelay(te);           

% init sequence
seq = SequenceSBB();
% seq = SequenceSBB(lims); for clinical scanners

% M0 pulse if required
% seq.addBlock(mr.makeDelay(seq_defs.Trec_M0)); add a delay for m0 
% seq.addBlock(imagingPulse);% add imaging block
% seq.addPseudoADCBlock();    % additional feature of the SequenceSBB class
% seq.addBlock(imagingDelay);

% Loop b1s
idx = 1;
for B1 = seq_defs.B1pa
    
    if idx > 1 % add relaxtion block after first measurement
        seq.addBlock(mr.makeDelay(seq_defs.Trec - te));% net recovery time
    end
    
    % saturation pulse
    currentOffsetPPM = seq_defs.offsets_ppm(idx);
    currentOffsetHz  = seq_defs.offsets_ppm(idx)*seq_defs.B0*gyroRatio_hz;
    fa_sat     = B1*gyroRatio_rad*seq_defs.tp; % flip angle of sat pulse
    
        % add pulses
        for np = 1:seq_defs.n_pulses
            
            % If B1 is 0 simulate delay instead of a saturation pulse
            if B1 ==0
                seq.addBlock(mr.makeDelay(seq_defs.tp));% net recovery time
            else
                
                satPulse   = mr.makeBlockPulse(fa_sat, 'Duration', seq_defs.tp,'freqOffset', currentOffsetHz);
                % =='system', lims== should be added for clinical scanners


                seq.addBlock(satPulse);
            end
            
            % delay between pulses
            if np < seq_defs.n_pulses 
                seq.addBlock(mr.makeDelay(seq_defs.td)); 
            end
        end
        
    % Spoiling
    % additional feature of the SequenceSBB class
    % seq.addSpoilerGradients(lims);
    
    % Imaging pulse
    seq.addBlock(imagingPulse);
    seq.addBlock(imagingDelay);
    seq.addPseudoADCBlock(); % additional feature of the SequenceSBB class

    idx = idx+1;
end

def_fields = fieldnames(seq_defs);
for n_id = 1:numel(def_fields)
    seq.setDefinition(def_fields{n_id}, seq_defs.(def_fields{n_id}));
end
seq.write(seq_fn);

%% Finally, lets run the dictionary generation

num_workers = 8; % Threads 
disp(['Assuming ',num2str(num_workers), ' threads are available'])
dict = generate_mrf_cest_dictionary(seq_fn, yaml_fn, dict_fn, num_workers);
