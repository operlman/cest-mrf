% MRF pipeline demonstration
% Kai Herz and Or Perlman 2021 
% kai.herz@tuebingen.mpg.de; operlman@mgh.harvard.edu

%% Step 1 - Generate a dictionary

% Installing packages 
install_pulseq_cest_mrf 

% Generate dictionary
t1=clock;
generate_dict
t2=clock;
disp(['Dictionary generation took ',num2str(etime(t2,t1)/60), ' min'])

disp('Three files should now be available in the dictionary_generation folder:')
disp('1. Acquisition protocol (acq_protocol.seq)')
disp('2. Simulated biological scenario (scenario.yaml)')
disp('3. Dictionary (dict.mat)')

clearvars

%% Dot product matching

% Load example raw data (phantom scanned at 9.4T preclinical scanner)
load acquired_data.mat acquired_data

% Load the generated dictionary
load dict.mat dict

t1=clock;
quant_maps = dot_prod_matching(dict, acquired_data);
t2=clock;
disp(['Dot product matching took ',num2str(etime(t2,t1)), ' sec'])

% Plotting dot-product matched maps

% Setting graphic parameters
set(0, 'DefaultAxesLineWidth', 1.2, 'DefaultAxesFontSize', 12, ...
          'DefaultAxesFontWeight', 'bold', 'DefaultAxesFontname','Times New Roman',...
          'DefaultLineLineWidth', .2, 'DefaultLineMarkerSize', 8);
set(0,'defaultfigurecolor',[1 1 1])

figure
ax1 = subplot(121);
imagesc(quant_maps.fs.*110e3/3, [0, 120]);
colormap(ax1,parula)
colorbar
axis off
title('[L-arg] (mM)')

ax2 = subplot(122);
imagesc(quant_maps.ksw, [0 500]);
colormap(ax2, hot)
colorbar
axis off
title('k_{sw} (Hz)')

%% Deep reconstruction

disp("=================================")
disp("The deep reconstruction is performed using a Python code.")
disp("It requires having python installed with the following packages:")
disp("numpy, scipy ,matplotlib, and torch.")
disp("Suggested installation routes:")
disp("1) Use pip (https://pip.pypa.io/en/stable/)")
disp("   ---  OR ----    ")
disp("2) Use Anaconda (https://www.anaconda.com/products/individual-d)")
disp("   * a YAML file that allows creating the relevant environment")
disp("   is available in this folder: 'conda_environment.yml'")
disp("   ---  OR ----    ")
disp("3) Docker (https://www.docker.com/).")
disp("   * A docker-image with the required packages can be obtained by:")
disp("   'docker pull operlman/pytroch_scipy_matplotlib_scikit-image'")
disp("=================================")
disp("Once the packages are installed, run ---> 'deep_reco.py' <----")
disp("The script will use the file dict.mat, generated in the previous steps, as")
disp("well as the file acquired_data.mat, available in this folder.")

