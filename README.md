# CEST/MT MR-Fingerprinting


**The purpose of this code** is to demonstrate the CEST/MT-MRF parameter quantification pipeline.

The folder contains data ('acquired_data.mat'), obtained at 9.4T, where 3 vials of 50 mM L-arginine at pH 4, 4.5, and 5 were scanned using a CEST-MRF protocol of 30 iterations. 

# The code contains two main parts:

## A. MATLAB part
run the file:  **run_demo.m** 

The code will guide you through the different steps to:

(1) Generate a CEST-MRF dictionary.
 *   External packages (pulseq, yamlmatlab) will be automatically installed.
 *   Parallel computation is performed while using the open Pulseq standard: https://pulseq-cest.github.io/
 *   For more info, see: Herz, K, Mueller, S, Perlman, O, et al. Pulseq-CEST: Towards multi-site multi-vendor compatibility and reproducibility of CEST experiments using an open-source sequence standard. Magn Reson Med.; https://doi.org/10.1002/mrm.28825
 
(2) Perform dot-product matching.

(3) Install the packages for deep reconstruction.

## B. Python part: 
The deep reconstruction is performed using Python code.
It requires having python installed with the following packages:
numpy, scipy, matplotlib, and torch.
Suggested installation routes:

1) pip (https://pip.pypa.io/en/stable/)

2) Anaconda (https://www.anaconda.com/products/individual-d)
* a YAML file that allows creating the relevant environment is available in this folder: 'conda_environment.yml'

3) Docker (https://www.docker.com/).
* A docker-image with the required packages can be obtained by: 

'docker pull operlman/pytroch_scipy_matplotlib_scikit-image'
   
Once the packages are installed, run **deep_reco.py**

The script will use the file dict.mat, generated in the previous steps, as well as the file acquired_data.mat, available in this folder.

### This repository is associated with a review paper called: "Deep Learning Based MR Fingerprinting for Semisolid Magnetization Transfer and Chemical Exchange Saturation Transfer Quantification". Additional details will be provided soon. 

**In the meanwhile, if you use this code in a scientific publication please cite one or more of the following relevant papers:**
1) Perlman, O, Herz, K, Zaiss, M, Cohen, O, Rosen, MS, Farrar, CT. CEST MR-Fingerprinting: Practical considerations and insights for acquisition schedule design and improved reconstruction. Magn Reson Med. 2020; 83: 462– 478. https://doi.org/10.1002/mrm.27937 

2) Herz, K, Mueller, S, Perlman, O, et al. Pulseq-CEST: Towards multi-site multi-vendor compatibility and reproducibility of CEST experiments using an open-source sequence standard. Magn Reson Med. 2021; 86: 1845– 1858. https://doi.org/10.1002/mrm.28825 
 
3) Perlman, O., Ito, H., Herz, K., Shono, N., Nakashima, H., Zaiss, M., Chiocca, E.A., Cohen, O., Rosen, M.S. and Farrar, C.T., 2020. AI boosted molecular MRI for apoptosis detection in oncolytic virotherapy. bioRxiv.

4) Cohen, O, Huang, S, McMahon, MT, Rosen, MS, Farrar, CT. Rapid and quantitative chemical exchange saturation transfer (CEST) imaging with magnetic resonance fingerprinting (MRF). Magn Reson Med. 2018; 80: 2449– 2463. https://doi.org/10.1002/mrm.27221 
