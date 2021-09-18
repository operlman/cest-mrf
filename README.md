The purpose of this code is to demonstrate the CEST/MT-MRF parameter quantification pipeline.

The folder contains data ('acquired_data.mat'), obtained at 9.4T,
where 3 vials of 50 mM L-arginine at pH 4, 4.5, and 5 were scanned
using a CEST-MRF protocol of 30 iterations. (See paper X).

The code contains two main parts:

A. MATLAB part
=============
run the file: --> run_demo.m <--
The code will guide you through the different steps to:
(1) Generate a CEST-MRF dictionary.
    - External packages (pulseq, yamlmatlab) will be automatically installed.
    - Parallel computation is performed while using the open Pulseq standard.
    For more info, see:
    Herz, K, Mueller, S, Perlman, O, et al. Pulseq-CEST: Towards multi-site multi-vendor
    compatibility and reproducibility of CEST experiments using an open-source sequence standard. 
    Magn Reson Med.; https://doi.org/10.1002/mrm.28825
(2) Perform dot-product matching.
(3) Install the packages for deep reconstruction.

B. Python part: 
===============
The deep reconstruction is performed using Python code.
It requires having python installed with the following packages:
numpy, scipy, matplotlib, and torch.
Suggested installation routes:")
1) Use pip (https://pip.pypa.io/en/stable/)")
   ---  OR ----
2) Use Anaconda (https://www.anaconda.com/products/individual-d)
   * a YAML file that allows creating the relevant environment
   is available in this folder: 'conda_environment.yml'
   ---  OR ----
3) Docker (https://www.docker.com/).
   * A docker-image with the required packages can be obtained by:
   'docker pull operlman/pytroch_scipy_matplotlib_scikit-image'
Once the packages are installed, run ---> 'deep_reco.py' <----
The script will use the file dict.mat, generated in the previous steps, as
well as the file acquired_data.mat, available in this folder.

If you use this code in a scientific publication, please cite: 
X.
