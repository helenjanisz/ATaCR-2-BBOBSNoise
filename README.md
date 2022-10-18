# ATaCR-2-BBOBSNoise
A set of codes for organizing and visualizing Broadband Ocean Bottom Seismometer (BBOBS) spectra output from running the ATaCR Package. Will reproduce figures and workflow from Janiszewski et al. (in review) https://eartharxiv.org/repository/view/3319/. 

For any use of these datasets or codes, please cite ____ 

Brief details for each of these codes is given below. Contact hajanisz@hawaii.edu for questions. Please note, these codes are designed to read in matfiles that are output from the ATaCR package. They also require a metadata table (defaults are either in excel or matlab format) to keep station information organized. The versions of these tables compatible with these scripts are given here for example formatting; however, the user is responsible for making tables like these of their own.

The scripts should also be able to be adapted to read in the text file versions of the data that we have made available to go along with our paper, but will require some editing. Text Files can be found here: Janiszewski, Helen et al. (2022), BBOBS_Noise_Properties_Review, Dryad, Dataset, https://doi.org/10.25349/D90042

This also contains scripts to format output ATaCR data to be compatible as input data for the spectral angle analysis available here: https://github.com/brennanbrunsvik/Ocean-bottom-seismometer-noise-clustering

OBS_TableParams.m - script that is called by multiple other scripts. This script is used for consistent formatting, keeping experiment names organized, specifying input and output directories, and more. Users would need to either edit or make additions to the experiment information and station color and symbol parameters for new analyses.

mk_spect_txt.m, mk_spec_corr_txt.m, mk_cpa_txt.m - codes that read in matfiles and export text files for spectra, coherence, phase, and admittance found at Janiszewski, Helen et al. (2022), BBOBS_Noise_Properties_Review, Dryad, Dataset, https://doi.org/10.25349/D90042

OBS_Working4Paper.xlsx/mat - excel and matlab versions of metadata table. Contains same information as TableS2, but some header information is different for script compatibility.

smoothSpectrum_octave.m - used under MIT License from here: https://github.com/IoSR-Surrey/MatlabToolbox

Make_Dat_SepctraNoCorr.m, Make_Dat_SepctraTiltCorr.m, Make_Dat_SepctraCorr.m - scripts for taking output data from ATaCR and making it compatible as input data for spectra angle analysis (https://github.com/brennanbrunsvik/Ocean-bottom-seismometer-noise-clustering). 

