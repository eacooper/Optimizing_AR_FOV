This folder contains the data and plotting/analysis code for Experiment 2. The data are in the data subfolder, with one file per participant. The data file format is as follows:

- dat = data structure
- dat.stim = stimulus information for each trial
	- col 1 = stimulus type (1 = oval, 2 = rect, 3 = AR)
	- col 2 = configuration (1 = convergent, 2 = divergent)
	- col 3 = monocular field size (30 or 40 deg)
	- col 4 = monocular region size (4.5 or 69 deg)
	- col 5 = AR background scene number
	- col 6 = AR icon set number
- dat.key_down = key press
- dat.key_time = time point that key press was sampled from


plotBinoRegion.m
- stores trial by trial data as trial_trial_data_for_GLM.csv
- analyzes the data and make plots (Figure 5 prop fading and Figure 7 div-conv difference)


glm_binShapevsAR.m
- binarize proportion fading into fade or no fade, then run logistic regression using MATLAB's fitglme function 
- also runs follow up analysis to explore interaction effect
