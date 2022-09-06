data
- dat.stim = stimulus information for each trial
- dat.key_down = key press
- dat.key_time = time point that key press was sampled from


plotBinoRegion.m
- analyze the data (from data folder) and make plots (prop fading and div-conv difference)


glm_preprocess.m
- take data and save as a table (.csv) for glm


glm_binShapevsAR.m
- binarize proportion fading into fade or no fade, then run logistic regression using MATLAB's fitglme function 
- also runs follow up analysis to explore interaction effect
