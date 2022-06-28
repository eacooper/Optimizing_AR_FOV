data file
- dat = data
- dat.stim = stimulus info for each trial
- dat.order = whether convergent was presented first or second
- dat.resp = subject response (they selected the first or second stim to be better)


plotBinoRegion.m
- analyze the data (from data folder) and make plots
- store the proportion data as result.mat

glm_ShapevsAR.m
- run logistic regression using MATLAB's fitglme function 
- also runs follow up analysis to explore interaction effects

do_chisq.m
- load result.mat and run chi square test
