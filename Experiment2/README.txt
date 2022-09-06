This folder contains the data and plotting/analysis code for Experiment 2. The data are in the data subfolder, with one file per participant. The data file format is as follows:

- dat = data structure
- dat.stim = stimulus info for each trial
	- col 1 = stimulus type (1 = oval, 2 = rect, 3 = AR)
	- col 3 = monocular field size (30 or 40 deg)
	- col 4 = monocular region size (4.5 or 9 deg)
	- col 5 = AR scene number
	- col 6 = AR icon number 
- dat.order = the order in which convergent configuration was shown (whether convergent was presented first or second)
- dat.resp = subject response (participants selected the first or second stimulus that they saw was better)

plotBinoRegion
- analyzes the data (from data folder) and make plots (Figure 7A)

glm_ShapevsAR
- process data for logistic regression
- run logistic regression using MATLAB's fitglme function 
- also runs follow up analysis to explore interaction effects

doChiSq
- follow up chi square test to see if people's preference for convergent over divergent and vice versa was more than chance level
