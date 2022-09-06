These scripts calculate and plot the perceptual fading predictions described in the Discussion section, which combines the geometric analysis and the data from Exp 1

avgdata.mat
- average data from observers in Experiment 1, used to fit a piecewise function (obtained from plot_diff_xaxis function in Experiment 1 folder)
	- col 1 = stimulus type (1=simple, 3 = AR)
	- col 2 = stimulus configuration (1 = convergent, 2 = divergent)
	- col 3 = binocular overlap size
	- col 4 = total FOV
	- col 5 = percentage of monocular region/monocular field

extrapFit
- extrapolate and fit a picewise function to the average data from Experiment 1 and plot it (Figure 12A)

predict_fading
- combine the perceptual fitting and the geometric analysis to predict the amount of fading for the three different camera configurations across different distances (Figure 12B)