%Use a ChiSq test to ask whether or not there's equal number of people preferring convergent vs. divergent
clear all;

%load the proportion convergent chosen data
load('results_for_ChiSq.mat');
%binarize participant to convergent preferer
binarydata = propData; %if prop>0.5, this person prefers convergent more (1), else divergent prefer (0)
binarydata(:,3:22) = propData(:,3:22)>0.5;

stimtype = {'simple',[],'AR'};
statresult = propData(:,1:2); %stimulus info for stat

for stims = 1:size(binarydata,1)

    conv = sum(binarydata(stims,3:end));
    
    obsCounts = [conv 20-conv]; %out of 20 people, count how many prefer conv.
    expCounts = [10 10];
    [h,p,st] = chi2gof([0 1],'NBins',2,...
                        'Frequency',obsCounts, ...
                        'Expected',expCounts);
    stmtxt = [stimtype{binarydata(stims,1)} ',binoregion:' num2str(binarydata(stims,2))];
    disp(stmtxt)
    statresult(stims,3:6) = [conv 20-conv st.chi2stat p];
    disp(statresult(stims,:))
end
