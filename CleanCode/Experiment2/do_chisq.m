
clear all;
%load data
load('result.mat');

%binarize participant as prefer convergent (>0.5) or divergent(<=0.5)
binarydata = propData; 
binarydata(:,3:22) = propData(:,3:22)>0.5;

stimtype = {'simple',[],'AR'};
statresult = propData(:,1:2); %stimulus info for stat

%out of 20 people, how many prefer conv.
for stims = 1:size(binarydata,1)

    conv = sum(binarydata(stims,3:end));
    
    obsCounts = [conv 20-conv]; %observed count
    expCounts = [10 10]; %expect by chance, half of participant would prefer each configuration
    [h,p,st] = chi2gof([0 1],'NBins',2,...
                        'Frequency',obsCounts, ...
                        'Expected',expCounts);
    stmtxt = [stimtype{binarydata(stims,1)} ',binoregion:' num2str(binarydata(stims,2))];
    %disp(stmtxt)
    statresult(stims,3:6) = [conv 20-conv st.chi2stat p];
end

disp('stimType  binoSize  Nconv  Ndiv  chisq  p')
disp(statresult)
