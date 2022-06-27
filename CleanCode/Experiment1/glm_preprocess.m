%process data for logistic regression
clear all;
close all;
addpath('./data');

files = dir(fullfile('./data/*.mat'));
files = {files.name};
%stim info:[stimtype configuration monocular_field monocular_region scene icon];
central_measure = 0;%what to plot: mean = 0, median = 1; 
result = [];
stimtype = [1 2 3]; %oval, rect, AR
configtype = [1 2]; %1 = convergent, 2 = divergent;
mono_field = [40 30];
mono_zone = [9 4.5];
load(files{1});
stimOrder = dat.stim;
alldata = {};
%sort all the subj's data
for s = 1:size(files,2)
    load(files{s});
    subjdata = [dat.stim dat.key_down];
    %can use ismember to sort this each trial is unique (no exact repeat)
    [~,ind] = ismember(stimOrder,subjdata(:,1:6),'rows');
    reorder_subjdata = subjdata(ind,:);
    alldata = [alldata reorder_subjdata(:,7:end)];
end
result = [];
result_count = [];
nokey_result = [];

%calculate the proportion
tres = 150;
stim_data = [];
for n = 1:size(files,2)
    raw_resp = alldata{n};
    subjdata = [];

    for trial = 1:size(raw_resp,1)
        
        see_lune = find(raw_resp(trial,:) ==1);
        
        no_lune = find(raw_resp(trial,:) ==0);
        
        others = sum(raw_resp(trial,:) ==3) + sum(raw_resp(trial,:)==2) + sum(double(isnan(raw_resp(trial,:)))); %number of times no key response
        
        nume = size(see_lune,2);
        denom = (size(see_lune,2)+size(no_lune,2));

        prop = nume / denom;
        
        
        if others<tres %only add the trials that are valid, include binocular region size
            stim_data = [stim_data ;n,stimOrder(trial,1:4),stimOrder(trial,3)-stimOrder(trial,4),nume,denom,prop];
        end
    end
    
end

stim_data = sortrows(stim_data,[1 2 3 4]);

%% save as csv for glme, make each trial into individual rows
expandslist = [];
stimtype = {'oval','rect','AR'};
configtype = {'conv','div'};
subjid = {'S1','S2','S3','S4','S5','S6',...
    'S7','S8','S9','S10','S11','S12',...
    'S13','S14','S15','S16','S17','S18',...
    'S19','S20'};

subjlist = subjid(stim_data(:,1))';
stimlist = stimtype(stim_data(:,2))';
configlist = configtype(stim_data(:,3))';
mflist = stim_data(:,4);
mrlist = stim_data(:,5);
binolist = stim_data(:,6);
numpoint = stim_data(:,8);
numFade = stim_data(:,7);
proplist = stim_data(:,7)./stim_data(:,8);


T = table(subjlist,stimlist,configlist,mflist,mrlist,binolist,numFade,numpoint,proplist);
T.Properties.VariableNames = {'Subj','StimType','ConvDiv','Mfov','MonoReg','BinoReg','NumFade','TimePt','PropFade'};
writetable(T,'./cont_resp_trial_data.csv','Delimiter',',');

