%% do glme on oval vs rectangle stimulus
clear all;

%combine all participant data
addpath('./data');
N = 20;

%load the first participant data and use its trial order to sort all other
%participant data
load('POL_1.mat');
stimOrder = dat.stim;
alldata = [];
%load each participant's data
for s = 1:N
    load(['POL_',num2str(s),'.mat']);
    %dat.order = whether convergent is first or second
    choseConverge = dat.resp==dat.order; %if 1, chose convergent(dat.order = the order of the convergent stimulus presented 
    subjdata = [ones(size(dat.stim,1),1)*s,dat.stim choseConverge];
    
    alldata = [alldata; subjdata];
end

%for storing data
subjid = {};%col1
stimtype = {}; %col2
monofield = []; %col4
monozone = []; %col5
resp = [];%col8

%coding the stimulus type as simple vs AR
for n = 1:size(alldata,1)
    subjid{n,1} = ['S',num2str(alldata(n,1))];
    
    if alldata(n,2)==1
        ST = 'simple';
    elseif alldata(n,2)==2
        ST = 'simple';
    elseif alldata(n,3)==3
        ST = 'AR';
    end
    
    stimtype{n,1} = ST;
    
    monofield(n,1) = alldata(n,4); %monocular FOV
    monozone(n,1) = alldata(n,5);  %monocular region
    binozone(n,1) = alldata(n,4) - alldata(n,5); %binocular overlap region
    resp(n,1) = alldata(n,8);
end

stimtype = categorical(stimtype);
stimtype = reordercats(stimtype,{'simple','AR'});

T = table(subjid,stimtype,monofield,monozone,binozone,resp);

%main GLME
glme_full = fitglme(T,'resp ~ 1 + stimtype + binozone + stimtype*binozone + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

%% follow up interaction analysis
%separate by stimulus (simple vs AR)
T_simple = T(T.stimtype =='simple',:);
T_AR = T(T.stimtype =='AR',:);

glme_shape = fitglme(T_simple,'resp ~ 1 + binozone + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

glme_AR = fitglme(T_AR,'resp ~ 1 + binozone + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

% separate by 2 levels of bino region (SMALL VS. LARGE)
binozone2 = {};
binozone2(T.binozone<30,1)= {'small bino'};
binozone2(T.binozone>=30,1)= {'large bino'};
T.binozone2 = binozone2;

T_smallb = T(strcmp(T.binozone2,'small bino'),:);
T_largeb = T(strcmp(T.binozone2,'large bino'),:);

glme_sbino = fitglme(T_smallb,'resp ~ 1 + stimtype + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

glme_lbino = fitglme(T_largeb,'resp ~ 1 + stimtype + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

