%% binarized data, compare between shape and AR
clear all; close all;
% load in data
T = readtable('cont_resp_trial_data.csv');

% fields are: Subj    StimType    ConvDiv    Mfov    MonoReg    BinoReg    NumFade    TimePt    PropFade

% Calculate prop fade with full precision
T.PropFade = T.NumFade./T.TimePt;

% binarize PropFade
T.PropFadeBin = T.PropFade;
T.PropFadeBin(T.PropFade < 0.5) = 0;
T.PropFadeBin(T.PropFade >= 0.5) = 1;

% change oval and rect into shape
T.StimType(strcmp('oval',T.StimType))={'shape'};
T.StimType(strcmp('rect',T.StimType))={'shape'};
T.TotalFOV = T.BinoReg + 2*T.MonoReg;
keyboard

% run glme
glme_full = fitglme(T,'PropFadeBin ~ 1 + TotalFOV + StimType + ConvDiv + TotalFOV:StimType + TotalFOV:ConvDiv + ConvDiv:StimType + (1|Subj)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

% run glme
glme_full = fitglme(T,'PropFadeBin ~ 1 + BinoReg + StimType + ConvDiv + BinoReg:StimType + BinoReg:ConvDiv + ConvDiv:StimType + (1|Subj)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

%% follow up analysis

% subset data into just simple and just AR to examine interactions with StimType
T_shape = T(strcmp(T.StimType,'shape'),:);
T_AR = T(strcmp(T.StimType,'AR'),:);  
% follow up logistic to examine conv/div difference for simple and AR
glme_shape_convdiv = fitglme(T_shape,'PropFadeBin ~ 1 + ConvDiv + (1|Subj)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);  
glme_AR_convdiv = fitglme(T_AR,'PropFadeBin ~ 1 + ConvDiv + (1|Subj)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);  
% follow up logistic to examine binoreg difference for simple and AR
glme_shape_binoreg = fitglme(T_shape,'PropFadeBin ~ 1 + BinoReg + (1|Subj)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);  
glme_AR_binoreg = fitglme(T_AR,'PropFadeBin ~ 1 + BinoReg + (1|Subj)', ...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
