%% Figure 7A
% This script plots the simple vs AR stimuli, as in Figure 7A. We plot the average and 95% CI of the prop fading time
% as a function of binocular region size. Individual subject data are also
% plotted

clear all;
close all;

%% PROCESSING DATA
addpath('./data');
%load the first participant data and use its trial order to sort all other
%participant data, since each row is unique, can just use default sortrows
%to get the same ordering
load('POL_1.mat');
stimOrder = sortrows(dat.stim);
alldata = stimOrder;

%total number of participants
N = 20;

% for each participant
for s = 1:N
    
    %load each participant's data
    load(['POL_',num2str(s),'.mat']);
    % add response (whether the chosen stimuli was convergent or not) to the stimulus information
    subjdata = [dat.stim, dat.resp==dat.order]; 
    %sortrows to match P1's stim order
    reorder_subjdata = sortrows(subjdata); %rows sorted based on stim
    
    %sanity check for trial order
    if reorder_subjdata(:,1:6) ~= stimOrder
        disp('something is wrong')
    end
    
    alldata = [alldata reorder_subjdata(:,7)];

end

propData = []; %create matrix to store stimulus + all participant's proportion
r = 1;
for type = [1 3] %simple vs AR
    for mono_field = [30 40]
        for mono_zone = [4.5 9]

            if type ==1    %Oval and Rect in data (1 and 2) are recoded as 1 (simple shape)         
                ind = find((alldata(:,1)==1 | alldata(:,1)==2) & alldata(:,3)==mono_field & alldata(:,4)==mono_zone);
            else 
                ind = find(alldata(:,1)==type & alldata(:,3)==mono_field & alldata(:,4)==mono_zone);
            end
            subdata = alldata(ind,:);
            
            %calculate the proportion of times that convergent was chosen for each subj by dividing by the
            %number of repeats
            propConv = sum(subdata(:,7:end),1)./length(ind);
            bino_reg = mono_field - mono_zone; %calculate the binocular region size
            propData(r,:) = [type,bino_reg,propConv];
            
            r = r+1;
        end
    end
end

%save the proportion for chi sq test later
save('results_for_ChiSq.mat','propData');


%% MAKE THE PLOT
figure();
hold on;
markers = {'rs-','bd-'};
colors = [[228 30 38]*0.5;51 127 186];
markers2 = {'rs','bd',};
fillc = colors./255;
offset = [2 4 6 8]';
c = 1;

for st = [1 3] %simple vs AR
    subdata = propData(propData(:,1)==st,:);
    subdata = sortrows(subdata,2); %sort based on binocular region size
    %each column is a subj. each row is a trial
    %get the mean and standard deviation across subjects for each trial
    means = mean(subdata(:,3:end),2);
    stds = std(subdata(:,3:end),[],2); 
    z = 1.96;
    n = sum(~isnan(subdata(:,3:end)),2); %Ns should be 20ppl since there's no missing data (due to thresholding) from this expt.
    if n ~=20
        disp('something is wrong')
    end
    %x and y value for ploting the shaded standard error
    xvals = [subdata(:,2);flip(subdata(:,2))];
    yvals = [means+z*stds./sqrt(n); flip(means-z*stds./sqrt(n))];
    fill(xvals,yvals,fillc(c,:),'LineStyle','none','HandleVisibility','off','FaceAlpha',0.05);
    %plot the mean
    plot(subdata(:,2),means,markers{c},'Color',fillc(c,:),'MarkerFaceColor',fillc(c,:));
    % for each bino overlap, plot raw data
    plot(subdata(1,2)*ones(1,20)-0.5+1*rand(1,20),subdata(1,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    plot(subdata(2,2)*ones(1,20)-0.5+1*rand(1,20),subdata(2,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    plot(subdata(3,2)*ones(1,20)-0.5+1*rand(1,20),subdata(3,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    plot(subdata(4,2)*ones(1,20)-0.5+1*rand(1,20),subdata(4,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    
    c = c+1;
end
xlabel('Binocular region (deg)');
ylabel('Prop convergent chosen');
ylim([0 1.2]);
xlim([20 36.5]);
yticks([0 0.25 0.5 0.75 1]);
plot([20 36],[0.5 0.5],'k--','HandleVisibility','off');
legend('simple','AR','Location','Northeast','Orientation','horizontal');
title('Experiment 2');
axis square;
box on;

