%% Figure 7A
% xaxis as binocular region size instead of monocular FOV/region 
% simple vs AR stimuli, plotting the average and 95% CI, prop fading time
clear all;
close all;

addpath('./data');
%load the first participant data and use its trial order to sort all other
%participant data
load('POL_1.mat');
stimOrder = sortrows(dat.stim);
alldata = stimOrder;

N = 20;
%load each participant's data
for s = 1:N
    load(['POL_',num2str(s),'.mat']);
    % dat.order = whether convergent is first or second
    % add response to this participant's stim, then sort based on P1's trial order
    subjdata = [dat.stim, dat.resp==dat.order];
    reorder_subjdata = sortrows(subjdata); %rows sorted based on stim
    
    %sanity check for trial order
    if reorder_subjdata(:,1:6) ~= stimOrder
        disp('something is wrong')
    end
    
    alldata = [alldata reorder_subjdata(:,7)];

end

propData = []; %stimulus + all participant proportion
r = 1;
for type = [1 3] %stimulus types (Oval,Rect, AR) %change stimtype 2 -->1 as simple shape
    for mono_field = [30 40]
        for mono_zone = [4.5 9]
            if type ==1            
                ind = find((alldata(:,1)==1 | alldata(:,1)==2) & alldata(:,3)==mono_field & alldata(:,4)==mono_zone);
            else 
                ind = find(alldata(:,1)==type & alldata(:,3)==mono_field & alldata(:,4)==mono_zone);
            end
            subdata = alldata(ind,:);
            %prop across repeats for each subj, combine the shapes
            propConv = sum(subdata(:,7:end),1)./length(ind);
            bino_reg = mono_field - mono_zone;
            propData(r,:) = [type,bino_reg,propConv];
            
            r = r+1;
        end
    end
end

%save the proportion data for chi sq. test
save('result.mat','propData');

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
    means = mean(subdata(:,3:end),2);
    stds = std(subdata(:,3:end),[],2,'omitnan'); %the standard deviation of the data, there shouldn't be any NaNs
    z = 1.96;
    n = sum(~isnan(subdata(:,3:end)),2); %Ns should be 20ppl since there's no missing data from this expt.
    if n ~=20
        disp('something is wrong')
    end
    xvals = [subdata(:,2);flip(subdata(:,2))];
    yvals = [means+z*stds./sqrt(n); flip(means-z*stds./sqrt(n))];
    %mean + error bar
    fill(xvals,yvals,fillc(c,:),'LineStyle','none','HandleVisibility','off','FaceAlpha',0.05);
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

