%% Figure 5, Figure 6 and Figure 7B

% xaxis as binocular region, plotting the average and 95% CI, prop convergent
% chosen as in Figure 5, and compute the difference between convergent and
% divergent as in Figure 7B

clear all;
close all;

%% PROCESSING DATA
%get data files
addpath('./data');
files = dir(fullfile('./data/*.mat'));
files = {files.name};
%stim info  --> [stimtype configuration monocular_field monocular_region scene icon];

% load stimulus order for first subject to use for sorting trials
load(files{1});
stimOrder = dat.stim;

% initialize a cell array for storing the data
alldata = {};

%sort all the subj's data by trials
for s = 1:size(files,2)
    load(files{s});
    subjdata = [dat.stim dat.key_down];
    [~,ind] = ismember(stimOrder,subjdata(:,1:6),'rows'); % get trials ordering that matches subject 1
    reorder_subjdata = subjdata(ind,:); % re-order trials
    alldata = [alldata reorder_subjdata(:,7:end)]; % combine with other data
end

% initialize matrices to store results
result = [];
result_count = [];
trial_by_trial_data = [];

%calculate the proportion of time seeing fading
tres = 150; %threshold for invalid key response (number of sample points)

% for each subject
for n = 1:size(files,2)
    
    raw_resp = alldata{n}; % grab their data
    
    % for each trial
    for trial = 1:size(raw_resp,1)
        %find key presses that indicate fading
        see_lune = find(raw_resp(trial,:) == 1);
        %find key presses that indicate no fading
        no_lune = find(raw_resp(trial,:) == 0);
        %find key presses that are invalid (wrong key, both keys at once, no key)
        others = sum(raw_resp(trial,:) == 3) + sum(raw_resp(trial,:)== 2) + sum(double(isnan(raw_resp(trial,:))));
        %calculate the proportion of time fading / (fading + no fading)
        nume = size(see_lune,2);
        denom = (size(see_lune,2)+size(no_lune,2));
        prop = nume / denom;
        
        % if there are too many invalid sample points (based on tres), throw away that trial (make into
        % NaN and do calculation by omitting nans)
        % otherwise, store the proportion fading and the total response
        % time

        if others<tres
            result(trial,n) = prop; %each column is a cond,row = subj
            result_count(trial,n) = denom;
            trial_by_trial_data = [trial_by_trial_data ;n,stimOrder(trial,1:4),stimOrder(trial,3)-stimOrder(trial,4),nume,denom,prop];
        else
            result(trial,n) = NaN;
            result_count(trial,n) = NaN;
        end
    end
    
end

trial_by_trial_data = sortrows(trial_by_trial_data,[1 2 3 4]);

% save as csv for glme, make each trial into individual rows
expandslist = [];
stimtype = {'oval','rect','AR'};
configtype = {'conv','div'};
subjid = {'S1','S2','S3','S4','S5','S6',...
    'S7','S8','S9','S10','S11','S12',...
    'S13','S14','S15','S16','S17','S18',...
    'S19','S20'};

subjlist = subjid(trial_by_trial_data(:,1))';       % list of subject numbers
stimlist = stimtype(trial_by_trial_data(:,2))';     % stimulus types (oval, rect, AR)
configlist = configtype(trial_by_trial_data(:,3))'; % conv/div configuration
mflist = trial_by_trial_data(:,4);                  % monocular field size
mrlist = trial_by_trial_data(:,5);                  % monocular region size
binolist = trial_by_trial_data(:,6);                % binocular region size
numFade = trial_by_trial_data(:,7);                 % number of fading button presses
numpoint = trial_by_trial_data(:,8);                % number of valid button presses
proplist = trial_by_trial_data(:,9);                % proportion fading

T = table(subjlist,stimlist,configlist,mflist,mrlist,binolist,numFade,numpoint,proplist);
T.Properties.VariableNames = {'Subj','StimType','ConvDiv','Mfov','MonoReg','BinoReg','NumFade','TimePt','PropFade'};
writetable(T,'./trial_by_trial_data_for_GLM.csv','Delimiter',',');


%combining stim info and the processed data
stim_data = [stimOrder result result_count];

%create stim list averaged over repeats for each subj
%we will combine the single trial for Oval (1) and Rect (2), and also average the 5 scenes for AR (3)
slist = []; 
fovlist = []; 

for st = [1 3] %stimulus types (oval, rect, AR) -- we're going to combine stimulus types 1 and 2
    for cf = [1 2] %display configuration (convergent/divergent)
        for mf = [30 40] %monocular field size 
            for mz = [9 4.5] %monocular zone size
                if st ==1 %combine oval and rect simple stimuli, and make stimtype = 1 in subdata
                    subdata = stim_data((stim_data(:,1)==1 | stim_data(:,1)==2) & stim_data(:,2)==cf &stim_data(:,3)==mf & stim_data(:,4)==mz,:);
                    subdata(:,1) = ones(size(subdata,1),1);
                else
                    subdata = stim_data(stim_data(:,1)==st & stim_data(:,2)==cf &stim_data(:,3)==mf & stim_data(:,4)==mz,:);
                end
                
                %average across repeats (exact repeats for simple stimuli and different AR scenes for the AR stimuli)
                avgdata = mean(subdata,1,'omitnan');
                
                %average across subjects
                avgprop = mean(avgdata(7:26),'omitnan');

                %store the data, and convert monocular field and monocular
                %region into binocular overlap size
                slist = [slist; st cf mf-mz avgprop avgdata(7:26)];
                
                % calculate and store 3 ways to summarize FOV for Fig 6
                bino_size = mf-mz;
                total_size = mf+mz;
                ratio_prt = 100*mz/mf;
                
                fovlist = [fovlist; st,cf,bino_size,total_size,ratio_prt,...
                                mf,mz,avgprop,avgdata(7:26)];
                            

                                
            end
        end
    end
end

% save the data with fov calculations for use in the guidelines
save('../PerceptualPredictions/avgdata.mat','fovlist');

%% MAKE THE PLOTS
%mean + errorbar

% figure 5
figure(1);
hold on;
markers = {'r^-','ro--','b^-','bo--'}; %red = simple, blue = AR stim
colors = [228 30 38;0.5*[228 30 38];51 127 186 ;0.5*[51 127 186]];
markers2 = {'r^','ro','b^','bo'};
fillc = colors./255;
offset = [2 4 6 8]';
c = 1;

for st = [1 3] %simple vs. AR
    for config = [1 2] %display configuration
        subdata = slist(slist(:,1)==st&slist(:,2)==config,:);
        subdata = sortrows(subdata,3);

        %each column is a subj. each row is a trial
        stds = std(subdata(:,5:end),[],2,'omitnan'); %the standard deviation of the data, omitting NaN
        z = 1.96;
        n = sum(~isnan(subdata(:,5:end)),2); %Ns may be different betweene conditions since there might be missing data due to invalid responses
        %xval of errorbar
        xvals = [subdata(:,3);flip(subdata(:,3))];
        yvals = [subdata(:,4)+(z*stds./sqrt(n)); flip(subdata(:,4)-(z*stds./sqrt(n)))];
        
        if st ==1
            subplot(1,2,1)
            hold on;
        else
            subplot(1,2,2)
            hold on;
        end
%       errorbar
        fill(xvals,yvals,fillc(c,:),'LineStyle','none','HandleVisibility','off','FaceAlpha',0.05);
%       average data
        plot(subdata(:,3),subdata(:,4),markers{c},'Color',fillc(c,:),'MarkerFaceColor',fillc(c,:));
%       for each bino overlap, plot raw data
        plot(subdata(1,3)*ones(1,20)-0.5+1*rand(1,20),subdata(1,5:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
        plot(subdata(2,3)*ones(1,20)-0.5+1*rand(1,20),subdata(2,5:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
        plot(subdata(3,3)*ones(1,20)-0.5+1*rand(1,20),subdata(3,5:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
        plot(subdata(4,3)*ones(1,20)-0.5+1*rand(1,20),subdata(4,5:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');

        c = c+1;
    end
    if st ==1
        title('Simple');
        legend('conv','div','Location','Northeast','Orientation','horizontal');
    else
        title('AR');
        legend('conv','div','Location','Northeast','Orientation','horizontal');
    end
    xlabel('Binocular region (deg)');
    ylabel('Prop time fading');
    ylim([0 1.25]);
    xlim([20 36.5]);
    yticks([0 0.25 0.5 0.75 1]);
    axis square;
    box on;

end

% this is the re-plotting of the data for comparison to Exp 2 (figure 7B)
figure(2);
hold on;
markers = {'rs-','bd-'};
colors = [[228 30 38]*0.5;51 127 186];
markers2 = {'rs','bd',};
fillc = colors./255;
offset = [2 4 6 8]';
c = 1;

for st = [1 3]
    convdata = slist(slist(:,1)==st&slist(:,2)==1,:); % grab data just for convergent trials
    divdata = slist(slist(:,1)==st&slist(:,2)==2,:);    % grad data just for divergence trials
    
    diffdata(:,1:2) = [divdata(:,1) divdata(:,3)]; % stimtype and bino region
    diffdata(:,3:22) = divdata(:,5:end)-convdata(:,5:end); % difference between convergent and divergent for each of the 20 subj
    diffdata = sortrows(diffdata,2);
    
    %error bar
    stds = std(diffdata(:,3:end),[],2,'omitnan'); %the standard deviation of the data, omitting NaN
    means = mean(diffdata(:,3:end),2,'omitnan');
    z = 1.96;
    n = sum(~isnan(diffdata(:,3:end)),2); %Ns are a bit different since there's some missing data for the condition
    xvals = [diffdata(:,2);flip(diffdata(:,2))];
    yvals = [means+(z*stds./sqrt(n)); flip(means-(z*stds./sqrt(n)))];
 
    fill(xvals,yvals,fillc(c,:),'LineStyle','none','HandleVisibility','off','FaceAlpha',0.05);
    plot(diffdata(:,2),means,markers{c},'Color',fillc(c,:),'MarkerFaceColor',fillc(c,:));
    % for each bino overlap, plot raw data
    plot(diffdata(1,2)*ones(1,20)-0.5+1*rand(1,20),diffdata(1,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    plot(diffdata(2,2)*ones(1,20)-0.5+1*rand(1,20),diffdata(2,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    plot(diffdata(3,2)*ones(1,20)-0.5+1*rand(1,20),diffdata(3,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    plot(diffdata(4,2)*ones(1,20)-0.5+1*rand(1,20),diffdata(4,3:end),markers2{c},'MarkerEdgeColor',fillc(c,:),'MarkerFaceColor',fillc(c,:),'MarkerSize',2,'HandleVisibility','off');
    
    c = c+1;
    

    xlabel('Binocular region (deg)');
    ylabel('Div - conv prop time fading');
    ylim([-1 1.5]);
    xlim([20 36.5]);
    yticks([-1 -0.5 0 0.5 1]);
    axis square;
    box on;
end
plot([20 36],[0 0],'k--','HandleVisibility','off');

legend('simple','AR','Location','Northeast','Orientation','horizontal');
title('Experiment 1:Div-Conv');


% figure 6

%% MAKE THE PLOT
fig = figure();
markers = {'r^-','ro--','b^-','bo--'};
colors = [228 30 38;0.5*[228 30 38];51 127 186 ;0.5*[51 127 186]];
markers2 = {'r^','ro','b^','bo'};
fillc = colors./255;
offset = [2 4 6 8]';
firstind = 9;

slist = fovlist;

for p = 1:3 %for the three types of xaxis measure (bino region, total fov, prop mono)
    c = 3;
    d = p+2; %average data index in slist
    for st = 3% just the AR stim
        
        for config = [1 2]
            
            subdata = slist(slist(:,1)==st&slist(:,2)==config,:);
            subdata = sortrows(subdata,d);
            
            stds = std(subdata(:,firstind:end),[],2,'omitnan'); 
            z = 1.96;
            n = sum(~isnan(subdata(:,firstind:end)),2); 
            xvals = [subdata(:,d);flip(subdata(:,d))];
            yvals = [subdata(:,firstind-1)+z*stds./sqrt(n); flip(subdata(:,firstind-1)-z*stds./sqrt(n))];

            if st ==3 && p ==1
                subplot(2,2,1)
                hold on;
            elseif st ==3 && p ==2
                subplot(2,2,2)
                hold on;
            elseif st ==3 && p ==3
                subplot(2,2,3)
                hold on;
            end
            fill(xvals,yvals,fillc(c,:),'LineStyle','none','HandleVisibility','off','FaceAlpha',0.05);
            plot(subdata(:,d),subdata(:,firstind-1),markers{c},'Color',fillc(c,:),'MarkerFaceColor',fillc(c,:));
            c = c+1;
        end
        if p ==3
            legend('conv','div','Location','South');
        end
        

        if p ==1
           xlabel('Binocular region (deg)');
        elseif p ==2
            xlabel('Total FOV (deg)');
        elseif p ==3
            xlabel('Ratio mz/mf (%)');
        end
        ylabel('Prop time fading');
        ylim([0 1.25]);
        yticks([0 0.25 0.5 0.75 1]);
        axis square;
        box on;
    end
end

%% fit a line to convergent and divergent data.
%bino_size
data = slist(slist(:,1)==3,[3 8]);
m1= fitlm(data(:,1),data(:,2));
%totalFOV
data = slist(slist(:,1)==3,[4 8]);
m2= fitlm(data(:,1),data(:,2));
%ratio
data = slist(slist(:,1)==3,[5 8]);
m3 = fitlm(data(:,1),data(:,2));

%add in the regression line obtained from fitlm above
x = 20:1:35.5;
y = m1.Coefficients{2,1}.*x+1.106;
subplot(2,2,1)
hold on;
plot(x,y,'-','HandleVisibility','off');
title(['R2=',num2str(m1.Rsquared.Ordinary)]); 

x = 34.5:1:50;
y = m2.Coefficients{2,1}.*x+0.73231;
subplot(2,2,2)
hold on;
plot(x,y,'-','HandleVisibility','off');
title(['R2=',num2str(m2.Rsquared.Ordinary)]); 

x = 10:1:30;
y = m3.Coefficients{2,1}.*x+0.3035;
subplot(2,2,3)
hold on;
plot(x,y,'-');
title(['R2=',num2str(m3.Rsquared.Ordinary)]); 
legend('conv','div','regfit','Location','South');
