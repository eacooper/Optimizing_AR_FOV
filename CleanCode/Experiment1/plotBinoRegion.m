%% FIgure 5 and Figure 7B

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

load(files{1});
stimOrder = dat.stim;
alldata = {};

%sort all the subj's data by trials
for s = 1:size(files,2)
    load(files{s});
    subjdata = [dat.stim dat.key_down];
    [~,ind] = ismember(stimOrder,subjdata(:,1:6),'rows');
    reorder_subjdata = subjdata(ind,:);
    alldata = [alldata reorder_subjdata(:,7:end)];
end

result = [];
result_count = [];
nokey_result = [];

%calculate the proportion of time seeing fading
tres = 150; %threshold for invalid key response
a = [];
for n = 1:size(files,2)
    raw_resp = alldata{n};
    for trial = 1:size(raw_resp,1)
        %find key presses that indicate fading
        see_lune = find(raw_resp(trial,:) ==1);
        %find key presses that indicate no fading
        no_lune = find(raw_resp(trial,:) ==0);
        %find key presses that are invalid (neither fading/no fading)
        others = sum(raw_resp(trial,:) ==3) + sum(raw_resp(trial,:)==2) + sum(double(isnan(raw_resp(trial,:)))); %number of times no key response
        %calculate the proportion of time fading / (fading + no fading)
        nume = size(see_lune,2);
        denom = (size(see_lune,2)+size(no_lune,2));
        prop = nume / denom;
        
        a = [a ; nume denom prop];
        numobs = denom; %number of observations (exclude invalid key responses)
        checkint = floor(numobs)~=ceil(numobs);
        

        nokey_result(trial,n) = others;
        
        %if there are too many invalid sample points (based on tres), throw away that trial (make into
        %NaN and do calculation by omitting nans)

        if others<tres
            result(trial,n) = prop; %each column is a cond,row = subj
            result_count(trial,n) = numobs;
        else
            result(trial,n) = NaN;
            result_count(trial,n) = NaN;
        end
    end
    
end

%combining stim info and the processed data
stim_data = [stimOrder result result_count];

stim_data = sortrows(stim_data,[1 2 3]);


%stim list without repeats and average for each subj across repeats
%(single trial for Oval and Rect, 5 scenes for AR)
slist = []; 

for st = [1 3] %stimulus types (oval, rect, AR)
    for cf = [1 2] %display configuration
        for mf = [30 40] %monocular field size 
            for mz = [9 4.5] %monocular zone size
                if st ==1 %combine oval and rect simple stimuli, and make stimtype =1 in subdata
                    subdata = stim_data((stim_data(:,1)==1 | stim_data(:,1)==2) & stim_data(:,2)==cf &stim_data(:,3)==mf & stim_data(:,4)==mz,:);
                    subdata(:,1) = ones(size(subdata,1),1);
                else
                    subdata = stim_data(stim_data(:,1)==st & stim_data(:,2)==cf &stim_data(:,3)==mf & stim_data(:,4)==mz,:);
                end
                
                %average across repeats (exact repeats for simple stimuli and different AR scenes for the AR stimuli)
                avgdata = mean(subdata,1,'omitnan');
                sumdata = sum(subdata,1,'omitnan');
                
                %average across subjects
                avgprop = mean(avgdata(7:26),'omitnan');

                
                %store the data, and convert monocular field and monocular
                %region into binocular overlap size
                slist = [slist; st cf mf-mz avgprop avgdata(7:26)];
                                
            end
        end
    end
end

%% MAKE THE PLOTS
%mean + errorbar
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
        yvals = [subdata(:,4)+z*stds./sqrt(n); flip(subdata(:,4)-z*stds./sqrt(n))];
        
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


figure(2);
hold on;
markers = {'rs-','bd-'};
colors = [[228 30 38]*0.5;51 127 186];
markers2 = {'rs','bd',};
fillc = colors./255;
offset = [2 4 6 8]';
c = 1;
for st = [1 3]
    convdata = slist(slist(:,1)==st&slist(:,2)==1,:);
    divdata = slist(slist(:,1)==st&slist(:,2)==2,:);
    diffdata(:,1:2) = [divdata(:,1) divdata(:,3)]; % stimtype and bino region
    diffdata(:,3:22) = divdata(:,5:end)-convdata(:,5:end); %20 subj data
    diffdata = sortrows(diffdata,2);
    
    %error bar
    stds = std(diffdata(:,3:end),[],2,'omitnan'); %the standard deviation of the data, omitting NaN
    means = mean(diffdata(:,3:end),2,'omitnan');
    z = 1.96;
    n = sum(~isnan(diffdata(:,3:end)),2); %Ns are a bit different since there's some missing data for the condition
    xvals = [diffdata(:,2);flip(diffdata(:,2))];
    yvals = [means+z*stds./sqrt(n); flip(means-z*stds./sqrt(n))];
 
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
