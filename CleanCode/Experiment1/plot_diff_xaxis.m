%% Figure 6

% plot prop time fading against 3 independent variable (x axis): bino
% region, total FOV, and proportion monocular

clear all;
close all;
addpath('./data');
files = dir(fullfile('./data/*.mat'));
files = {files.name};
%stim info:[stimtype configuration monocular_field monocular_region scene icon];

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
%calculate the proportion, if there's too many invalid sample points (based on tres), throw away that trial (make into
%NaN and do calculation by omitting nans)
tres = 150;
a = [];
for n = 1:size(files,2)
    raw_resp = alldata{n};
    for trial = 1:size(raw_resp,1)
        
        see_lune = find(raw_resp(trial,:) ==1);
        
        no_lune = find(raw_resp(trial,:) ==0);
        
        others = sum(raw_resp(trial,:) ==3) + sum(raw_resp(trial,:)==2) + sum(double(isnan(raw_resp(trial,:)))); %number of times no key response
        
        nume = size(see_lune,2);
        denom = (size(see_lune,2)+size(no_lune,2));

        prop = nume / denom;
        a = [a ; nume denom prop];
        numobs = denom;
        checkint = floor(numobs)~=ceil(numobs);
        

        nokey_result(trial,n) = others;
        
        if others<tres
            result(trial,n) = prop; %each column is a cond,row = subj    alldata_prop = ;
            result_count(trial,n) = numobs;
        else
            result(trial,n) = NaN;
            result_count(trial,n) = NaN;
        end
    end
    
end

%combining stim and data, and average across the different types of stims
%for each subj
stim_data = [stimOrder result result_count];

stim_data = sortrows(stim_data,[1 2 3]);

%stim list without repeats and average for each subj across repeats
%(single trial for Oval and Rect, 5 scenes for AR)
colors = {'b.','r.'}; %convergent=blue or divergent=red
cl = {'b','r'};
offset = [0.1 -0.1];
slist = []; %just bino for plotting
slist_mfmz = []; %full set for R

for st = [1 3]
    for cf = [1 2]
        for mf = [30 40]
            for mz = [9 4.5]
                %combine the 2 shapes
                if st ==1
                    subdata = stim_data((stim_data(:,1)==1 | stim_data(:,1)==2) & stim_data(:,2)==cf &stim_data(:,3)==mf & stim_data(:,4)==mz,:);
                    subdata(:,1) = ones(size(subdata,1),1);
                 
                else
                    subdata = stim_data(stim_data(:,1)==st & stim_data(:,2)==cf &stim_data(:,3)==mf & stim_data(:,4)==mz,:);
                end
                
                %average across similar trials (i.e different AR scenes)
                avgdata = mean(subdata,1,'omitnan');
                sumdata = sum(subdata,1,'omitnan');
                
                avgprop = mean(avgdata(7:26),'omitnan');

                bino_size = mf-mz;
                total_size = mf+mz;
                ratio_prt = 100*mz/mf;
                slist = [slist; st,cf,bino_size,total_size,ratio_prt,...
                                mf,mz,avgprop,avgdata(7:26)];
                

            end
        end
    end
end

save('../CombinedApproach/avgdata.mat','slist');

%% MAKE THE PLOT
fig = figure();
markers = {'r^-','ro--','b^-','bo--'};
colors = [228 30 38;0.5*[228 30 38];51 127 186 ;0.5*[51 127 186]];
markers2 = {'r^','ro','b^','bo'};
fillc = colors./255;
offset = [2 4 6 8]';
firstind = 9;
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

%fit through convergent and divergent data.
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
y = -0.017458.*x+1.106;
subplot(2,2,1)
hold on;
plot(x,y,'p-','HandleVisibility','off');
title('adjR2=0.542'); %bino region

x = 34.5:1:50;
y = -0.0028631.*x+0.73231;
subplot(2,2,2)
hold on;
plot(x,y,'p-','HandleVisibility','off');
title('adjR2=-0.148'); %total FOV

x = 10:1:30;
y = 0.015709.*x+0.3035;
subplot(2,2,3)
hold on;
plot(x,y,'p-');
title('adjR2=0.823'); %ratio
legend('conv','div','regfit','Location','South');

