%% Figure 12A
%proportion fading as a function of proportion monocular in Experiment 1
clear all;
close all;

% load data from Exp 1, variable is slist
load('avgdata.mat');
%	- col 1 = stimulus type (1=simple, 3 = AR)
%	- col 2 = stimulus configuration (1 = convergent, 2 = divergent)
%	- col 3 = binocular overlap size
%	- col 4 = total FOV
%	- col 5 = percentage of monocular region/monocular field

%average the convergent and divergent
new = [slist(9:12,8),slist(13:16,8)];
new = mean(new,2);
new = [slist(9:12,5) new];

c = [51 127 186]./255; %marker color
%plot the data points
data = [new; [0 0]]; %add in a theoretical point for [0 0], 0 monocular -> 0 fading
x = data(:,1)/100; %turn percentage into proportions
y = data(:,2);
%fit a square root
p = polyfit(sqrt(x),y,1);
x2 = [0:0.01:50.65]/100; %piecewise1, turn percentage into proportions 
y2 = p(1).*sqrt(x2); %Eq. 7
x3 = [50.65:1:100]/100; %%piecewise2, turn percentage into proportions 
y3 = ones(size(x3)); %Eq. 7
%% make the plot
figure();
subplot(1,2,1);
hold on;
%plot the data
plot(x,y,'bo','Color',c,'MarkerFaceColor',c);
%plot the piecewise function
plot(x2,y2,'-','Color',[0.5 0.5 0.2],'MarkerFaceColor','y');
plot(x3,y3,'-','Color',[0.5 0.5 0.2],'MarkerFaceColor','y');
ylim([0 1.1]);
xlim([0 1]);
ylabel('prop fading');
xlabel('proportion ratio');
axis square;



