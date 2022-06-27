%% Figure 12A

clear all;
close all;
%cols = st,cf,bino_size,total_size,ratio_prt,mf,mz,avgprop,avgdata(7:26)];
load('avgdata.mat');
%average the convergent and divergent
new = [slist(9:12,8),slist(13:16,8)];
new = mean(new,2);
new = [slist(9:12,5) new];

figure();
subplot(1,2,1);
hold on;
c = [51 127 186]./255;
data = [new; [0 0]]; %add in a theoretical point for [0 0], 0 monocular -> 0 fading
x = data(:,1);
y = data(:,2);
plot(x,y,'bo','Color',c,'MarkerFaceColor',c);
ylim([0 1.1]);
xlim([0 100]);

x2 = 0:0.01:50.65;
y2 = 0.1405.*sqrt(x2); %Eq. 7
x3 = 50.65:100;
y3 = ones(size(x3)); %Eq. 7
plot(x2,y2,'-','Color',[0.5 0.5 0.2],'MarkerFaceColor','y');
plot(x3,y3,'-','Color',[0.5 0.5 0.2],'MarkerFaceColor','y');
ylabel('prop fading');
xlabel('ratio (%)');
axis square;


