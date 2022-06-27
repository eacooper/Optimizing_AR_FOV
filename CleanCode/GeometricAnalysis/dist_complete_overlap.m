%% Figure 11B

%for nasal shift with different Mf sizes, what is d0 (distance with
%complete overlap)

clear all;
close all;
IPD = 6; %in cm

d = 100; %in cm, distance to focal plane
conv_color = [76 175 74]./255;

n = 1;
N = 21;
for Mf = [20 40 60]
    for asym_size = 0:1:20 %amount of asymmetry between left and right half of Mf, nasal angle - temp angle
        %negative = temporal shift, positive = nasal shift
        if asym_size == 0
            nasal_half = Mf/2;
            temp_half = Mf/2;
        else %nasal - temp > 0 --> nasal shift; nasal - temp < 0 --> temporal shift
            nasal_half = Mf/2+asym_size/2;
            temp_half = Mf/2-asym_size/2;
        end
        
        d = IPD/(tand(nasal_half) - tand(temp_half));
        d = d/100; %convert to m
        data_conv(n,:) = [Mf asym_size 1/d];
        n = n+1;
    end
    
end

f = figure();
subplot(1,2,1);
hold on;
yyaxis left;

plot(data_conv(1:N,2),data_conv(1:N,3),':','Color',conv_color, 'LineWidth', 1.5);
plot(data_conv(N+1:2*N,2),data_conv(N+1:2*N,3),'-','Color',conv_color ,'LineWidth',1);
plot(data_conv(2*N+1:3*N,2),data_conv(2*N+1:3*N,3),'--','Color',conv_color ,'LineWidth',1.5);

legend('f_m=20','f_m=40','f_m=60','Location','Best');
xlabel('asymmetry (deg)');
ylabel('distance(D)');
ylim([0 4]);
set(gca, 'Ydir', 'reverse')

yyaxis right;
ylabel('distance(cm)');
ylim([0 4]);
ylist = 0:8;
yticks(ylist);
yticklabels({'infinity','100','50','','25'});
axis square;
set(gca, 'Ydir', 'reverse')
title('Overlap Distances');


