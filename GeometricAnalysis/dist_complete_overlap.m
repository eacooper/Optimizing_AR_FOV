%% Figure 11B

%for nasal shifted frusta with different monocular field (Mf) sizes, what is d0 (distance with
%complete overlap) Eq. 6
clear all;
close all;

%fixed parameters 
IPD = 6; %interpupillary distance in cm 
d = 100; %distance to focal plane in cm

n = 1; %index for storing calculated result

for Mf = [20 40 60] %monocular field sizes

    %asym_size refers to the amount of asymmetry between left and right half of Mf (nasal angle - temp angle)
    %0 = symmetric, negative = temporal shift, positive = nasal shift

    for asym_size = 0:1:20 % positive = nasal shift

        % nasal and temporal halves of the FOV angle
        nasal_half = Mf/2+asym_size/2;
        temp_half = Mf/2-asym_size/2;
        
        d = IPD/(tand(nasal_half) - tand(temp_half));
        d = d/100; %convert to meters
        data_conv(n,:) = [Mf asym_size 1/d];
        n = n+1;
    end
    
end


N = 21; % index for splitting the data into different monocular field sizes
conv_color = [76 175 74]./255; % color

f = figure();
hold on;
%plot for each monocular field size, by varying the amount of nasal shift (asymmetry) what's the distance at which there's 100%
%overlap
plot(data_conv(1:N,2),data_conv(1:N,3),':','Color',conv_color, 'LineWidth', 1.5);
plot(data_conv(N+1:2*N,2),data_conv(N+1:2*N,3),'-','Color',conv_color ,'LineWidth',1);
plot(data_conv(2*N+1:3*N,2),data_conv(2*N+1:3*N,3),'--','Color',conv_color ,'LineWidth',1.5);

legend('f_m=20','f_m=40','f_m=60','Location','Best');
xlabel('asymmetry (deg)');
ylabel('distance(D)');
ylim([0 4]);
axis square;
set(gca, 'Ydir', 'reverse')
title('Overlap Distances');


