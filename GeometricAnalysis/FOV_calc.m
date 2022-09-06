%% Figure 11A, Eq. 1-5
% plot the amount of total FOV size and binocular region size for the given
% camera frusta and content distance
clear all;
close all;

%fixed parameters 
IPD = 6; %in cm
Mf = 40; %each eye's monocular field size, in deg

n = 1; %index for storing calculated result

%% do the calculations
for D = 0.1:0.1:4 %distance in diopters
    
    d = (1/D)*100; %distance in cm

    %asym_size refers to the amount of asymmetry between left and right half of Mf (nasal angle - temp angle)
    %0 = symmetric, negative = temporal shift, positive = nasal shift
    
    for asym_size = [-10 0 10] 

        % nasal and temporal halves of the FOV angle
        nasal_half = Mf/2+asym_size/2;
        temp_half = Mf/2-asym_size/2;
        
        %use the x coordinates of Ll Lr Rl Rr (see Figure 10)
        Lr = d * tand(nasal_half) -IPD/2;
        Rl = -d * tand(nasal_half) +IPD/2;
        
        Ll = -d * tand(temp_half) -IPD/2;
        Rr = d * tand(temp_half) +IPD/2;
        
        %the binocular overlap will always be <= total FOV, so it's bounded
        %by the edges with the smaller x coordinate whereas total FOV is bounded by
        %the edges with the larger x coordinate
        max_edge = max([Lr Rr]);
        min_edge = min([Lr Rr]);

        %finding the angles for total FOV and binocular region
        T_fov = 2*atand(max_edge/d);
        B_reg = 2*atand(min_edge/d);
        
        %proportion of the total FOV that's binocular
        propBino = B_reg / T_fov;
        
        if asym_size == 0  %symmetrical
            data_sym(n,:) = [D d T_fov B_reg propBino];
        elseif asym_size>0 %nasal - temp > 0 --> nasal shift
            data_conv(n,:) = [D d T_fov B_reg propBino];
        elseif asym_size<0 %nasal - temp < 0 --> temporal shift
            data_div(n,:) = [D d T_fov B_reg propBino];
        end

    end
    
    n = n+1;
end

%% make the plot
f = figure('Units', 'centimeters', 'Position', [0.1, 3, 15, 15], 'PaperPositionMode','Auto');

hAx(1) = gca;
hold on;
%plot total FOV
plot(data_sym(:,1),data_sym(:,3),'-','Color',[154 80 159]./255,'LineWidth',2);
plot(data_conv(:,1),data_conv(:,3),'-','Color',[76 175 74]./255,'LineWidth',2);
plot(data_div(:,1),data_div(:,3),'-','Color',[76 175 74]./255*0.5,'LineWidth',2);
plot([0 4],[Mf,Mf],'k--','HandleVisibility','off'); 
%plot binocular region
plot(data_sym(:,1),data_sym(:,4),'--','Color',[154 80 159]./255,'LineWidth',2);
plot(data_conv(:,1),data_conv(:,4),'--','Color',[76 175 74]./255,'LineWidth',2);
plot(data_div(:,1),data_div(:,4),'--','Color',[76 175 74 ]./255*0.5,'LineWidth',2);

ylim([0 80]);
axis square;
ylabel('angular size (deg)');
yticks([0 20 40 60 80]);
legend('sym f_c','conv f_c','div f_c','sym r_b','conv r_b','div r_b','Location','Best','NumColumns',2);
xlabel('distance(D)');
set(hAx(1), 'Xdir', 'reverse')

% add another axis for distance in cm
hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
set(hAx(2),'ytick',[]);
xlim([0,4]);
xlist = [0 1 2 3 4];
xticks(xlist);
xticklabels({'infinity','100','50','33','25'})
xlabel('distance(cm)');
set(hAx(2), 'Xdir', 'reverse')
title('Different Frustum Configuration');
axis square;





