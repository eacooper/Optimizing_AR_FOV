%% Figure 12B

%use geometric analysis predict the proportion monocular at various distances for
%the three configurations (symmetrical, nasal shift, temporal shift)
%then use Eq 8 to predict the amount of fading for different configurations

clear all;
close all;

%fixed parameters 
IPD = 6; %interpupillary distance in cm
Mf = 40; %each eye's field size, in deg

n = 1; %index for storing calculated result

for D = 0.1:0.01:4 %diopters
    d = (1/D)*100; %distance in cm

    for asym_size = [-10 0 10] %amount of asymmetry between left and right half of Mf, nasal angle - temp angle)

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
        B_reg = 2*atand(min_edge/d); %if negative, no binocular overlap
        propBino = B_reg / T_fov;
        
        M_reg = (T_fov - B_reg)/2;
        Mratio = M_reg/Mf; %monocular region to monocular field ratio

        %calculate the amount of fading based on the proportion of
        %monocular region
        if Mratio > 0.5
            Fade = 1;
        else
            Fade = 1.405*sqrt(Mratio);
        end
        
        %set boundaries
        if Fade >1
            Fade =1;
        elseif Fade <0
            Fade = 0;
        end
        
        if asym_size == 0
            data_sym(n,:) = [D d T_fov B_reg Fade Mratio];
        elseif asym_size>0 %nasal - temp > 0 --> nasal shift
            data_conv(n,:) = [D d T_fov B_reg Fade Mratio];
        elseif asym_size<0 %nasal - temp < 0 --> temporal shift
            data_div(n,:) = [D d T_fov B_reg Fade Mratio];
        end

    end
    
    n = n+1;
end

%% make the plot
f = figure('Units', 'centimeters', 'Position', [0.1, 3, 15, 15], 'PaperPositionMode','Auto');
hAx(1) = gca;
hold on;
%plot the different configurations
plot(data_sym(:,1),data_sym(:,5),'-','Color',[154 80 159]./255,'LineWidth',2);
plot(data_conv(:,1),data_conv(:,5),'-','Color',[76 175 74]./255,'LineWidth',2);
plot(data_div(:,1),data_div(:,5),'-','Color',[76 175 74]./255*0.5,'LineWidth',2);
%add a reference line for the monocular field size
plot([0 4],[Mf,Mf],'k--','HandleVisibility','off');

legend('sym','conv','div','Location','Best');
xlabel('diopters (1/meter)');
ylim([0 1.2]);
axis square;
ylabel('prop fading');
yticks([0 0.25 0.5 0.75 1]);
xlabel('distance(D)');
set(hAx(1), 'Xdir', 'reverse')
%add another set of axis for distance in cm
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
