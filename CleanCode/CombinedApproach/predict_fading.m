%% Figure 12B

%use Eq. 8 to predict the amount of fading at various distances for
%the three configurations (symmetrical, nasal shift, temporal shift)

clear all;
close all;
IPD = 6; %in cm
Mf = 40; %each eye's field size, in deg

n = 1;
for D = 0.1:0.01:4 %diopters
    d = (1/D)*100; %distance in cm
    for asym_size = [-10 0 10] %amount of asymmetry between left and right half of Mf, nasal angle - temp angle)
        %negative = temporal shift, positive = nasal shift
        if asym_size == 0
            nasal_half = Mf/2;
            temp_half = Mf/2;
        else %nasal - temp > 0 --> nasal shift; nasal - temp < 0 --> temporal shift
            nasal_half = Mf/2+asym_size/2;
            temp_half = Mf/2-asym_size/2;
        end
        
        %use the x coordinates of Ll Lr Rl Rr
        %binocular = nasal half Lr - Rl (if postive, there's binocular)
        %total FOV = temporal half Rr - Ll
        Lr = d * tand(nasal_half) -IPD/2;
        Rl = -d * tand(nasal_half) +IPD/2;
        
        Ll = -d * tand(temp_half) -IPD/2;
        Rr = d * tand(temp_half) +IPD/2;
        
        max_edge = max([Lr Rr]);
        min_edge = min([Lr Rr]);
        %finding the angles
        T_fov = 2*atand(max_edge/d);
        B_reg = 2*atand(min_edge/d); %if negative, no binocular overlap
        propBino = B_reg / T_fov;
        
        u = 90-nasal_half;
        nearpt = tand(u) * (IPD/2);
        M_reg = (T_fov - B_reg)/2;
        Mratio = M_reg/Mf; %monocular region to monocular field ratio
        Fade = 1.405*sqrt(Mratio);
        
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

f = figure('Units', 'centimeters', 'Position', [0.1, 3, 15, 15], 'PaperPositionMode','Auto');
subplot(1,2,1);
hAx(1) = gca;
hold on;
plot(data_sym(:,1),data_sym(:,5),'-','Color',[154 80 159]./255,'LineWidth',2);
plot(data_conv(:,1),data_conv(:,5),'-','Color',[76 175 74]./255,'LineWidth',2);
plot(data_div(:,1),data_div(:,5),'-','Color',[76 175 74]./255*0.5,'LineWidth',2);
plot([0 4],[Mf,Mf],'k--','HandleVisibility','off'); 
legend('sym','conv','div','Location','Best');
xlabel('diopters (1/meter)');
ylim([0 1.2]);
axis square;

ylabel('prop fading');
yticks([0 0.25 0.5 0.75 1]);

xlabel('distance(D)');
set(hAx(1), 'Xdir', 'reverse')


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
