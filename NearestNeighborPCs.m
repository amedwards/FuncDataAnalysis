function percent_in_cat = NearestNeighborPCs(C,idx,pc_fdmat,meanfd_fdmat,daytime,yaxisname,category,group_names_orstruct,daysbefore,daysafter,limofinterest,dayfd,id,color)

% load('gooddata_hero.mat')
% load('kmeans4_fromharmscr.mat')
% load('pc_fd_hero.mat')
% pc_fdmat = fdmat;
% load('meanfd_fdmat_hero.mat')
% meanfd_fdmat = fdmat;
neighborhoods = length(unique(idx));
categories = length(unique(category));
numincategory = zeros(neighborhoods,1);
idsincategory = zeros(neighborhoods,1);

if isstruct(group_names_orstruct)
    numgroups = size(group_names_orstruct,2);
    for g=1:numgroups
        categories(g) = size(group_names_orstruct(g).group_names,1);
    end
else
    numgroups = 1;
end

percent_in_cat = zeros(neighborhoods,max(categories),numgroups);
babypercent_in_cat = zeros(neighborhoods,max(categories),numgroups);
cultures_in_cat = zeros(neighborhoods,max(categories),numgroups);
ids_in_cat = zeros(neighborhoods,max(categories),numgroups);
mortperc_in_cat = zeros(neighborhoods,max(categories),numgroups);

for n=1:neighborhoods
    numincategory(n) = sum(idx==n);
%     idsincategory(n) = length(unique(id(idx==n)));
end
[~,plotorder] = sort(numincategory,'descend');
plotorder = sort(plotorder);
for p = 1:neighborhoods % Neighbor Number
    bn = plotorder(p);
    if neighborhoods<=6
        subplot(2,3,p)
    else
        subplot(3,4,p)
    end
    pc_xaxis = linspace(daysbefore,daysafter,length(pc_fdmat));
    meanfd_xaxis = linspace(daysbefore,daysafter,length(meanfd_fdmat));

    % Resample the fd objects so that they line up with the data points
    meanfd_fdmat_resampled = interp1(meanfd_xaxis,meanfd_fdmat,daytime);
    pc_fdmat_resampled = interp1(pc_xaxis,pc_fdmat,daytime);
    
    % Resample the individual baby objects so that they line up with the
    % data points
    if ~isempty(dayfd)
        ind_baby_fdmat_resampled = eval_fd(daytime,dayfd,int2Lfd(0));
    end

    % Plot the principal components multiplied by the coefficients for a single
    % baby
%     plot(daytime,C(bn,1:size(pc_fdmat_resampled,2)).*pc_fdmat_resampled);
    hold on
    % Plot the mean
    multiplot = size(yaxisname,1)>1;
    
    colors = {[1 0 0];[0 1 0]; [0 0 1]};
    
    if max(meanfd_fdmat_resampled)<1
        multiplier = 1000;
    else
        multiplier = 1;
    end
   
    if multiplot
        h = plot(daytime,meanfd_fdmat_resampled*multiplier,'LineWidth',1,'Color',color);
    else
        plot(daytime,meanfd_fdmat_resampled*multiplier,'LineWidth',1,'Color',[0.6 0.6 0.6]);
    end

    % Add the mean to the principal component sum
    pcsum = sum(C(bn,1:size(pc_fdmat_resampled,2)).*pc_fdmat_resampled,2);
    if isrow(pcsum)
        pcsum = pcsum';
    end
    if isrow(meanfd_fdmat_resampled)
        meanfd_fdmat_resampled = meanfd_fdmat_resampled';
    end
    pc_plus_mean = pcsum + meanfd_fdmat_resampled;

    % Plot all the babies in the neighborhood
%     if ~isempty(dayfd)
%         plot(daytime,ind_baby_fdmat_resampled(:,idx==bn));
%     end
    hold on
    % Plot the sum of the principal components plus the mean for a particular
    % neighborhood
    if multiplot
        j = plot(daytime,pc_plus_mean*multiplier,'LineWidth',3,'Color',color);
    else
        plot(daytime,pc_plus_mean*multiplier,'LineWidth',3,'Color',[0 0 0]);
    end
    
    % Find number of babies in each category
    cultures_in_hood = sum(idx==bn); % cultures in neighborhood
%     babies_in_hood = length(unique(id(idx==bn)));
    for g=1:numgroups
        for c=1:categories(g)
%             incategory = and((idx==bn),group_names_orstruct(g).category==c);
%             ids_in_cat(bn,c,g) = length(unique(id(incategory)));
            cultures_in_cat(bn,c,g) = sum(and((idx==bn),group_names_orstruct(g).category==c)); % cultures in category
            percent_in_cat(bn,c,g) = round(cultures_in_cat(bn,c,g)/cultures_in_hood*100);
%             babypercent_in_cat(bn,c,g) = round(cultures_in_cat(bn,c,g)/babies_in_hood*100);
%             mortperc_in_cat(bn,c,g) = round(ids_in_cat(bn,c,g)/babies_in_hood*100);
        end
    end
    
    ax = gca;
    ax.LineWidth = 2;
    ax.Box = 'off';
    ax.FontSize = 14;
    ax.FontName = 'Arial';
    xlabel('Days to Event')
    if multiplot
        if size(yaxisname,1)==2
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2}];
        elseif size(yaxisname,1)==3
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2} '/' yaxisname{3}];
        elseif size(yaxisname,1)==4
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2} '/' yaxisname{3} '/' yaxisname{4}];
        elseif size(yaxisname,1)==5
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2} '/' yaxisname{3} '/' yaxisname{4} '/' yaxisname{5}];
        elseif size(yaxisname,1)==6
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2} '/' yaxisname{3} '/' yaxisname{4} '/' yaxisname{5} '/' yaxisname{6}];
        elseif size(yaxisname,1)==7
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2} '/' yaxisname{3} '/' yaxisname{4} '/' yaxisname{5} '/' yaxisname{6} '/' yaxisname{7}];
        elseif size(yaxisname,1)==8
            yaxisnamenew = [yaxisname{1} '/' yaxisname{2} '/' yaxisname{3} '/' yaxisname{4} '/' yaxisname{5} '/' yaxisname{6} '/' yaxisname{7} '/' yaxisname{8} ];
        end
    else
        yaxisnamenew = yaxisname;
    end
    ylabel(yaxisnamenew)
    PC_title_string = '';
%     for pc = 1:size(C,2)
%         PC_title_string = horzcat(PC_title_string,[', PC' num2str(pc) ': ' num2str(round(C(bn,pc),2))]);
%     end

    if ~isstruct(group_names_orstruct)
        title({['Cluster ' num2str(bn) PC_title_string]...
            [num2str(cultures_in_hood) ' of ' num2str(length(idx)) ' cultures'],...
            [num2str(percent_in_cat(bn,1)) '% ' strtrim(group_names_orstruct(1,:))]})
    else
        if numgroups==1
            title({['Cluster ' num2str(bn) PC_title_string]...
            [num2str(cultures_in_hood) ' of ' num2str(length(idx)) ' cultures'],...
            [num2str(percent_in_cat(bn,1,g)) '% ' strtrim(group_names_orstruct(g).group_names(2,:))]})
        elseif numgroups==2
            title({['Cluster ' num2str(bn) PC_title_string ': ' num2str(cultures_in_hood) ' of ' num2str(length(idx)) ' cultures'],...
            [num2str(percent_in_cat(bn,1,1)) '% ' strtrim(group_names_orstruct(1).group_names(2,:))],...
            [num2str(percent_in_cat(bn,1,2)) '% ' strtrim(group_names_orstruct(2).group_names(2,:))]})
        elseif numgroups==3
            title({['Cluster ' num2str(bn) PC_title_string ': ' num2str(cultures_in_hood) ' of ' num2str(length(idx)) ' cultures'],...
            [num2str(percent_in_cat(bn,1,1)) '% ' strtrim(group_names_orstruct(1).group_names(2,:))],...
            [num2str(percent_in_cat(bn,1,2)) '% ' strtrim(group_names_orstruct(2).group_names(2,:))],...
            [num2str(percent_in_cat(bn,1,3)) '% ' strtrim(group_names_orstruct(3).group_names(2,:))]})
        elseif numgroups==4
            title({['Cluster ' num2str(bn) PC_title_string ': ' num2str(cultures_in_hood) ' of ' num2str(length(idx)) ' cultures'],...
            [num2str(percent_in_cat(bn,1,1)) '% ' strtrim(group_names_orstruct(1).group_names(2,:)) ', ' num2str(percent_in_cat(bn,1,2)) '% ' strtrim(group_names_orstruct(2).group_names(2,:))],...
            [num2str(percent_in_cat(bn,1,3)) '% ' strtrim(group_names_orstruct(3).group_names(2,:)) ', ' num2str(percent_in_cat(bn,1,4)) '% ' strtrim(group_names_orstruct(4).group_names(2,:))]})
        end
    end


%     if size(C,2)==1
%         legend('1st Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
%     elseif size(C,2)==2
%         legend('1st Principal Component','2nd Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
%     elseif size(C,2)==3
%         legend('1st Principal Component','2nd Principal Component','3rd Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
%     elseif size(C,2)==4
%         legend('1st Principal Component','2nd Principal Component','3rd Principal Component','4th Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
%     end

    ylim(limofinterest)
%     hold off
end

if multiplot
    if size(yaxisname,1)==2
        legend(['Mean Curve ' yaxisname{1}],['Characteristic Curve ' yaxisname{1}],['Mean Curve ' yaxisname{2}],['Characteristic Curve ' yaxisname{2}],'Location','northwest')
    elseif size(yaxisname,1)==3
        legend(['Mean Curve ' yaxisname{1}],['Characteristic Curve ' yaxisname{1}],['Mean Curve ' yaxisname{2}],['Characteristic Curve ' yaxisname{2}],['Mean Curve ' yaxisname{3}],['Characteristic Curve ' yaxisname{3}],'Location','northwest')
    end
else
    legend('Mean Curve','Characteristic Curve ','Location','northwest')
end