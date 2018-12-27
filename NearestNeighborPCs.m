function percent_in_cat = NearestNeighborPCs(C,idx,pc_fdmat,meanfd_fdmat,daytime,yaxisname,category,group_names,daysbefore,daysafter,limofinterest)

% load('gooddata_hero.mat')
% load('kmeans4_fromharmscr.mat')
% load('pc_fd_hero.mat')
% pc_fdmat = fdmat;
% load('meanfd_fdmat_hero.mat')
% meanfd_fdmat = fdmat;
neighborhoods = length(unique(idx));
categories = length(unique(category));
numincategory = zeros(neighborhoods,1);
percent_in_cat = zeros(neighborhoods,categories);
babies_in_cat = zeros(neighborhoods,categories);
for n=1:neighborhoods
    numincategory(n) = sum(idx==n);
end
[~,plotorder] = sort(numincategory,'descend');
plotorder = sort(plotorder);
for p = 1:neighborhoods % Neighbor Number
    bn = plotorder(p);
    subplot(2,3,p)
    pc_xaxis = linspace(daysbefore,daysafter,length(pc_fdmat));
    meanfd_xaxis = linspace(daysbefore,daysafter,length(meanfd_fdmat));

    % Resample the fd objects so that they line up with the data points
    meanfd_fdmat_resampled = interp1(meanfd_xaxis,meanfd_fdmat,daytime);
    pc_fdmat_resampled = interp1(pc_xaxis,pc_fdmat,daytime);

    % Plot the principal components multiplied by the coefficients for a single
    % baby
    plot(daytime,C(bn,1:size(pc_fdmat_resampled,2)).*pc_fdmat_resampled);
    hold on
    % Plot the mean
    plot(daytime,meanfd_fdmat_resampled);

    % Add the mean to the principal component sum
    pcsum = sum(C(bn,1:size(pc_fdmat_resampled,2)).*pc_fdmat_resampled,2);
    if isrow(pcsum)
        pcsum = pcsum';
    end
    if isrow(meanfd_fdmat_resampled)
        meanfd_fdmat_resampled = meanfd_fdmat_resampled';
    end
    pc_plus_mean = pcsum + meanfd_fdmat_resampled;

    % Plot the sum of the principal components plus the mean for a particular
    % neighborhood
    plot(daytime,pc_plus_mean,'LineWidth',3);
    
    % Find number of babies in each category
    babies_in_hood = sum(idx==bn); % babies in neighborhood
    for c=1:categories
        babies_in_cat(bn,c) = sum(and((idx==bn),category==c)); % babies in category
        percent_in_cat(bn,c) = round(babies_in_cat(bn,c)/babies_in_hood*100);
    end

    xlabel('Days Until Event')
    xlabel('Days of Life')
    ylabel(yaxisname)
    PC_title_string = '';
%     for pc = 1:size(C,2)
%         PC_title_string = horzcat(PC_title_string,[', PC' num2str(pc) ': ' num2str(round(C(bn,pc),2))]);
%     end
    if size(group_names,1)>2
        title({['Neighbor: ' num2str(bn) PC_title_string]...
        [num2str(babies_in_hood) ' out of ' num2str(length(idx)) ' babies are in this neighborhood. '],...
        [num2str(percent_in_cat(bn,1)) '%: ' group_names(1,:) ', ' num2str(percent_in_cat(bn,2)) '%: ' group_names(2,:) ', ' num2str(percent_in_cat(bn,3)) '%: ' group_names(3,:) ]})
    else
        title({['Neighbor: ' num2str(bn) PC_title_string]...
            [num2str(babies_in_hood) ' out of ' num2str(length(idx)) ' babies are in this neighborhood. '],...
            [num2str(percent_in_cat(bn,1)) '%: ' group_names(1,:) ', ' num2str(100 - percent_in_cat(bn,1)) '%: ' group_names(2,:) ]})
    end
    if size(C,2)==1
        legend('1st Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
    elseif size(C,2)==2
        legend('1st Principal Component','2nd Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
    elseif size(C,2)==3
        legend('1st Principal Component','2nd Principal Component','3rd Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
    elseif size(C,2)==4
        legend('1st Principal Component','2nd Principal Component','3rd Principal Component','4th Principal Component','Mean Curve','Characteristic Curve','Location','northwest')
    end
    ylim(limofinterest)
    hold off
end