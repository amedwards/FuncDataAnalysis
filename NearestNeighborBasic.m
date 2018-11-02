function percent_in_cat = NearestNeighborBasic(C,idx,daytime,yaxisname,category,group_names,daysbefore,daysafter)

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
for n = 1:neighborhoods % Neighbor Number
    bn = plotorder(n);
    subplot(2,3,n)
    
    if size(C,2)==1
        linetoplot=ones(length(daytime),1)*C(bn);
        plot(daytime,linetoplot,'LineWidth',3);
        xlabel('Days Until Event')
    elseif size(C,2)==2
        linetoplot=(daytime-daysafter-daysbefore)*C(bn,1)+C(bn,2);
        plot(daytime,linetoplot,'LineWidth',3);
        xlabel('Days Until Event')
    else
        linetoplot=C(bn,:);
        plot(linetoplot,'LineWidth',3);
        xlabel('Coefficients')
    end

    
    hold on
    
    % Find number of babies in each category
    babies_in_hood = sum(idx==bn); % babies in neighborhood
    for c=1:categories
        babies_in_cat(bn,c) = sum(and((idx==bn),category==c)); % babies in category
        percent_in_cat(bn,c) = round(babies_in_cat(bn,c)/babies_in_hood*100);
    end


    ylabel(yaxisname)
    PC_title_string = '';
    for pc = 1:size(C,2)
        PC_title_string = horzcat(PC_title_string,[', PC' num2str(pc) ': ' num2str(round(C(bn,1),2))]);
    end
    title({['Neighbor: ' num2str(bn) PC_title_string]...
        [num2str(babies_in_hood) ' out of ' num2str(length(idx)) ' babies are in this neighborhood. '],...
        [num2str(percent_in_cat(bn,1)) '%: ' group_names(1,:) ', ' num2str(100 - percent_in_cat(bn,1)) '%: ' group_names(2,:)]})
    if size(C,2)<=2
        legend('Group Mean','Location','northwest')
    elseif size(C,2)==3
        legend('1st Principal Component','2nd Principal Component','3rd Principal Component','Group Mean','Sum of Fit','Location','northwest')
    elseif size(C,2)==4
        legend('1st Principal Component','2nd Principal Component','3rd Principal Component','4th Principal Component','Group Mean','Sum of Fit','Location','northwest')
    end
    ylim([-3 7])
    hold off
end