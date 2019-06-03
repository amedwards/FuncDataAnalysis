
% Load Working Dataset
% load('X:\Amanda\BSI\OneDayBeforeAndAfterCulture.mat')
vname = {'Positive','GramPositive','GramNegative','Fungus','HR','RespRate','SPO2','ECGDerivedRespv2','MeanRRInterval','StdofRRInt','DiastolicBP','SystolicBP','RelRiskPos','RelRiskGP','RelRiskGN','RelRiskFungal','PulsePressure'};
dataday(:,17,:) = dataday(:,12,:)-dataday(:,11,:);

% Choose Parameters of Interest
varofinterest = {'StdofRRInt','DiastolicBP','PulsePressure'}; %{'HR','ECGDerivedRespv2','MeanRRInterval'};

% Count Data Quantities
[varlog,varnums] = ismember(vname,varofinterest);
nv = sum(varlog);
n = size(dataday,3);

% Outcome Classification
pg = nansum(squeeze(dataday(:,1,:)),1); % Positive Culture
pg = double(pg>0)'; % Negative Culture
gramneg = nansum(squeeze(dataday(:,3,:)),1);
gramneg = double(gramneg>0)';
fungus = nansum(squeeze(dataday(:,4,:)),1);
fungus = double(fungus>0)';

% Deal with Time
tt = -1:1/4320:1; % 20 second chunks
daysbefore = -1; %-0.0417; % 1 hr before culture
daysafter = 0; % 0.0417; % 1 hr after culture
start_tt = find(round(tt,4)==daysbefore);
end_tt = find(round(tt,4)==daysafter);
ttwindow = tt(start_tt:end_tt);
daytime = ttwindow;

% Initialize Empty Arrays
tempdata = zeros(n,nv);
% goodindices = zeros(n,nv);
% percentavail = zeros(n,nv);

for v=1:nv
    % --------------------- Data Preparation ----------------------------
    column = find(varnums==v);
    % Pull the data for variable v
    vdata = squeeze(dataday(:,column,:));
%     % How much data does each subject have?
%     percentavail(:,v) = sum(~isnan(vdata(start_tt:end_tt,:)),1)/length(daytime)';
%     gooddata = percentavail(:,v)>0;
%     goodindices(:,v) = gooddata;
    vdata = nanmean(vdata(start_tt:end_tt,:),1);
    tempdata(:,v) = vdata';
end

nancheck = any(isnan(tempdata),2);
goodindices = ~nancheck;

category = fungus(goodindices);
data = tempdata(goodindices,:);

if nv==2
    T = table(data(:,1),data(:,2),'VariableNames',varofinterest);
elseif nv==3
    T = table(data(:,1),data(:,2),data(:,3),'VariableNames',varofinterest);
end

% Regression Tree
tree = fitrtree(T,category,'MaxNumSplits',8,'MinLeafSize',50);
view(tree,'Mode','graph')
% Check Regression Tree
[pX,pIdx] = datasample(data,500);
if nv==2
    T2 = table(pX(:,1),pX(:,2),'VariableNames',varofinterest);
elseif nv==3
    T2 = table(pX(:,1),pX(:,2),pX(:,3),'VariableNames',varofinterest);
end
label = predict(tree,T2);
[X,Y,~,AUC] = perfcurve(category(pIdx),label,1);
plot(X,Y)
xlabel('False positive rate'); ylabel('True positive rate');
title(['AUC: ' num2str(round(AUC,2))])
fprintf(['AUC: ' num2str(round(AUC,2)) '\n'])

% Scatter Plot
figure()
if nv==3
    scatter3(data(category==1,1),data(category==1,2),data(category==1,3),36,'r','filled')
    hold on
    scatter3(data(category==0,1),data(category==0,2),data(category==0,3),36,'bo')
    xlabel(varofinterest(1))
    ylabel(varofinterest(2))
    zlabel(varofinterest(3))
elseif nv==2    
    scatter(data(category==1,1),data(category==1,2),36,'r','filled')
    hold on
    scatter(data(category==0,1),data(category==0,2),36,'bo')
    xlabel(varofinterest(1))
    ylabel(varofinterest(2))
end
