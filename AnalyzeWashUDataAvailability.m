%  --------------- Load Data and Set Up Workspace  ----------------

addpath ('X:\Amanda\FuncDataAnalysis\Ramsay') % this adds Ramsay's scripts to the path
clear

% Choose which dataset to load
% load('X:\Amanda\FuncDataAnalysis\UVA_FPCA_Stats.mat')
% load('X:\Amanda\FuncDataAnalysis\CU_FPCA_Stats.mat') % time stamps are not corrected for chronologic age yet
load('X:\Amanda\FuncDataAnalysis\WUSTL_FPCA_Stats_plos.mat') % time stamps are not corrected for chronologic age yet

% Load Pete IDs
% load('X:\Amanda\HistogramsAndHeatmaps\peteids.mat')

weeks = 6;
maxday = 7*weeks; 
dayrange = [1,maxday];
daytime = 1:maxday;
nbasis = ceil(maxday/2);

% Keep only the 502 Pete babies
if exist('peteids','var')
    peterows = ismember(pid,peteids); % ex. data for baby 1089 is stored in row 1081 out of 6755
else
    peterows = ones(n,1);
end

percentavail = zeros(n,nv);
daysdataperbaby = zeros(n,nv);
vdata_interp = zeros(maxday,nv,size(dataday,3));
fdstruct = struct();
pcastruct = struct();
thresh = 0.5;
nharm  = 4;

for v=1:nv
    % Pull the data for variable v
    vdata = squeeze(dataday(daytime,v,:));
    
    % Plot the number of babies with data available each day
    figure(); plot(sum(~isnan(vdata),2))
    xlabel('Day of age'); ylabel('Number of babies from which we have data that day'); title(vname{v})
    
    % How much data does each baby have?
    daysdataperbaby(:,v) = sum(~isnan(vdata),1);
    percentavail(:,v) = sum(~isnan(vdata),1)/maxday';
    
    % Interpolate/Extrapolate for missing data
    vdata_interp(:,v,:) = fillmissing(vdata,'linear','SamplePoints',daytime,'EndValues','nearest');
    
    % Find babies which have a certain percentage of data and are Pete babies
    gooddata = percentavail(:,v)>=thresh;
    gooddata = and(gooddata,peterows);
    gooddata = squeeze(vdata_interp(:,v,gooddata));
end

bar([daysdataperbaby(:,1),plos])
ylabel('Days')
legend('Days of HR data','Length of Stay')
xlabel('Baby ID')

figure()
bar([daysdataperbaby(:,2),plos])
ylabel('Days')
legend('Days of SPO2% data','Length of Stay')
xlabel('Baby ID')

figure()
hist(percentavail)
legend('HR','SPO2')
xlabel('Percentage of 6 weeks available')
ylabel('Number of Babies')