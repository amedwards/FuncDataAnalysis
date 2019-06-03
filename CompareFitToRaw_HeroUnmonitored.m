load('UnmonitoredRawDataFilledIn.mat') % This was taken from runFPCA from gooddata around line 502
load('UnmonitoredPCOutputCurves.mat') % This was taken from PlotReducedBabyPCLines in FuncDataAnalysis/Figures/RawBabyBsplines/Version2
load('UnmonitoredAllRowNames.mat')

daymortality = 7;
if daymortality == 7
    load('Unmonitored7DayMortality.mat') % category from runFPCA
    load('UnmonitoredAllPiHat7DayMortality.mat') % allpihat from runFPCA
elseif daymortality == 30
    load('Unmonitored30DayMortality.mat') % category from runFPCA
    load('UnmonitoredAllPiHat30DayMortality.mat') % allpihat form runFPCA
end

A = UnmonitoredPCOutputCurves;
B = UnmonitoredRawDataFilledIn;
cultures = size(A,1);
RMSE = zeros(cultures,1);
for c=1:cultures
    RMSE(c) = sqrt(mean((A(c,:)-B(c,:)).^2));
end

scatter(allpihat(category==2,24),RMSE(category==2),[],'b.')
hold on
scatter(allpihat(category==1,24),RMSE(category==1),[],'r.')
xlabel('Predicted Risk')
ylabel('Poorness of Fit')
title([ num2str(daymortality) ' Day Mortality'])