numclusters = 1:10;
percvarexplained = [72,9,5,3,2,2,1,1,1,1];
plot(numclusters,percvarexplained,'LineWidth',2,'Color','k')
title('Scree Plot')
xlabel('Principal Component')
ylabel('Percent of Variability Explained')
ylim([0 100])