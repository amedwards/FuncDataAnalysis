load('rng_stats.mat')
% Categories:
% 1. 30 day survival
% 2. Pos/Neg blood culture
% 3. Negsep
% 4. Non-Cons Bacteria
% 5. Gram Neg vs Everything Else
% 6. 7 Day Survival
category = 1;
% Models:
% 1.  'BWT';
% 2.  'BWT + 2 Neighborhoods';
% 3.  'BWT + 3 Neighborhoods';
% 4.  'BWT + 4 Neighborhoods';
% 5.  'BWT + 5 Neighborhoods';
% 6.  'BWT + 6 Neighborhoods';
% 7.  'BWT + 7 Neighborhoods';
% 8.  'BWT + 8 Neighborhoods';
% 9.  'BWT + 9 Neighborhoods';
% 10. 'BWT + 10 Neighborhoods';
% 11. 'BWT + 11 Neighborhoods';
% 12. 'BWT + 12 Neighborhoods';
% 13. 'BWT + Last Value';
% 14. 'BWT + Last Value + 2 Neighborhoods';
% 15. 'BWT + Last Value + 3 Neighborhoods';
% 16. 'BWT + Last Value + 4 Neighborhoods';
% 17. 'BWT + Last Value + 5 Neighborhoods';
% 18. 'BWT + Last Value + 6 Neighborhoods';
% 19. 'BWT + Last Value + 7 Neighborhoods';
% 20. 'BWT + Last Value + 8 Neighborhoods';
% 21. 'BWT + Last Value + 9 Neighborhoods';
% 22. 'BWT + Last Value + 10 Neighborhoods';
% 23. 'BWT + Last Value + 11 Neighborhoods';
% 24. 'BWT + Last Value + 12 Neighborhoods';
% 25. 'BWT/DOA';
% 26. 'BWT/DOA + 2 Neighborhoods';
% 27. 'BWT/DOA + 3 Neighborhoods';
% 28. 'BWT/DOA + 4 Neighborhoods';
% 29. 'BWT/DOA + 5 Neighborhoods';
% 30. 'BWT/DOA + 6 Neighborhoods';
% 31. 'BWT/DOA + 7 Neighborhoods';
% 32. 'BWT/DOA + 8 Neighborhoods';
% 33. 'BWT/DOA + 9 Neighborhoods';
% 34. 'BWT/DOA + 10 Neighborhoods';
% 35. 'BWT/DOA + 11 Neighborhoods';
% 36. 'BWT/DOA + 12 Neighborhoods';
% 37. 'BWT/DOA + Last Value';
% 38. 'BWT/DOA + Last Value + 2 Neighborhoods';
% 39. 'BWT/DOA + Last Value + 3 Neighborhoods';
% 40. 'BWT/DOA + Last Value + 4 Neighborhoods';
% 41. 'BWT/DOA + Last Value + 5 Neighborhoods';
% 42. 'BWT/DOA + Last Value + 6 Neighborhoods';
% 43. 'BWT/DOA + Last Value + 7 Neighborhoods';
% 44. 'BWT/DOA + Last Value + 8 Neighborhoods';
% 45. 'BWT/DOA + Last Value + 9 Neighborhoods';
% 46. 'BWT/DOA + Last Value + 10 Neighborhoods';
% 47. 'BWT/DOA + Last Value + 11 Neighborhoods';
% 48. 'BWT/DOA + Last Value + 12 Neighborhoods';
% 49. 'Last Value';
% 50. 'Last Value + 6 Neighborhoods';
meanAUC = zeros(50,1);
meanAIC = zeros(50,1);
meanBIC = zeros(50,1);
stdAUC = zeros(50,1);
stdAIC = zeros(50,1);
stdBIC = zeros(50,1);
AUCdiffmean = zeros(50,1);
AUCdiffSE = zeros(50,1);
AICdiffmean = zeros(50,1);
AICdiffSE = zeros(50,1);
BICdiffmean = zeros(50,1);
BICdiffSE = zeros(50,1);
for model =1:50
    startvalue = (category-1)*50+model;
    rows = startvalue:300:30000;
    meanAIC(model) = mean(savedallAIC(rows));
    meanBIC(model) = mean(savedallBIC(rows));
    meanAUC(model) = mean(savedallAUC(rows));
    stdAIC(model) = std(savedallAIC(rows));
    stdBIC(model) = std(savedallBIC(rows));
    stdAUC(model) = std(savedallAUC(rows));
    
    baselinemodel = 37;
    startvaluebaseline = (category-1)*50+baselinemodel;
    rowsbaseline = startvaluebaseline:300:30000;
    
    AUCdiffs = savedallAUC(rows)-savedallAUC(rowsbaseline);
    AUCdiffmean(model) = mean(AUCdiffs);
    AUCdiffSE(model) = std(AUCdiffs)/sqrt(length(AUCdiffs));
    
    BICdiffs = savedallBIC(rows)-savedallBIC(rowsbaseline);
    BICdiffmean(model) = mean(BICdiffs);
    BICdiffSE(model) = std(BICdiffs)/sqrt(length(BICdiffs));
    
    AICdiffs = savedallAIC(rows)-savedallAIC(rowsbaseline);
    AICdiffmean(model) = mean(AICdiffs);
    AICdiffSE(model) = std(AICdiffs)/sqrt(length(AICdiffs));
end

