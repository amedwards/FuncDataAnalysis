load('prob_of_outcome_baseline') % Birthweight only
base = prob_of_outcome_baseline;
load('prob_of_outcome_lasthero') % Birthweight and last hero score
base = prob_of_outcome_lasthero;
load('prob_of_outcome_fullmodel') % Birthweight and b-spline coefficients
full = prob_of_outcome_fullmodel;
% load('prob_of_outcome_fullmodel') % Birthweight and b-spline coefficients
% full = prob_of_outcome;


scatter(base(diedin30days(goodindices)),full(diedin30days(goodindices)),8,'r','o')
hold on
scatter(base(~diedin30days(goodindices)),full(~diedin30days(goodindices)),8,'b','o')
axis('square')
xlim([0 0.7])
ylim([0 0.7])
xline = 0:0.1:1;
plot(xline,xline)
xlabel('Mortality Probability: Logistic Regression on Birthweight + Last Hero Score')
ylabel('Mortality Probability: Logistic Regression on Birthweight + B-splines')
bigbabies = bwt(pnum(goodindices))>1000;
scatter(base(bigbabies),full(bigbabies),8,'k','.')
daysofage = pdate-bd(pnum);
oldbabies = daysofage(goodindices)>30;
scatter(base(oldbabies),full(oldbabies),8,'m','.')

figure();
groups = 10;
probs = full;

n = length(probs);
n_group = round(n/groups);
diedin30daysdataset = diedin30days(goodindices);
actualprobofdeath = zeros(groups,1);
predictedprobofdeath = zeros(groups,1);
[psorted,I] = sort(probs);
for g=1:groups
    starti = (g-1)*n_group+1;
    endi = (g)*n_group;
    if g==groups
        endi = n;
    end
    group = psorted(starti:endi);
    key = I(starti:endi);
    numingroupdied = sum(diedin30daysdataset(key));
    numingroup = endi-starti;
    actualprobofdeath(g) = numingroupdied/numingroup;
    predictedprobofdeath(g) = mean(group);
end

plot(actualprobofdeath,predictedprobofdeath,'bo-');
hold on
xlabel('Actual Probability of Death in Decile')
ylabel('Predicted Probability of Death in Decile')

probs = base;
n = length(probs);
n_group = round(n/groups);
diedin30daysdataset = diedin30days(goodindices);
actualprobofdeath = zeros(groups,1);
predictedprobofdeath = zeros(groups,1);
[psorted,I] = sort(probs);
for g=1:groups
    starti = (g-1)*n_group+1;
    endi = (g)*n_group;
    if g==groups
        endi = n;
    end
    group = psorted(starti:endi);
    key = I(starti:endi);
    numingroupdied = sum(diedin30daysdataset(key));
    numingroup = endi-starti;
    actualprobofdeath(g) = numingroupdied/numingroup;
    predictedprobofdeath(g) = mean(group);
end
plot(actualprobofdeath,predictedprobofdeath,'ro-');
plot(xline,xline,'k')


legend('Birthweight + B-splines','Birthweight + Last Hero','Identity line')