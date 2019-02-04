function excel_p = makePriorProbabilityTable(goodlist,xtselected,a,u1,bwt,ttd,clusterGMM)
% load('X:\Amanda\FuncDataAnalysis\Hero\randomsegments2.mat')
goodlist = logical(goodlist);
daysofage = floor(xtselected(goodlist,1)/24);

pbw = bwt(a(u1(goodlist)));

pttd = ttd(u1(goodlist)); % Time to death in hours
diedin30days = pttd<=720;
diedin7days = pttd<=24*7;

allvar = [pbw,daysofage,diedin30days,diedin7days];

bw_cat_max = 400:200:1600; % Birthweight category maxima
bw_cat_min = bw_cat_max-200; % Birthweigh category minima
num_bw_cats = length(bw_cat_max); % Number of birthweight categories
doa_cat_max = 4:4:224; % Days of age category maxima
doa_cat_min = doa_cat_max-5; % Days of age category minima
num_doa_cats = length(doa_cat_max); % Number of days of age categories

pp_matrix = zeros(size(allvar,1),size(allvar,2)+1); % Prior Probability Matrix

for b=1:num_bw_cats
    for d=1:num_doa_cats
        incat = allvar(:,1)<bw_cat_max(b) & allvar(:,1)>=bw_cat_min(b) & allvar(:,2)<doa_cat_max(d) & allvar(:,2)>=doa_cat_min(d);
        relevantvar = allvar(incat,:);
        for c=1:sum(incat)
            relvar_leftout = relevantvar;
            relvar_leftout(c,:) = []; % Leave out the sample we are finding the prior probability for
            indicesincat = find(incat);
            if ~isempty(relvar_leftout)
                pp_matrix(indicesincat(c),:) = [nanmean(relvar_leftout,1) sum(incat)];
            else
                nanmatrix = ones(1,length(relevantvar))*nan;
                pp_matrix(indicesincat(c),:) = [nanmatrix sum(incat)];
            end
        end
    end
end

hoodlabels = clusterGMM;
maxh = max(hoodlabels);
par_in_hood = zeros(maxh,size(allvar,2));
pp_in_hood = zeros(maxh,size(allvar,2));
for h=1:maxh
    par_in_hood(h,:) = mean(allvar(hoodlabels==h,:),1);
    pp_in_hood(h,:) = nanmean(pp_matrix(hoodlabels==h,1:end-1),1);
end

ob = [par_in_hood(:,1:2) par_in_hood(:,3:end)*100]; % Observed
e = [pp_in_hood(:,1:2) pp_in_hood(:,3:end)*100]; % Expected

chisquared = nansum(((ob-e).^2)./e);
chisquared(isnan(chisquared)) = 0;
p = chi2cdf(chisquared,5,'upper');

excel_observed = ob';
excel_priorprob = e';
excel_p = p';