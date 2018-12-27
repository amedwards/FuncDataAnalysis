daysofage = pdate-bd(pnum);
diedin30days = ddate(pnum)<pdate+30;
diedin7days = ddate(pnum)<pdate+7;

pbw = bwt(pnum);
pega = ega(pnum);

allvar = [pbw,daysofage,pega,clinsep,negsep,sep,diedin30days,diedin7days,cons,pabx];
gooddata = allvar(goodindices,:);

bw_cat_max = 400:100:1500; % Birthweight category maxima
bw_cat_min = bw_cat_max-100; % Birthweigh category minima
num_bw_cats = length(bw_cat_max); % Number of birthweight categories
doa_cat_max = 3:3:222; % Days of age category maxima
doa_cat_min = doa_cat_max-3; % Days of age category minima
num_doa_cats = length(doa_cat_max); % Number of days of age categories

pp_matrix = zeros(sum(goodindices),size(allvar,2)+1); % Prior Probability Matrix

for b=1:num_bw_cats
    for d=1:num_doa_cats
        incat = gooddata(:,1)<bw_cat_max(b) & gooddata(:,1)>=bw_cat_min(b) & gooddata(:,2)<doa_cat_max(d) & gooddata(:,2)>doa_cat_min(d);
        relevantvar = gooddata(incat,:);
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

hoodlabels = resultsLWEA(:,5);
maxh = max(hoodlabels);
par_in_hood = zeros(maxh,size(allvar,2));
pp_in_hood = zeros(maxh,size(allvar,2));
for h=1:maxh
    par_in_hood(h,:) = mean(gooddata(hoodlabels==h,:),1);
    pp_in_hood(h,:) = nanmean(pp_matrix(hoodlabels==h,1:end-1),1);
end

ob = [par_in_hood(:,1:3) par_in_hood(:,4:end)*100]; % Observed
e = [pp_in_hood(:,1:3) pp_in_hood(:,4:end)*100]; % Expected
chisquared = sum(((ob-e).^2)./e);
p = chi2cdf(chisquared,5,'upper');

excel_observed = par_in_hood';
excel_priorprob = pp_in_hood';
excel_p = p';