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

pp_matrix = zeros(size(goodindices,1),size(allvar,2)+1); % Prior Probability Matrix

for b=1:num_bw_cats
    for d=1:num_doa_cats
        incat = pbw<bw_cat_max(b) & pbw>=bw_cat_min(b) & daysofage<doa_cat_max(d) & daysofage>doa_cat_min(d);
        relevantvar = allvar(incat,:);
        for c=1:sum(incat)
            relvar_leftout = relevantvar;
            relvar_leftout(c,:) = []; % Leave out the sample we are finding the prior probability for
            indicesincat = find(incat);
            pp_matrix(indicesincat(c),:) = [nanmean(relvar_leftout,1) sum(incat)];
        end
    end
end
