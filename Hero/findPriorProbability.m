% load('X:\Amanda\FuncDataAnalysis\Hero\randomsegments2.mat')
goodlist = logical(goodlist);
daysofage = floor(xtselected(goodlist,1)/24);

pbw = bwt(a(u1(goodlist)));

pttd = ttd(u1(goodlist)); % Time to death in hours
diedin30days = pttd<=720;

allvar = [pbw,diedin30days];

bw_cat_max = 400:200:1600; % Birthweight category maxima
bw_cat_min = bw_cat_max-200; % Birthweigh category minima
num_bw_cats = length(bw_cat_max); % Number of birthweight categories
doa_cat_max = 5:5:225; % Days of age category maxima
doa_cat_min = doa_cat_max-5; % Days of age category minima
num_doa_cats = length(doa_cat_max); % Number of days of age categories

pp_matrix = zeros(size(allvar,1),size(allvar,2)+1); % Prior Probability Matrix

for b=1:num_bw_cats
    for d=1:num_doa_cats
        incat = pbw<bw_cat_max(b) & pbw>=bw_cat_min(b) & daysofage<doa_cat_max(d) & daysofage>=doa_cat_min(d);
        relevantvar = allvar(incat,:);
        for c=1:sum(incat)
            relvar_leftout = relevantvar;
            relvar_leftout(c,:) = []; % Leave out the sample we are finding the prior probability for
            indicesincat = find(incat);
            pp_matrix(indicesincat(c),:) = [nanmean(relvar_leftout,1) sum(incat)];
        end
    end
end

pp_matrix