function out = findBW_DOA_priorprob(bwt,doa,category,id)

event = category==1;

allvar = [bwt,doa,event];

bw_cat_max = 400:200:1600; % Birthweight category maxima
% bw_cat_max = 500:200:1500; % Birthweight category maxima
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
        relevantids = id(incat);
        for c=1:sum(incat)
            id_to_leave_out = relevantids(c);
            indices_to_leave_out = (relevantids==id_to_leave_out);
            relvar_with_some_leftout = relevantvar;
            relvar_with_some_leftout(indices_to_leave_out,:) = []; % Leave out the sample we are finding the prior probability for
            indicesincat = find(incat);
            if ~isempty(relvar_with_some_leftout)
                pp_matrix(indicesincat(c),:) = [nanmean(relvar_with_some_leftout,1) sum(incat)-sum(indices_to_leave_out)];
            else
                nanmatrix = ones(1,size(relevantvar,2))*nan;
                pp_matrix(indicesincat(c),:) = [nanmatrix sum(incat)-sum(indices_to_leave_out)];
            end
        end
    end
end

out = pp_matrix(:,3); % prior probability of event given BW and DOA