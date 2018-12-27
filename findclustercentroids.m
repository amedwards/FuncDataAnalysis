function cgroup = findclustercentroids(clusterlist,harmscr)
clustergroups = size(clusterlist,2);
cgroup = struct;
for cg=1:clustergroups
    clusters = unique(clusterlist(:,cg));
    for c=1:max(clusters)
        a = harmscr(clusterlist(:,cg)==c,:);
        cgroup(cg).centroid(c,:) = mean(a,1);
    end
end
        
