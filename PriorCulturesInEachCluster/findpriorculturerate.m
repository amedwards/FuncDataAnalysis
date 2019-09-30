monitored = 1;

if monitored
    load('selecteddatamonitored.mat')
else
    load('selecteddataunmon.mat')
end
pnum = pnum(goodindices);
pdate = pdate(goodindices);
if monitored
    clusters = clusters(:,5);
end
precedingcultures = zeros(6,1);
numwprecedingcultures = zeros(6,1);

uniquenums = unique(pnum);
for u=1:length(uniquenums)
    num = uniquenums(u);
    relevantrows = pnum==num;
    clusternums = clusters(relevantrows);
    culturedates = pdate(relevantrows);
    for d=1:length(culturedates)
        timediff = culturedates(d)-culturedates;
        numprecedingcultures = sum(timediff<5&timediff>0);
        precedingcultures(clusternums(d)) = precedingcultures(clusternums(d))+numprecedingcultures;
        if numprecedingcultures
            numwprecedingcultures(clusternums(d)) = numwprecedingcultures(clusternums(d))+1;
        end
    end
end

culturesincluster = [sum(clusters==1);sum(clusters==2);sum(clusters==3);sum(clusters==4);sum(clusters==5);sum(clusters==6)];
percentwithprecedingcultures = numwprecedingcultures./culturesincluster;