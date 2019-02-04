% Create random 5 day samples
load('X:\Amanda\FuncDataAnalysis\Hero\rctrawhero.mat')
seeds = 2000;
hours = 0:24*5-1;
s = rng;
rng(s);
goodlist = ones(seeds,1);
xselected = zeros(seeds,length(hours));
xtselected = zeros(seeds,length(hours));

u1 = randi(length(a),1,seeds);
duplicates = 0;
w = waitbar(0,['Working on file 0 of ' num2str(length(u1))]);
for u=1:length(u1)
    waitbar(u/length(u1),w,['Working on file ' num2str(u) ' of ' num2str(length(u1))]);
    i = u1(u);
	ai = a(i);
    xti = xt(i);
    strike = 0;
    for h=1:length(hours)
        ii = find((a==ai)&(xt==xti+hours(h))); % find the next hour for the same patient
        if isempty(ii)
            strike = strike + 1;
            xselected(u,h) = nan;
            xtselected(u,h) = nan;
        else
            if length(ii)>1
                xselected(u,h) = x(ii(1));
                [u,h]
            else
                xselected(u,h) = x(ii);
                xtselected(u,h) = xt(ii);
            end
        end
    end
    if strike>12 % 13 strikes, you're out
        goodlist(u) = 0;
    end
end
waitbar(1,w,'Complete!')
% goodlist
% xselected
save('randomsegments3','goodlist','xselected','xtselected','s','u1','a')