%  --------------- Load Data and Set Up Workspace  ----------------

% Add Ramsay's scripts to the path
addpath ('X:\Amanda\FuncDataAnalysis\Ramsay')
clear
close all

rng(4); % this sets the starting parameters the same way each time for reproducibility

% DISPLAY OR NON-DISPLAY HERO
display = 1;
% Dataset:
% 4 = Display Hero Only, exclude bad pg's, only blood cultures
% 11= Non-Display Hero Only, exclue bad pg's, only blood cultures
% Algorithm: 
% 4 = randomk with LWEA (what is used for the Display group)
% 6 = bagged trees classifier from Display Hero (used for non-display Hero)
if display
    dataset = 4;
    alg = 4; 
else
    dataset = 11; 
    alg = 6; 
end

% Select Grouping:
% 6 - 30 day survival
% 7 - Pos/Neg blood culture
% 8 - negsep
% 9 - non-CONS bacteria
% 12 - gram negative vs everything else
% 14 - 7 day survival

% If we want to compute the outcomes for multiple categories at once, set multicategory to 1
multicategory = 1;
if multicategory
%     grouping = [7,15,12,16];
%     grouping = [6,7,8,9,12,14];
    grouping = [14,6];
else
    grouping = 14;
end

% Choose our data window
daysbefore = -5;
daysafter = 0;

% Percent threshold for the amount of data that is required for inclusion
thresh = 0.50;

% Load HeRO-Specific Dataset and set HeRO-specific parameters
load('X:\Amanda\FuncDataAnalysis\Hero\allcxhero.mat')
load('X:\Amanda\FuncDataAnalysis\Hero\ddate.mat')
load('X:\Amanda\FuncDataAnalysis\Hero\otherfiguredata.mat')
hero(hero==-1) = nan;
vdata = hero';
n = size(vdata,2);
vname = 'HeRO';
limofinterest = [0 7];
if dataset==4
    toinclude = pg>0 & control & c==1; % pg>0 removes cmenu 16-20, control => Hero DISPLAY patients (control is labeled backwards), c==1 indicates blood culture
elseif dataset==11
    load('X:\Amanda\FuncDataAnalysis\Hero\harmfdParHero.mat')
    toinclude = pg>0 & ~control & c==1; % pg>0 removes cmenu 16-20, control = Hero NON-DISPLAY patients (control is labeled backwards), c==1 indicates blood culture
end

% Find the running mean before the data window begins
start_tt_mean = 1;
end_tt_mean = find(tt==daysbefore);

% Keep only the window of data we want
start_tt = find(round(tt,4)==daysbefore);
end_tt = find(round(tt,4)==daysafter);
ttwindow = tt(start_tt:end_tt);
maxday = max(ttwindow);
minday = min(ttwindow);
dayrange = [minday,maxday];
daytime = ttwindow;

% Select Functional Data Parameters
reductionfactor = 6;
nbasis = ceil(length(daytime)/reductionfactor);
nharm = 4;

% ------------------- Create Basis Functions -----------------------

basis = create_bspline_basis(dayrange,nbasis,5); % create a bspline basis - need to up the order to 5 for smoothing penalty requirement if we want derivatives

% Initialize empty arrays and structs
fdstruct = struct();
pcastruct = struct();
lastvalues = zeros(n,1);

% --------------------- Data Preparation ----------------------------

% Plot the number of babies with data available each day
% figure(); plot(sum(~isnan(vdata),2))
% xlabel('Day of age'); ylabel('Number of babies from which we have data that day'); title(vname{v})

% How much data does each baby have?
percentavail = sum(~isnan(vdata(start_tt:end_tt,:)),1)/length(daytime);

% Interpolate/Extrapolate for missing data
vdata_interp_all = fillmissing(vdata,'linear','SamplePoints',tt,'EndValues','nearest');

% Subtract off the mean value from the days preceeding the dataset
priormean = nanmean(vdata_interp_all(start_tt_mean:end_tt_mean,:),1);
vdata_interp = vdata_interp_all(start_tt:end_tt,:);

% Find babies which have a certain percentage of data
gooddata = percentavail>=thresh;
gooddata = gooddata'&toinclude;
goodindices = gooddata;

% Subtract off the data mean
datamean = nanmean(squeeze(vdata_interp));

gooddata = squeeze(vdata_interp(:,goodindices));

% ------------------ Fit Basis Functions ---------------------------

% Fit basis functions for each baby
fdstruct.dayfd = smooth_basis(daytime, gooddata, basis);
dayfd_fdnames = {'Day','ID',vname};
fdstruct.dayfd = putnames(fdstruct.dayfd, dayfd_fdnames);

% ------------------------- PCA ------------------------------------

% If pca_fd is passed 3 arguments, it will use the principal components
% from the display group
if dataset==4
    pcastruct.daypcastr = pca_fd(fdstruct.dayfd,nharm);
elseif dataset == 11
    pcastruct.daypcastr = pca_fd(fdstruct.dayfd,nharm,harmfdPar);
end
figure(); 
[A,B] = plot_pca_fd(pcastruct.daypcastr,1);
meanfd_fdmat = A;
pc_fdmat = B;

% ----------------------- Outcomes --------------------------------
for g=1:size(grouping,2)
    switch grouping(g)
        case 6
            group_names = ['All                 '; 'mortality in 30 days'; 'Survival            '];
            diedinXdays = ddate(pnum)<pdate+30;
            category = double(diedinXdays(goodindices));
            category1 = find(category==1); % Died within 30 days of time 0
            category2 = find(category==0); % Survived for 30 days after time 0
            category(category==0) = 2; % Switch category label so that survival is category 2
        case 7
            group_names = ['All             ';'Positive Culture';'Negative Culture'];
            category = pg(goodindices);
            % Switch the category labels
            if max(category)==2
                category(category==2) = 0;
                category(category==1) = 2;
                category(category==0) = 1;
                category1 = find(category==1); % Positive culture
                category2 = find(category==2); % Negative culture
            else
                category1 = find(category==1); % Positive culture
                category2 = find(category==0); % Negative culture
                category(category==0) = 2;
            end
        case 8
            group_names = ['All       ';'Negsep    ';'Not negsep'];
            category = double(negsep(goodindices));
            category1 = find(category==1); % Negsep
            category2 = find(category==0); % Not negsep (Clinsep + Sep)
            category(category==0) = 2;
        case 9
            group_names = ['All              ';'Non-CONS Bacteria';'Negative or CONS '];
            nonconsbac = pres>2&pres<16;
            category = double(nonconsbac(goodindices));
            category1 = find(category==1); % Non-CONS Bacteria
            category2 = find(category==0); % Either CONS or no bacteria
            category(category==0) = 2;
        case 12
            group_names = ['All            ';'Gram negative  ';'Everything else'];
            gramneg = pres>6 & pres<15;
            category = double(gramneg(goodindices));
            category1 = find(category==1); % Gram Negative
            category2 = find(category==0); % Not Gram Negative (either Negative culture or Gram Positive)
            category(category==0) = 2;
        case 14
            group_names = ['All                '; 'mortality in 7 days'; 'Survival           '];
            if exist('pdate')
                diedinXdays = ddate(pnum)<pdate+7;
            end
            category = double(diedinXdays(goodindices));
            category1 = find(category==1); % Died within 30 days of time 0
            category2 = find(category==0); % Survived for 30 days after time 0
            category(category==0) = 2; % Switch category label so that survival is category 2
    end
    % Set up a design matrix having a column for the grand mean and a 
    % column for each category. Add a dummy contsraint observation.
    p = size(group_names,1);
    zmat = zeros(size(gooddata,2),p);
    zmat(:,1) = 1;
    zmat(category1,2) = 1;
    zmat(category2,3) = 1;

    c_struct(g).group_names = group_names;
    c_struct(g).category = category;
    c_struct(g).category1 = category1;
    c_struct(g).category2 = category2;
    c_struct(g).p = p;
    c_struct(g).zmat = zmat;
end
% Attach a row of 0, 1, 1 to force gender/site effects to sum to zero, 
% & define first regression function as grand mean for both genders/sites
zend = ones(1,p);
zend(1) = 0;
zmat = [zmat; zend];

% Revise YFDOBJ by adding a zero function
coef   = getcoef(fdstruct.dayfd);  
coefend = [coef,zeros(nbasis,1)];  
fdstruct.dayfd = putcoef(fdstruct.dayfd, coefend);

% --------------- Choose what we will cluster with -------------------

harmscr = pcastruct.daypcastr.harmscr;    

% --------------- Predict Probability of Outcome -----------------
switch alg
    case 4
        objects = size(harmscr,1);
        niterk = 100; % number of iterations with different values of k neighborhoods
        mink = 2;
        maxk = round(sqrt(objects));
        randk = randi([mink maxk],niterk,1); % generate 100 random integers between mink and maxk
        allclusters = zeros(objects,niterk);
        for a=1:length(randk)
            [idx] = kmeans(harmscr,randk(a));
            allclusters(:,a) = idx;
        end
        addpath('EnsembleClustering')
        [resultsLWEA,resultsLWGP] = LWEA_and_LWGP(allclusters,category,1);
        C_LWEA = findclustercentroids(resultsLWEA,harmscr);
        C_LWGP = findclustercentroids(resultsLWGP,harmscr); 
        for clusts = 1:11
            figure()
            if exist('pnum')
                percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),squeeze(pc_fdmat),meanfd_fdmat,daytime,vname,category,c_struct,daysbefore,daysafter,limofinterest,fdstruct.dayfd,id(pnum(goodindices)));
            else
                percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),squeeze(pc_fdmat),meanfd_fdmat,daytime,vname,category,c_struct,daysbefore,daysafter,limofinterest,fdstruct.dayfd);
            end
        end
        prob_of_outcome = percent_in_cat(resultsLWEA(:,clusts));
        clusters = resultsLWEA;
    case 6
        load('X:\Amanda\FuncDataAnalysis\Hero\trainedModelHero.mat')
        clusters = trainedModelHero.predictFcn(harmscr);

        C_LWEA = findclustercentroids(clusters,harmscr);
        figure()
        percent_in_cat = NearestNeighborPCs(C_LWEA.centroid,clusters,squeeze(pc_fdmat),meanfd_fdmat,daytime,vname,category,c_struct,daysbefore,daysafter,limofinterest,fdstruct.dayfd,id(pnum(goodindices)));
        prob_of_outcome = percent_in_cat(clusters);
end

% ------------ Create regression model comparison -------------------
if alg ==4 && dataset==4
    daysofage = pdate-bd(pnum);
    windowvdata = squeeze(vdata_interp_all(start_tt:end_tt,:));
    lastvalues = windowvdata(end,:)';

    priorprob = findBW_DOA_priorprob(bwt(pnum(goodindices)),daysofage(goodindices),category,id(pnum(goodindices)));

    n = sum(goodindices);

    regmodel(1).params = bwt(pnum(goodindices));
    regmodel(1).name = 'BWT';

    regmodel(2).params = [bwt(pnum(goodindices)), clusters(:,1)==1];
    regmodel(2).name = 'BWT + 2 Neighborhoods';

    regmodel(3).params = [bwt(pnum(goodindices)), clusters(:,2)==1,clusters(:,2)==2];
    regmodel(3).name = 'BWT + 3 Neighborhoods';

    regmodel(4).params = [bwt(pnum(goodindices)), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
    regmodel(4).name = 'BWT + 4 Neighborhoods';

    regmodel(5).params = [bwt(pnum(goodindices)), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
    regmodel(5).name = 'BWT + 5 Neighborhoods';

    regmodel(6).params = [bwt(pnum(goodindices)), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(6).name = 'BWT + 6 Neighborhoods';

    regmodel(7).params = [bwt(pnum(goodindices)), clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(7).name = 'BWT + 7 Neighborhoods';

    regmodel(8).params = [bwt(pnum(goodindices)), clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(8).name = 'BWT + 8 Neighborhoods';

    regmodel(9).params = [bwt(pnum(goodindices)), clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(9).name = 'BWT + 9 Neighborhoods';

    regmodel(10).params = [bwt(pnum(goodindices)), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(10).name = 'BWT + 10 Neighborhoods';

    regmodel(11).params = [bwt(pnum(goodindices)), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(11).name = 'BWT + 11 Neighborhoods';

    regmodel(12).params = [bwt(pnum(goodindices)), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(12).name = 'BWT + 12 Neighborhoods';

    regmodel(13).params = [bwt(pnum(goodindices)), lastvalues(goodindices)];
    regmodel(13).name = 'BWT + Last Value';

    regmodel(14).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,1)==1];
    regmodel(14).name = 'BWT + Last Value + 2 Neighborhoods';

    regmodel(15).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,2)==1,clusters(:,2)==2];
    regmodel(15).name = 'BWT + Last Value + 3 Neighborhoods';

    regmodel(16).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
    regmodel(16).name = 'BWT + Last Value + 4 Neighborhoods';

    regmodel(17).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
    regmodel(17).name = 'BWT + Last Value + 5 Neighborhoods';

    regmodel(18).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(18).name = 'BWT + Last Value + 6 Neighborhoods';

    regmodel(19).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(19).name = 'BWT + Last Value + 7 Neighborhoods';

    regmodel(20).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(20).name = 'BWT + Last Value + 8 Neighborhoods';

    regmodel(21).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(21).name = 'BWT + Last Value + 9 Neighborhoods';

    regmodel(22).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(22).name = 'BWT + Last Value + 10 Neighborhoods';

    regmodel(23).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(23).name = 'BWT + Last Value + 11 Neighborhoods';

    regmodel(24).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(24).name = 'BWT + Last Value + 12 Neighborhoods';

    regmodel(25).params = priorprob;
    regmodel(25).name = 'BWT/DOA';

    regmodel(26).params = [priorprob, clusters(:,1)==1];
    regmodel(26).name = 'BWT/DOA + 2 Neighborhoods';

    regmodel(27).params = [priorprob, clusters(:,2)==1,clusters(:,2)==2];
    regmodel(27).name = 'BWT/DOA + 3 Neighborhoods';

    regmodel(28).params = [priorprob, clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
    regmodel(28).name = 'BWT/DOA + 4 Neighborhoods';

    regmodel(29).params = [priorprob, clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
    regmodel(29).name = 'BWT/DOA + 5 Neighborhoods';

    regmodel(30).params = [priorprob, clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(30).name = 'BWT/DOA + 6 Neighborhoods';

    regmodel(31).params = [priorprob, clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(31).name = 'BWT/DOA + 7 Neighborhoods';

    regmodel(32).params = [priorprob, clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(32).name = 'BWT/DOA + 8 Neighborhoods';

    regmodel(33).params = [priorprob, clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(33).name = 'BWT/DOA + 9 Neighborhoods';

    regmodel(34).params = [priorprob, clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(34).name = 'BWT/DOA + 10 Neighborhoods';

    regmodel(35).params = [priorprob, clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(35).name = 'BWT/DOA + 11 Neighborhoods';

    regmodel(36).params = [priorprob, clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(36).name = 'BWT/DOA + 12 Neighborhoods';

    regmodel(37).params = [priorprob, lastvalues(goodindices)]; %[bwt(pnum(goodindices)), daysofage(goodindices), lastvalues(goodindices)];
    regmodel(37).name = 'BWT/DOA + Last Value';

    regmodel(38).params = [priorprob, lastvalues(goodindices), clusters(:,1)==1];
    regmodel(38).name = 'BWT/DOA + Last Value + 2 Neighborhoods';

    regmodel(39).params = [priorprob, lastvalues(goodindices), clusters(:,2)==1,clusters(:,2)==2];
    regmodel(39).name = 'BWT/DOA + Last Value + 3 Neighborhoods';

    regmodel(40).params = [priorprob, lastvalues(goodindices), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
    regmodel(40).name = 'BWT/DOA + Last Value + 4 Neighborhoods';

    regmodel(41).params = [priorprob, lastvalues(goodindices), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
    regmodel(41).name = 'BWT/DOA + Last Value + 5 Neighborhoods';

    regmodel(42).params = [priorprob, lastvalues(goodindices), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(42).name = 'BWT/DOA + Last Value + 6 Neighborhoods';

    regmodel(43).params = [priorprob, lastvalues(goodindices), clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(43).name = 'BWT/DOA + Last Value + 7 Neighborhoods';

    regmodel(44).params = [priorprob, lastvalues(goodindices), clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(44).name = 'BWT/DOA + Last Value + 8 Neighborhoods';

    regmodel(45).params = [priorprob, lastvalues(goodindices), clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(45).name = 'BWT/DOA + Last Value + 9 Neighborhoods';

    regmodel(46).params = [priorprob, lastvalues(goodindices), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(46).name = 'BWT/DOA + Last Value + 10 Neighborhoods';

    regmodel(47).params = [priorprob, lastvalues(goodindices), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(47).name = 'BWT/DOA + Last Value + 11 Neighborhoods';

    regmodel(48).params = [priorprob, lastvalues(goodindices), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(48).name = 'BWT/DOA + Last Value + 12 Neighborhoods';

    regmodel(49).params = lastvalues(goodindices);
    regmodel(49).name = 'Last Value';

    regmodel(50).params = [lastvalues(goodindices),clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(50).name = 'Last Value + 6 Neighborhoods';

    alld = zeros(length(regmodel),1);
    alldev = zeros(length(regmodel),1);
    allAIC = zeros(length(regmodel),1);
    allBIC = zeros(length(regmodel),1);
    allBriar = zeros(length(regmodel),1);
    allAUC = zeros(length(regmodel),1);
    allrownames = {};

    for c = 1:length(c_struct)
        for m=1:length(regmodel)
            [B,d,dev,AIC,BIC,Briar,AUC,pihat,X,Y,stats] = modelcomparison(double(regmodel(m).params),c_struct(c).category,n);
            allB(m).B = B;
            alld(m) = d;
            alldev(m) = dev;
            allAIC(m) = AIC;
            allBIC(m) = BIC;
            allBriar(m) = Briar;
            allAUC(m) = AUC;
            allstats(m).stats = stats;
            ROC(m).X = X;
            ROC(m).Y = Y;
            allrownames{m} =  regmodel(m).name;
        end

        alldev = round(alldev,1);
        allAIC = round(allAIC,1);
        allBIC = round(allBIC,1);
        allBriar = round(allBriar,5);
        allAUC = round(allAUC,2);

        disp(c_struct(c).group_names(2,:))
        Table = table(alld,alldev,allAIC,allBIC,allBriar,allAUC,'VariableNames',{'DOF','Dev','AIC','BIC','Briar','AUC'},'rowNames',allrownames);
        disp(Table)
    end

    priorprobmodel = 25;
    basemodel = 37;
    bettermodel = 30;
    bestmodel = 42;

    figure
    plot(ROC(priorprobmodel).X,ROC(priorprobmodel).Y)
    hold on
    plot(ROC(basemodel).X,ROC(basemodel).Y)
    plot(ROC(bettermodel).X,ROC(bettermodel).Y)
    plot(ROC(bestmodel).X,ROC(bestmodel).Y)

    xlabel('False positive rate') 
    ylabel('True positive rate')
    title('ROC for Classification by Logistic Regression')
    legend([regmodel(priorprobmodel).name ': AUC ' num2str(allAUC(priorprobmodel))],...
        [regmodel(basemodel).name  ': AUC ' num2str(allAUC(basemodel))],...
        [regmodel(bettermodel).name  ': AUC ' num2str(allAUC(bettermodel))],...
        [regmodel(bestmodel).name  ': AUC ' num2str(allAUC(bestmodel))])
    axis square
end

if dataset==11
    daysofage = pdate-bd(pnum);
    windowvdata = squeeze(vdata_interp_all(start_tt:end_tt,:));
    lastvalues = windowvdata(end,:)';

    priorprob = findBW_DOA_priorprob(bwt(pnum(goodindices)),daysofage(goodindices),category,id(pnum(goodindices)));

    n = sum(goodindices);

    regmodel(1).params = bwt(pnum(goodindices));
    regmodel(1).name = 'BWT';

    regmodel(2).params = [bwt(pnum(goodindices)), clusters==1,clusters==2,clusters==3,clusters==4,clusters==5];
    regmodel(2).name = 'BWT + 6 Neighborhoods';

    regmodel(3).params = priorprob; %[bwt(pnum(goodindices)), daysofage(goodindices)];
    regmodel(3).name = 'BWT/DOA';

    regmodel(4).params = [priorprob, lastvalues(goodindices)]; %[bwt(pnum(goodindices)), daysofage(goodindices), lastvalues(goodindices)];
    regmodel(4).name = 'BWT/DOA + Last Value';

    regmodel(5).params = [priorprob, clusters==1,clusters==2,clusters==3,clusters==4,clusters==5];
    regmodel(5).name = 'BWT/DOA + 6 Neighborhoods';

    regmodel(6).params = [bwt(pnum(goodindices)), lastvalues(goodindices)];
    regmodel(6).name = 'BWT + Last Value';

    regmodel(7).params = [bwt(pnum(goodindices)), lastvalues(goodindices), clusters==1,clusters==2,clusters==3,clusters==4,clusters==5];
    regmodel(7).name = 'BWT + Last Value + 6 Neighborhoods';

    regmodel(8).params = [priorprob, lastvalues(goodindices), clusters==1,clusters==2,clusters==3,clusters==4,clusters==5];
    regmodel(8).name = 'BWT/DOA + Last Value + 6 Neighborhoods';

    regmodel(9).params = lastvalues(goodindices);
    regmodel(9).name = 'Last Value';

    regmodel(10).params = [lastvalues(goodindices), clusters==1,clusters==2,clusters==3,clusters==4,clusters==5];
    regmodel(10).name = 'Last Value + 6 Neighborhoods';

    alld = zeros(length(regmodel),1);
    alldev = zeros(length(regmodel),1);
    allAIC = zeros(length(regmodel),1);
    allBIC = zeros(length(regmodel),1);
    allBriar = zeros(length(regmodel),1);
    allAUC = zeros(length(regmodel),1);
    allrownames = {};

    for c = 1:length(c_struct)
        for m=1:length(regmodel)
            [B,d,dev,AIC,BIC,Briar,AUC,pihat,X,Y,stats] = modelcomparison(double(regmodel(m).params),c_struct(c).category,n);
            allB(m).B = B;
            alld(m) = d;
            alldev(m) = dev;
            allAIC(m) = AIC;
            allBIC(m) = BIC;
            allBriar(m) = Briar;
            allAUC(m) = AUC;
            allstats(m).stats = stats;
            ROC(m).X = X;
            ROC(m).Y = Y;
            allrownames{m} =  regmodel(m).name;
        end

        alldev = round(alldev,1);
        allAIC = round(allAIC,1);
        allBIC = round(allBIC,1);
        allBriar = round(allBriar,5);
        allAUC = round(allAUC,2);

        disp(c_struct(c).group_names(2,:))
        Table = table(alld,alldev,allAIC,allBIC,allBriar,allAUC,'VariableNames',{'DOF','Dev','AIC','BIC','Briar','AUC'},'rowNames',allrownames);
        disp(Table)
    end

    priorprobmodel = 3;
    basemodel = 4;
    bettermodel = 5;
    bestmodel = 8;

    figure
    plot(ROC(priorprobmodel).X,ROC(priorprobmodel).Y)
    hold on
    plot(ROC(basemodel).X,ROC(basemodel).Y)
    plot(ROC(bettermodel).X,ROC(bettermodel).Y)
    plot(ROC(bestmodel).X,ROC(bestmodel).Y)

    xlabel('False positive rate') 
    ylabel('True positive rate')
    title('ROC for Classification by Logistic Regression')
    legend([regmodel(priorprobmodel).name ': AUC ' num2str(allAUC(priorprobmodel))],...
        [regmodel(basemodel).name  ': AUC ' num2str(allAUC(basemodel))],...
        [regmodel(bettermodel).name  ': AUC ' num2str(allAUC(bettermodel))],...
        [regmodel(bestmodel).name  ': AUC ' num2str(allAUC(bestmodel))])
    axis square
end


function [B,d,dev,AIC,BIC,Briar,AUC,pihat,X,Y,stats] = modelcomparison(modelparams,category,n)
   % Logistic Regression
    category(category==2) = 0;
    category = logical(category);
    [B,dev,stats] = glmfit(modelparams,category,'binomial','link','logit');
    pihat = glmval(B,modelparams,'logit');
    RSS = nansum(stats.resid.^2);
    d = size(modelparams,2);

    [X,Y,~,AUC] = perfcurve(category,pihat,1);

    % BIC from: http://www.stat.wisc.edu/courses/st333-larget/aic.pdf
    % BIC = n+n*log(2*pi)+n*log(RSS/n)+log(n)*(d+1);
    Briar = RSS/n;
    AIC = dev+2*d;
    BIC = dev+log(n)*d;
end