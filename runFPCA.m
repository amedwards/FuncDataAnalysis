%  --------------- Load Data and Set Up Workspace  ----------------

% Add Ramsay's scripts to the path
addpath ('X:\Amanda\FuncDataAnalysis\Ramsay')
clear
close all

rng(4); % this sets the starting parameters the same way each time for reproducibility

% Choose which dataset to load
% 1 = Sepsis/NEC
% 2 = HERO Sepsis Cases
% 3 = All HERO cultures
% 4 = Control Hero Only, exclude bad pg's, only blood cultures
% 5 = UVA, CU, WUSTL Stats - daily
% 6 = UVA, CU, WUSTL Stats - hourly
% 7 = Random subset of 5 day HERO trial samples
% 8 = Sepsis/NEC SPO2, XC SPO2, SD HR
% 9 = Sepsis/NEC SPO2, XC SPO2, SD HR - Off Vent ONLY
dataset = 4; 

% 1-Mean model, 2-Slope model, 3-Bspline model, 4-Last Hero Value, 5-Bspline limited model
model = 3;

% Choose whether you want to subtract off the PRECEEDING mean
subtractoffmean = 0;

% Choose whether to use 1: Principal Components or 0: Raw output from basis functions for clustering
usePCs = 1;

% Subtract (NOT preceeding) mean from data, then add to array to help with prediction
colofmeans = 0;

% Add in birthweight
add_bwt = 0;

% Add in days of age
add_daysofage = 0;

% Add in birthweight/days of age prior probability
add_bwtdoa = 0;

% Add in last hero score
add_lasthero = 0;

% Only use birthweight
only_bwt = 0;

% Only use birthweight + days of age
only_bwt_doa = 0;

% Algorithm: 
% 1 = logistic regression
% 2 = kmeans clustering
% 3 = Gaussian Mixture Model
% 4 = randomk with LWEA
alg = 4;

% Select Grouping:
% 1 - site
% 2 - gender
% 3 - ELBW/VLBW
% 4 - EGA
% 5 - control vs display (for Hero)
% 6 - 30 day survival
% 7 - Pos/Neg blood culture
% 8 - negsep
% 9 - non-CONS bacteria
% 10 - sepsis no vent, nec no vent, sepsis vent, nec vent
% 11 - sepsis (no vent + vent), nec (no vent + vent)
% 12 - gram negative vs everything else
% 13 - Candida/Yeast/Fungus vs everything else
% 14 - 7 day survival
grouping = 14; 

% Choose number of neighborhoods
neighborhoods = 6; % usually 6

% Choose our data window
daysbefore = -5;
daysafter = 0;

% Percent threshold for the amount of data that is required for inclusion
thresh = 0.50;

% Choose whether to weight risk score by distance
weightbydistance = 0;

% Turn on smoothing
smoothingon = 0;

% Use Varimax to rotate the principal components
usevarimax = 0;

% Use the 12 datapoint derivative (12 hour derivative for Hero)
derivative12hr = 0;

derivative = 0;



% Load Dataset
if dataset==1
    load('X:\Amanda\FuncDataAnalysis\SepsisNEC\nearevent.mat')
    
    % Convert Variables to Be Consistent with Previous Code
    vname = xname;
    inst = esite; % 1 = UVA, 2 = CU
    tt = xt'/24; % Convert hours to days
    dataday = permute(x,[2,3,1]); % Arrange the data so that it goes time,variable,subject
    n = size(dataday,3);
    
    % Only Run Analysis on Variables of Interest
    varofinterest = {'Mean HR','Mean SPO2-%'};
    limofinterest = [100 180; 80 100];
    [varlog,varnums] = ismember(vname,varofinterest);
    nv = sum(varlog);
    
    % Set basis parameters
    reductionfactor = 24;
    lambdabase = 10;
elseif dataset==2
	load('X:\Amanda\FuncDataAnalysis\Hero\herosepsis.mat')
    
    % Convert Variables to Be Consistent with Previous Code
    hero(hero==-1) = nan;
    dataday(:,1,:) = hero';
    n = size(dataday,3);
    
    % Only Run Analysis on Variables of Interest
    varofinterest = {'hero'};
    vname = varofinterest;
    limofinterest = [-1 7];
    varnums = 1;
    nv = 1;
    
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
elseif dataset==3
    load('X:\Amanda\FuncDataAnalysis\Hero\allcxhero.mat')
    load('X:\Amanda\FuncDataAnalysis\Hero\ddate.mat')
    
    % Convert Variables to Be Consistent with Previous Code
    hero(hero==-1) = nan;
    dataday(:,1,:) = hero';
    n = size(dataday,3);
    
    % Only Run Analysis on Variables of Interest
    varofinterest = {'hero'};
    vname = varofinterest;
    limofinterest = [-1 7];
    varnums = 1;
    nv = 1;
    
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
elseif dataset==4
    load('X:\Amanda\FuncDataAnalysis\Hero\allcxhero.mat')
    load('X:\Amanda\FuncDataAnalysis\Hero\ddate.mat')
    load('X:\Amanda\FuncDataAnalysis\Hero\otherfiguredata.mat')
    
    % Convert Variables to Be Consistent with Previous Code
    hero(hero==-1) = nan;
    dataday(:,1,:) = hero';
    n = size(dataday,3);
    toinclude = pg>0 & control & c==1; % pg>0 removes cmenu 16-20, control = Hero non-monitored patients, c==1 indicates blood culture
    % Only Run Analysis on Variables of Interest
    varofinterest = {'hero'};
    vname = varofinterest;
    limofinterest = [0 7];
    varnums = 1;
    nv = 1;
    
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
elseif dataset==5
    load('X:\Amanda\FuncDataAnalysis\UVA_CU_WUSTL_FPCA_Stats_day.mat')
    load('X:\Amanda\FuncDataAnalysis\GoodIDs.mat')

    % Keep only the data from qualifying babies who have at least 7 days of
    % data in the first 28 days of life
    dataday = datamerge(:,:,GoodIDs);
    pbd = pbd(GoodIDs);
    pbw = pbw(GoodIDs);
    pega = pega(GoodIDs);
    pgen = pgen(GoodIDs);
    pid = pid(GoodIDs);
    numsamps = numsamps(:,:,GoodIDs);
    inst = inst(GoodIDs);
    n = sum(GoodIDs);
    
    varofinterest = {'Mean HR','Mean SPO2-%'};
    [varlog,varnums] = ismember(vname,varofinterest);
    nv = sum(varlog);
    tt=1:84;
    limofinterest = [100 200; 80 100];
    daysbefore = 1;
    daysafter = 42;
    
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
elseif dataset==6
    load('X:\Amanda\FuncDataAnalysis\UVA_CU_WUSTL_FPCA_Stats_hr.mat')
    load('X:\Amanda\FuncDataAnalysis\GoodIDs.mat')

    % Keep only the data from qualifying babies who have at least 7 days of
    % data in the first 28 days of life
    dataday = datamerge(:,:,GoodIDs);
    pbd = pbd(GoodIDs);
    pbw = pbw(GoodIDs);
    pega = pega(GoodIDs);
    pgen = pgen(GoodIDs);
    pid = pid(GoodIDs);
    numsamps = numsamps(:,:,GoodIDs);
    inst = inst(GoodIDs);
    n = sum(GoodIDs);
    
    varofinterest = {'Mean HR','Mean SPO2-%'};
    [varlog,varnums] = ismember(vname,varofinterest);
    nv = sum(varlog);
    limofinterest = [100 200; 80 100];
    daysbefore = 1;
    daysafter = 24*2; % hours after start of life
    tt = daysbefore:daysafter; % 1st 6 weeks of life by hour
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
elseif dataset==7
    load('X:\Amanda\FuncDataAnalysis\Hero\rctrawhero.mat')
    load('X:\Amanda\FuncDataAnalysis\Hero\randomsegments3.mat')
    goodlist = logical(goodlist);
    dataday(:,1,:) = xselected(goodlist,:)';
    pbd = bd(a(u1(goodlist)));
    pbw = bwt(a(u1(goodlist)));
    pid = id(a(u1(goodlist)));
    gg = g(a(u1(goodlist)));
    
    pttd = ttd(u1(goodlist)); % Time to death in hours
    diedinXdays = pttd<=720;
    n = sum(goodlist);
    varofinterest = {'hero logit'};
    vname = varofinterest;
    varnums = 1;
    nv = 1;
    limofinterest = [-8 3];
    daysbefore = 1; % start of the dataset
    daysafter = 24*5; % hours after the start of the dataset
    tt = daysbefore:daysafter;
    
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
elseif dataset==8
    load('X:\Amanda\FuncDataAnalysis\SepsisNEC\nearevent.mat')
    
    % Convert Variables to Be Consistent with Previous Code
    vname = xname;
    inst = esite; % 1 = UVA, 2 = CU
    tt = xt'/24; % Convert hours to days
    dataday = permute(x,[2,3,1]); % Arrange the data so that it goes time,variable,subject
    n = size(dataday,3);
    
    % Only Run Analysis on Variables of Interest
    varofinterest = {'Max XC HR SPO2-%','Mean SPO2-%','SD HR'};
    limofinterest = [0 0.5; 0 12; 0 100];
    [varlog,varnums] = ismember(vname,varofinterest);
    nv = sum(varlog);
    
    % Set basis parameters
    reductionfactor = 24;
    lambdabase = 10;
elseif dataset==9
    load('X:\Amanda\FuncDataAnalysis\SepsisNEC\nearevent.mat')
    
    % Convert Variables to Be Consistent with Previous Code
    vname = xname;
    inst = esite; % 1 = UVA, 2 = CU
    tt = xt'/24; % Convert hours to days
    dataday = permute(x,[2,3,1]); % Arrange the data so that it goes time,variable,subject
    
    % Only keep babies who are off the ventilator
    dataday = dataday(:,:,etype<3);
    etype = etype(etype<3);
    
    n = size(dataday,3);
    
    % Only Run Analysis on Variables of Interest
    varofinterest = {'Max XC HR SPO2-%','SD HR','Mean SPO2-%'};
    limofinterest = [0 0.5; 0 12; 88 100];
    [varlog,varnums] = ismember(vname,varofinterest);
    nv = sum(varlog);
    
    % Set basis parameters
    reductionfactor = 24;
    lambdabase = 10;
end

if derivative12hr
    dataday = dataday(13:daysafter,:,:)-dataday(1:daysafter-12,:,:);
    tt = tt(13:end);
    if dataset<5
        daysbefore = -4;
    else
        daysbefore = 24; % Hours
        limofinterest = [-12 12; -12 12];
    end
end

if derivative
    dataday = dataday(daysbefore+1:daysafter,:,:)-dataday(1:daysafter-1,:,:);
    tt = tt(2:end);
    if dataset<5
        daysbefore = -4;
    else
        daysbefore = 2; % Hours
        limofinterest = [-12 12; -12 12];
    end
end

% Find the running mean before the data window begins
start_tt_mean = 1;
end_tt_mean = find(tt==daysbefore);

% Keep only the window of data we want
start_tt = find(tt==daysbefore);
end_tt = find(tt==daysafter);
ttwindow = tt(start_tt:end_tt);
maxday = max(ttwindow);
minday = min(ttwindow);
dayrange = [minday,maxday];
daytime = ttwindow;

% Select Functional Data Parameters
switch model
    case 1
        nbasis = 1;
        nharm = 1;
    case 2
        nbasis = 2;
        nharm = 2;
    case 3
        nbasis = ceil(length(daytime)/reductionfactor);
        nharm = 4;
    case 4
        nbasis = 1;
        nharm = 1;
    case 5
        nbasis = 5;
        nharm = 4;
end

% ------------------- Create Basis Functions -----------------------

Lbasis  = create_constant_basis(dayrange);  %  create a constant basis
Mbasis = create_monomial_basis(dayrange,nbasis); % create a monomial basis
if model==3 || model==5
    Bbasis = create_bspline_basis(dayrange,nbasis,5); % create a bspline basis - need to up the order to 5 for smoothing penalty requirement if we want derivatives
end

% Choose which basis function to use
if model==1 || model==2 || model==4
    basis = Mbasis;
else
    basis = Bbasis;
end

%  ----------  Set up the harmonic acceleration operator  --------------

Lcoef   = [0,(4/(length(daytime)))^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object

% Initialize empty arrays and structs
percentavail = zeros(n,nv);
percentavail28 = zeros(n,nv);
vdata_interp = zeros(length(ttwindow),nv,size(dataday,3));
vdata_interp_all = zeros(length(tt),nv,size(dataday,3));
fdstruct = struct();
pcastruct = struct();
datamean = zeros(nv,size(dataday,3));
lastvalues = zeros(n,nv);


for v=1:nv
    column = find(varnums==v);
    % --------------------- Data Preparation ----------------------------
    % Pull the data for variable v
    % vdata = squeeze(dataday(ismember(tt,daytime),column,:));
    vdata = squeeze(dataday(:,column,:));
    
    % Plot the number of babies with data available each day
    % figure(); plot(sum(~isnan(vdata),2))
    % xlabel('Day of age'); ylabel('Number of babies from which we have data that day'); title(vname{v})
    
    % How much data does each baby have?
    percentavail(:,v) = sum(~isnan(vdata(start_tt:end_tt,:)),1)/length(daytime)';
       
    % Interpolate/Extrapolate for missing data
    vdata_interp_all(:,v,:) = fillmissing(vdata,'linear','SamplePoints',tt,'EndValues','nearest');
    
    % If we want to use the model where we just take the last reading
    if model == 4
        windowvdata = squeeze(vdata_interp_all(start_tt:end_tt,v,:));
        vdata_interp_all(:,v,:) = ones(length(tt),1)*windowvdata(end,:);
    end
    
    % Subtract off the mean value from the days preceeding the dataset
    priormean = nanmean(squeeze(vdata_interp_all(start_tt_mean:end_tt_mean,v,:)),1);
    
    if subtractoffmean
        vdata_interp(:,v,:) = squeeze(vdata_interp_all(start_tt:end_tt,v,:))-priormean;
    else
        vdata_interp(:,v,:) = squeeze(vdata_interp_all(start_tt:end_tt,v,:));
    end
    
    % Find babies which have a certain percentage of data
    gooddata = percentavail(:,v)>=thresh;
    if dataset==4
        gooddata = gooddata&toinclude;
    end
    goodindices(v,:) = gooddata;
    
    % Subtract off the data mean
    datamean(v,:) = nanmean(squeeze(vdata_interp(:,v,:)));
    if colofmeans
        vdata_interp(:,v,:) = squeeze(vdata_interp(:,v,:))-datamean(v,:);
    end
    
    if dataset == 1
        fprintf(['UVA   Babies Included: ' num2str(sum(goodindices(v,inst==1))) '\n'])
        fprintf(['CU    Babies Included: ' num2str(sum(goodindices(v,inst==2))) '\n'])
    elseif dataset == 2
        fprintf(['Babies Included: ' num2str(sum(goodindices(v,:))) '\n'])
    end
    
    gooddata = squeeze(vdata_interp(:,v,goodindices(v,:)));
    
    % ------------------ Fit Basis Functions ---------------------------
    
    % Fit basis functions for each baby
    fdstruct(v).dayfd = smooth_basis(daytime, gooddata, basis);
    dayfd_fdnames = {'Day','ID',varofinterest{v}};
    fdstruct(v).dayfd = putnames(fdstruct(v).dayfd, dayfd_fdnames);
    
    % Plot individual basis curves and values
    % figure(); plotfit_fd(gooddata, daytime, fdstruct(v).dayfd)
    
    % Plot all basis curves at once
    figure(); plot(fdstruct(v).dayfd)
    
    % ----------------- Determine Level of Smoothing ------------------
    % Choose level of smoothing using the generalized cross-validation 
    %          criterion with smoothing function smooth_basis.

    if smoothingon
        % Set up range of smoothing parameters in log_10 units
        loglam = (-5:1)';
        nlam   = length(loglam);
        dfsave  = zeros(nlam,1);
        gcvsave = zeros(nlam,1);

        % Loop through smoothing parameters
        for ilam=1:length(loglam)
            lambda = lambdabase^loglam(ilam);
            display(['lambda = ',num2str(lambda)])
            fdParobj = fdPar(basis, harmaccelLfd, lambda);
            [fdobj, df, gcv] = smooth_basis(daytime, gooddata, fdParobj);
            dfsave(ilam)  = df;
            gcvsave(ilam) = sum(gcv);
        end

        % Display and plot degrees of freedom and GCV criterion
        disp('Log lambda    df          gcv')
        disp([loglam, dfsave, gcvsave])

        figure()
        subplot(2,1,1)
        plot(loglam, gcvsave, 'o-')
        ylabel('\fontsize{16} GCV Criterion')
        title(['\fontsize{16} ' vname{column} ' Smoothing'])
        subplot(2,1,2)
        plot(loglam, dfsave, 'o-')
        xlabel('\fontsize{16} log_{10} \lambda')
        ylabel('\fontsize{16} Degrees of Freedom')
    
    % ---------------------- Smooth our data -------------------------

        %  Do final smooth with minimum GCV value
        lambda = 1;  %  minimum GCV estimate
        fdParobj = fdPar(basis, harmaccelLfd, lambda);
        fdstruct(v).dayfd = smooth_basis(daytime, gooddata, fdParobj);

        % Plot data and fit
        % subplot(1,1,1)
        % plotfit_fd(gooddata, daytime, fdstruct(v).dayfd, vname{v})
    
    end
           
    % ------------------------- PCA ------------------------------------
    
    pcastruct(v).daypcastr = pca_fd(fdstruct(v).dayfd,nharm);
    % Rotate the PC's using the varimax function to explain the most variability
    if usevarimax
        pcastruct(v).daypcastr = varmx_pca_fd(pcastruct(v).daypcastr);
    end
    figure(); 
    [meanfd_fdmat,pc_fdmat] = plot_pca_fd(pcastruct(v).daypcastr,1);
    
%     % -------------------------- CCA ------------------------------------
%     if v==2 && nv>=2
%         bothindices = sum(goodindices,1)==2;
%         fd1indices = bothindices(logical(goodindices(1,:)));
%         fd2indices = bothindices(logical(goodindices(2,:)));
%         ccastruct.daypcastr = cca_fd(fdstruct(1).dayfd(fd1indices),fdstruct(2).dayfd(fd2indices),3);
%         figure()
%         plot(ccastruct.daypcastr.wtfdx)
%         hold on
%         plot(ccastruct.daypcastr.wtfdy)
%         
%         group_names = ['All   ';'Sepsis';'NEC   '];
%         category = etype(bothindices);
%         color = (category==1 | category==3); % Sepsis
%         figure()
%         for cwi = 1:3 % canonical weight index
%             subplot(1,3,cwi)
%             scatter(ccastruct.daypcastr.varx(color,cwi),ccastruct.daypcastr.vary(color,cwi)) % Sepsis
%             hold on
%             scatter(ccastruct.daypcastr.varx(~color,cwi),ccastruct.daypcastr.vary(~color,cwi)) % NEC
%             xlabel([varofinterest{1} ' Canonical Weight'])
%             ylabel([varofinterest{2} ' Canonical Weight'])
%             legend(group_names(2:3,:))
%             title(['Canonical Weight Pair ' num2str(cwi) ' with Correlation: ' num2str(ccastruct.daypcastr.corrs(cwi))])
%         end
% %         plot_cca(ccastruct.daypcastr,1);
%     end
    
    
    % -------------------- Functional ANOVA -------------------------
    % Names for groups
    switch grouping
        case 1
            if length(unique(inst))==3
                group_names = ['All  '; 'UVA  ';'CU   '; 'WUSTL'];
                % Site Indices
                category = inst(goodindices(v,:));
                category1 = find(category==1)'; % UVA
                category2 = find(category==2)'; % CU
                category3 = find(category==3)'; % WUSTL
                % Set up a design matrix having a column for the grand mean and a 
                % column for each gender/site. Add a dummy constraint observation.
                p = size(group_names,1);
                zmat = zeros(size(gooddata,2),p);
                zmat(:,1) = 1;
                zmat(category1,2) = 1;
                zmat(category2,3) = 1;
                zmat(category3,4) = 1;
            elseif length(unique(inst))==2
                group_names = ['All'; 'UVA';'CU '];
                % Site Indices
                category = inst(goodindices(v,:));
                category1 = find(category==1)'; % UVA
                category2 = find(category==2)'; % CU
                % Set up a design matrix having a column for the grand mean and a 
                % column for each gender/site. Add a dummy constraint observation.
                p = size(group_names,1);
                zmat = zeros(size(gooddata,2),p);
                zmat(:,1) = 1;
                zmat(category1,2) = 1;
                zmat(category2,3) = 1;
            else
                disp('Number of institutions is not 2 or 3. Edit code to reflect correct number of institutions.')
            end
        case 2
            group_names = ['All   ';'Female'; 'Male  '];
            % Gender Indices
            category = pgen(goodindices(v,:));
            category1 = find(category==1)'; % Female
            category2 = find(category==2)'; % Male
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy constraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 3
            group_names = ['All ';'ELBW'; 'VLBW'];
            % Gender Indices
            category = pbw(goodindices(v,:));
            category1 = find(category<1000)'; % ELBW
            category2 = find(category>=1000)'; % VLBW non ELBW
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site/weight. Add a dummy constraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 4
            group_names = ['All        '; '<27 Weeks  ';'27-30 Weeks';'31-34 Weeks';'>34 Weeks  '];
            % Site Indices
            category = pega(goodindices(v,:));
            category1 = find(category<27)'; % <27 Weeks
            category2 = find(category>=27&category<30)'; % 27-30 Weeks
            category3 = find(category>=30&category<34)'; % 31-34 Weeks
            category4 = find(category>=34)'; % >34 Weeks
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy constraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
            zmat(category3,4) = 1;
            zmat(category4,5) = 1;
        case 5
            group_names = ['All    ';'Control'; 'Display'];
            % Control/Display Index
            category = gg(goodindices(v,:));
            category1 = find(category==1)'; 
            category2 = find(category==2)';
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy constraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 6
            group_names = ['All         '; 'Death in 30d'; 'Survival    '];
            if exist('pdate')
                diedinXdays = ddate(pnum)<pdate+30;
            end
            category = double(diedinXdays(goodindices(v,:)));
            category1 = find(category==1); % Died within 30 days of time 0
            category2 = find(category==0); % Survived for 30 days after time 0
            category(category==0) = 2; % Switch category label so that survival is category 2
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy contsraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 7
            group_names = ['All             ';'Positive Culture';'Negative Culture'];
            category = pg(goodindices(v,:));
            % Switch the category labels
            category(category==2) = 0;
            category(category==1) = 2;
            category(category==0) = 1;
            category1 = find(category==1); % Positive culture
            category2 = find(category==2); % Negative culture
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 8
            group_names = ['All       ';'Negsep    ';'Not negsep'];
            category = double(negsep(goodindices(v,:)));
            category1 = find(category==1); % Negsep
            category2 = find(category==0); % Not negsep (Clinsep + Sep)
            category(category==0) = 2;
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 9
            group_names = ['All              ';'Non-CONS Bacteria';'Negative or CONS '];
            nonconsbac = pres>2&pres<16;
            category = double(nonconsbac(goodindices(v,:)));
            category1 = find(category==1); % Non-CONS Bacteria
            category2 = find(category==0); % Either CONS or no bacteria
            category(category==0) = 2;
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 10
            group_names = ['All           ';'Sepsis No Vent';'NEC No Vent   ';'Sepsis Vent   ';'NEC Vent      '];
            category = etype(goodindices(v,:));
            category1 = find(category==1);
            category2 = find(category==2);
            category3 = find(category==3);
            category4 = find(category==4);
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
            zmat(category3,4) = 1;
            zmat(category4,5) = 1;
        case 11
            group_names = ['All   ';'Sepsis';'NEC   '];
            category = etype(goodindices(v,:));
            category1 = find(category==1 | category==3);
            category2 = find(category==2 | category==4);
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 12
            group_names = ['All            ';'Gram negative  ';'Everything else'];
            gramneg = pres>6 & pres<15;
            category = double(gramneg(goodindices(v,:)));
            category1 = find(category==1); % Candida, yeast, fungus
            category2 = find(category==0); % Not Candida, yeast, fungus
            category(category==0) = 2;
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 13
            group_names = ['All                 ';'Candida Yeast Fungus';'Everything else     '];
            candida = pres==15;
            category = double(candida(goodindices(v,:)));
            category1 = find(category==1); % Candida, yeast, fungus
            category2 = find(category==0); % Not Candida, yeast, fungus
            category(category==0) = 2;
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 14
            group_names = ['All        '; 'Death in 7d'; 'Survival   '];
            if exist('pdate')
                diedinXdays = ddate(pnum)<pdate+7;
            end
            category = double(diedinXdays(goodindices(v,:)));
            category1 = find(category==1); % Died within 30 days of time 0
            category2 = find(category==0); % Survived for 30 days after time 0
            category(category==0) = 2; % Switch category label so that survival is category 2
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy contsraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
    end

    % Attach a row of 0, 1, 1 to force gender/site effects to sum to zero, 
    % & define first regression function as grand mean for both genders/sites
    zend = ones(1,p);
    zend(1) = 0;
    zmat = [zmat; zend];

    % Revise YFDOBJ by adding a zero function
    coef   = getcoef(fdstruct(v).dayfd);  
    coefend = [coef,zeros(nbasis,1)];  
    fdstruct(v).dayfd = putcoef(fdstruct(v).dayfd, coefend);  

    xfdcell = cell(1,p);
    for j=1:p
        xfdcell{j} = zmat(:,j);
    end
    
    if model ==3
        % Set up the basis for the regression functions
        nbetabasis = nbasis; %10;
        betabasis  = create_bspline_basis(dayrange, nbetabasis);

        % Set up the functional parameter object for the regression functions
        betafd    = fd(zeros(nbetabasis,p), betabasis);
        estimate  = 1;
        lambda    = 0;
        betafdPar = fdPar(betafd, harmaccelLfd, lambda, estimate);
        betacell = cell(1,p);
        for j=1:p
            betacell{j} = betafdPar;
        end

        % Compute regression coefficient functions and predicted functions
        fRegressStruct = fRegress(fdstruct(v).dayfd, xfdcell, betacell);
        betaestcell = fRegressStruct.betahat;
        yhatfdobj = fRegressStruct.yhat;

        % Plot regression functions
        for j=1:p
            figure()
            plot(getfd(betaestcell{j}))
            title(['\fontsize{16} ',group_names(j,:)])
        end

        % Plot predicted functions
        figure()
        plot(yhatfdobj)
        title('Predicted values')
    end
    
    % --------------- Choose what we will cluster with -------------------
    
    if usePCs % Principal Components
        harmscr = pcastruct(v).daypcastr.harmscr;    
    else % Raw Coefficients
        harmscr = getcoef(fdstruct(v).dayfd)';
        harmscr = harmscr(1:end-1,:);
    end
    
    if colofmeans
        harmscr = [harmscr, datamean(goodindices(v,:))'];
    end
    
    if add_bwt
        harmscr = [harmscr, bwt(pnum(goodindices(v,:)))];
    end
    
    if add_daysofage
        daysofage = pdate-bd(pnum);
        harmscr = [harmscr, daysofage(goodindices(v,:))];
    end
    
    if add_bwtdoa
        daysofage = pdate-bd(pnum);
        priorprob = findBW_DOA_priorprob(bwt(pnum(goodindices(v,:))),daysofage(goodindices(v,:)),category);
        harmscr = [harmscr, priorprob];
    end
    
    if add_lasthero
        windowvdata = squeeze(vdata_interp_all(start_tt:end_tt,v,:));
        lastscore = windowvdata(end,:);
        harmscr = [harmscr, lastscore(goodindices(v,:))'];
    end
    
    if only_bwt
        harmscr = bwt(pnum(goodindices(v,:)));
    end
    
    if only_bwt_doa
        daysofage = pdate-bd(pnum);
        harmscr = [bwt(pnum(goodindices(v,:))), daysofage(goodindices(v,:))];
    end
    
    % --------------- Predict Probability of Outcome -----------------
    switch alg
        case 1 % Use Logistic Regression
            [B,dev,stats] = mnrfit(harmscr,category);
            pihat = mnrval(B,harmscr,stats);
            prob_of_outcome = pihat(:,1);
            
        case 2 % Use kmeans
            % Centroids Hero 6 neighborhoods, uniform start, 500 iterations
            startcentroids = [-5.71565437617825,-0.621723197077041,-0.155262049420456,0.145934521267938;6.87448810789400,0.393735699329390,-0.202435121306784,0.0447056045919093;2.60630621130926,-1.77446345925976,0.396459429151139,-0.000306584044559931;-2.35153663888559,0.283072268370912,-0.125753064280199,-0.0613435691042324;2.51369991809380,0.899730355949572,-0.210786120701946,0.0522966919589589;-0.267403459070076,0.0598615299867668,0.0625344896221859,-0.00775541810239860];
            [idx,C, ~, D] = kmeans(harmscr,neighborhoods,'Replicates',100,'Start','uniform');
            if usePCs
                percent_in_cat = NearestNeighborPCs(C,idx,pc_fdmat,meanfd_fdmat,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter);
            else
                percent_in_cat = NearestNeighborBasic(C,idx,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,basis);
            end
            % ------ Compute probability of death
            if weightbydistance
                % Weight risk by distance to different clusters
                weights_per_PC = sum(D,2)./D;
                normalized_wt_per_PC = weights_per_PC./sum(weights_per_PC,2);
                prob_of_outcome = normalized_wt_per_PC*percent_in_cat(:,1)/100;
            else
                % Just grab cluster risk
                prob_of_outcome = percent_in_cat(idx);
            end
        case 3 % Use Gaussian Mixture Model
            options = statset('MaxIter',300);
            nlogLall = nan*ones(8,1);
            pvalues = zeros(8,4);
            clusters = zeros(length(harmscr),7);
            GMMBIC = zeros(8,1);
            for neighborhoods=1:8
                rng(4); % this sets the starting parameters the same way each time for reproducibility
                figure()
                GMM = fitgmdist(harmscr,neighborhoods,'Options',options);
                [clusterGMM,nlogL] = cluster(GMM,harmscr);
                if neighborhoods>1
                    clusters(:,neighborhoods-1) = clusterGMM;
                    P(neighborhoods-1).clusterprob = posterior(GMM,harmscr);
                end
                if usePCs
                    NearestNeighborPCs(GMM.mu,clusterGMM,pc_fdmat,meanfd_fdmat,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,limofinterest(v,:),[]);
                else
                    NearestNeighborBasic(GMM.mu,clusterGMM,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,basis);
                end
                for bn = 1:neighborhoods
                    % Find number of babies in each category
                    babies_in_hood = sum(clusterGMM==bn); % babies in neighborhood
                    for c=1:length(unique(category))
                        babies_in_cat(bn,c) = sum(and((clusterGMM==bn),category==c)); % babies in category
                        percent_in_cat(bn,c) = round(babies_in_cat(bn,c)/babies_in_hood*100);
                    end
                end
%                 addpath('X:\Amanda\FuncDataAnalysis\Hero')
%                 pvalues(neighborhoods,:) = makePriorProbabilityTable(goodlist,xtselected,a,u1,bwt,ttd,clusterGMM);
                nlogLall(neighborhoods) = nlogL;
            end
            bestnumhoods = find(nlogLall==min(nlogLall));
%             bestnumhoods = 3;
            rng(4);
            GMM = fitgmdist(harmscr,bestnumhoods,'Options',options);
            [clusterGMM,nlogL] = cluster(GMM,harmscr);
            figure()
            if usePCs
                NearestNeighborPCs(GMM.mu,clusterGMM,pc_fdmat,meanfd_fdmat,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,limofinterest(v,:),fdstruct(v).dayfd);
            else
                NearestNeighborBasic(GMM.mu,clusterGMM,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,basis);
            end
            for bn = 1:bestnumhoods
                % Find number of babies in each category
                babies_in_hood = sum(clusterGMM==bn); % babies in neighborhood
                for c=1:length(unique(category))
                    babies_in_cat(bn,c) = sum(and((clusterGMM==bn),category==c)); % babies in category
                    percent_in_cat(bn,c) = round(babies_in_cat(bn,c)/babies_in_hood*100);
                end
            end
            
            figure(25)
            plot(1:8,nlogLall)
            hold on
            xlabel('Number of Neighborhoods')
            ylabel('nLogL')
            hold off
            
%             figure(26)
%             plot(1:8,pvalues(:,3));
%             hold on
%             xlabel('Number of Neighborhoods')
%             ylabel('p value for 30 day mortality')
%             hold off
%             
%             figure(27)
%             plot(1:8,pvalues(:,4));
%             hold on
%             xlabel('Number of Neighborhoods')
%             ylabel('p value for 7 day mortality')
%             hold off
            
            if weightbydistance
                P = posterior(GMM,harmscr);
                prob_of_outcome = P*percent_in_cat(:,1)/100;
            else
                prob_of_outcome = percent_in_cat(clusterGMM);
            end
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
            [resultsLWEA,resultsLWGP] = LWEA_and_LWGP(allclusters,category);
            C_LWEA = findclustercentroids(resultsLWEA,harmscr);
            C_LWGP = findclustercentroids(resultsLWGP,harmscr); 
            if usePCs
                for clusts = 1:11
                    figure()
                    percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),pc_fdmat,meanfd_fdmat,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,limofinterest(v,:),fdstruct(v).dayfd);
                end
            else
                for clusts = 1:11
                    figure()
                    percent_in_cat = NearestNeighborBasic(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter,basis);
                end
            end
            prob_of_outcome = percent_in_cat(resultsLWEA(:,clusts));
            clusters = resultsLWEA;
    end
    
    % ------------------- Make ROC Curve ------------------------------
    figure();
    [X,Y,T,AUC] = perfcurve(category,prob_of_outcome,1,'Prior','empirical');
    plot(X,Y)
    xlabel('False positive rate'); ylabel('True positive rate');
    title(['AUC: ' num2str(round(AUC,2))])
    fprintf(['AUC: ' num2str(round(AUC,2)) '\n'])
    
    % ------------------- Score One Method -----------------------------
    
    qtyofoutcome1 = length(category1);
    [~,I] = sort(prob_of_outcome,'descend');
    percentcorrect = sum(category(I(1:qtyofoutcome1))==1)/qtyofoutcome1;
    incidentrate = qtyofoutcome1/length(category);
    improvement_over_incident_rate = percentcorrect/incidentrate;
    fprintf(['Percent correct: ' num2str(round(percentcorrect*100)) '%% \n']);
    fprintf(['Improvement over incident rate: ' num2str(improvement_over_incident_rate) ' \n'])
    
    % -------------- Find Confidence Intervals ----------------------
    
    if model == 3 && (dataset<3 || dataset==5)
        % Compute mapping from data y to coefficients in c
        basismat = eval_basis(daytime, basis);
        y2cMap = (basismat'*basismat)\basismat';
    
        % Compute residual matrix and get covariance of residuals
        yhatmat  = eval_fd(daytime, yhatfdobj);
        ymat     = eval_fd(daytime, fdstruct(v).dayfd);
        tempresmat = ymat(:,1:size(gooddata,2)) - yhatmat(:,1:size(gooddata,2));
        SigmaE   = cov(tempresmat');

        % Plot covariance surface for errors
        % figure()
        % contour(SigmaE)
        % hold on
        % plot(dayrange,dayrange,'--')
        % hold off
        % colorbar

        % Plot standard deviation of errors
        % figure()
        % stddevE = sqrt(diag(SigmaE));
        % plot(daytime, stddevE, '-')
        % xlabel('\fontsize{19} Day')
        % ylabel(['\fontsize{19} ' vname{v} ' units'])

        % Repeat regression, this time outputting results for confidence intervals
        stderrStruct = fRegress_stderr(fRegressStruct, y2cMap, SigmaE);
        betastderrcell = stderrStruct.betastderr;

        % Plot regression functions
        % for j=1:p
        %     figure()
        %     plot(daytime, eval_fd(daytime, betastderrcell{j}))
        %     title(['\fontsize{16} ',group_names(j,:)])
        % end

        % Plot regression functions with confidence limits
        for j=1:p
            figure()
            plotbeta(betaestcell{j}, betastderrcell{j}, daytime')
            xlabel('\fontsize{19} Day')
            ylabel(['\fontsize{19} ' vname{column}])
            title(['\fontsize{16} ',group_names(j,:)])
        end

        % Plot predicted functions with shaded error bars
        plotprederrbar(yhatfdobj, betastderrcell, daytime')
        % legend([group_names(2:end,:);group_names(1,:)]) % The "All" group is at the end
        if dataset == 1
            legend(group_names(2:end,:));
            if ~subtractoffmean
                ylim([limofinterest(v,1),limofinterest(v,2)]);
            end
        elseif dataset == 2 || dataset==5
            labelordervector = unique(category,'stable')'; % NOTE: This only works if the indices are 1,2,3,etc.
            groupnamesnoall = group_names(2:end,:);
            legend(groupnamesnoall(labelordervector,:));
            ylim([limofinterest(v,1),limofinterest(v,2)]);
        end
    end
    
    % ------------ Create regression model comparison -------------------
    
    if alg ==4
        daysofage = pdate-bd(pnum);
        windowvdata = squeeze(vdata_interp_all(start_tt:end_tt,v,:));
        lastvalues(:,v) = windowvdata(end,:);

        priorprob = findBW_DOA_priorprob(bwt(pnum(goodindices(v,:))),daysofage(goodindices(v,:)),category,id(pnum(goodindices(v,:))));

        n = sum(goodindices);

        regmodel(1).params = bwt(pnum(goodindices(v,:)));
        regmodel(1).name = 'BWT';

        regmodel(2).params = [bwt(pnum(goodindices(v,:))), clusters(:,1)==1];
        regmodel(2).name = 'BWT + 2 Neighborhoods';

        regmodel(3).params = [bwt(pnum(goodindices(v,:))), clusters(:,2)==1,clusters(:,2)==2];
        regmodel(3).name = 'BWT + 3 Neighborhoods';

        regmodel(4).params = [bwt(pnum(goodindices(v,:))), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
        regmodel(4).name = 'BWT + 4 Neighborhoods';

        regmodel(5).params = [bwt(pnum(goodindices(v,:))), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
        regmodel(5).name = 'BWT + 5 Neighborhoods';

        regmodel(6).params = [bwt(pnum(goodindices(v,:))), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
        regmodel(6).name = 'BWT + 6 Neighborhoods';

        regmodel(7).params = priorprob; %[bwt(pnum(goodindices(v,:))), daysofage(goodindices(v,:))];
        regmodel(7).name = 'BWT/DOA';

        regmodel(8).params = [priorprob, lastvalues(goodindices(v,:),v)]; %[bwt(pnum(goodindices(v,:))), daysofage(goodindices(v,:)), lastvalues(goodindices(v,:),v)];
        regmodel(8).name = 'BWT/DOA + Last Value';

        regmodel(9).params = [priorprob, clusters(:,1)==1];
        regmodel(9).name = 'BWT/DOA + 2 Neighborhoods';

        regmodel(10).params = [priorprob, clusters(:,2)==1,clusters(:,2)==2];
        regmodel(10).name = 'BWT/DOA + 3 Neighborhoods';

        regmodel(11).params = [priorprob, clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
        regmodel(11).name = 'BWT/DOA + 4 Neighborhoods';

        regmodel(12).params = [priorprob, clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
        regmodel(12).name = 'BWT/DOA + 5 Neighborhoods';

        regmodel(13).params = [priorprob, clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
        regmodel(13).name = 'BWT/DOA + 6 Neighborhoods';

        regmodel(14).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v)];
        regmodel(14).name = 'BWT + Last Value';

        regmodel(15).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,1)==1];
        regmodel(15).name = 'BWT + Last Value + 2 Neighborhoods';

        regmodel(16).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,2)==1,clusters(:,2)==2];
        regmodel(16).name = 'BWT + Last Value + 3 Neighborhoods';

        regmodel(17).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
        regmodel(17).name = 'BWT + Last Value + 4 Neighborhoods';

        regmodel(18).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
        regmodel(18).name = 'BWT + Last Value + 5 Neighborhoods';

        regmodel(19).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
        regmodel(19).name = 'BWT + Last Value + 6 Neighborhoods';

        regmodel(20).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,1)==1];
        regmodel(20).name = 'BWT/DOA + Last Value + 2 Neighborhoods';

        regmodel(21).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,2)==1,clusters(:,2)==2];
        regmodel(21).name = 'BWT/DOA + Last Value + 3 Neighborhoods';

        regmodel(22).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
        regmodel(22).name = 'BWT/DOA + Last Value + 4 Neighborhoods';

        regmodel(23).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
        regmodel(23).name = 'BWT/DOA + Last Value + 5 Neighborhoods';

        regmodel(24).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
        regmodel(24).name = 'BWT/DOA + Last Value + 6 Neighborhoods';

        regmodel(25).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
        regmodel(25).name = 'BWT/DOA + Last Value + 7 Neighborhoods';

        regmodel(26).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
        regmodel(26).name = 'BWT/DOA + Last Value + 8 Neighborhoods';

        regmodel(27).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
        regmodel(27).name = 'BWT/DOA + Last Value + 9 Neighborhoods';

        regmodel(28).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
        regmodel(28).name = 'BWT/DOA + Last Value + 10 Neighborhoods';
        
        regmodel(29).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
        regmodel(29).name = 'BWT/DOA + Last Value + 11 Neighborhoods';
        
        regmodel(30).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
        regmodel(30).name = 'BWT/DOA + Last Value + 12 Neighborhoods';
        
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
        
        regmodel(37).params = lastvalues(goodindices(v,:),v);
        regmodel(37).name = 'Last Value';
        
        regmodel(38).params = [lastvalues(goodindices(v,:),v), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
        regmodel(38).name = 'Last Value + 6 Neighborhoods';
        
    %     
    %     regmodel(20).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), P(1).clusterprob(:,1)];
    %     regmodel(20).name = 'BWT + Last Value + 2 Neigh. Prob.';
    %     
    %     regmodel(21).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), P(2).clusterprob(:,1),P(2).clusterprob(:,2)];
    %     regmodel(21).name = 'BWT + Last Value + 3 Neigh. Prob.';
    %     
    %     regmodel(22).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), P(3).clusterprob(:,1),P(3).clusterprob(:,2),P(3).clusterprob(:,3)];
    %     regmodel(22).name = 'BWT + Last Value + 4 Neigh. Prob.';
    % 
    %     regmodel(23).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), P(4).clusterprob(:,1),P(4).clusterprob(:,2),P(4).clusterprob(:,3),P(4).clusterprob(:,4)];
    %     regmodel(23).name = 'BWT + Last Value + 5 Neigh. Prob.';
    %     
    %     regmodel(24).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), P(5).clusterprob(:,1),P(5).clusterprob(:,2),P(5).clusterprob(:,3),P(5).clusterprob(:,4),P(5).clusterprob(:,5)];
    %     regmodel(24).name = 'BWT + Last Value + 6 Neigh. Prob.';

        alld = zeros(length(regmodel),1);
        alldev = zeros(length(regmodel),1);
        allAIC = zeros(length(regmodel),1);
        allBIC = zeros(length(regmodel),1);
        allBriar = zeros(length(regmodel),1);
        allAUC = zeros(length(regmodel),1);
        allpihat = zeros(n,length(regmodel));
        allrownames = {};

        for m=1:length(regmodel)
            [B,d,dev,AIC,BIC,Briar,AUC,pihat,X,Y,stats] = modelcomparison(regmodel(m).params,category,n);
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
            allpihat(:,m) = pihat;
            allrownames{m} =  regmodel(m).name;
        end

        alldev = round(alldev,1);
        allAIC = round(allAIC,1);
        allBIC = round(allBIC,1);
        allBriar = round(allBriar,5);
        allAUC = round(allAUC,2);

        Table = table(alld,alldev,allAIC,allBIC,allBriar,allAUC,'VariableNames',{'DOF','Dev','AIC','BIC','Briar','AUC'},'rowNames',allrownames);
        Table

        priorprobmodel = 7;
        basemodel = 8;
        bettermodel = 13;
        bestmodel = 24;
        
%         % Wald Chi-Square Test
%         Y = allB(bestmodel).B;
%         Y0 = allB(basemodel).B;
%         model = arima(2,0,0);
%         [fit,V] = estimate(model,Y,'Y0',Y0);
%         r = fit.AR{2};
        

        figure()
        scatter(allpihat(category==2,basemodel),allpihat(category==2,bestmodel),[],'b.')
        hold on
        scatter(allpihat(category==1,basemodel),allpihat(category==1,bestmodel),[],'r.')
        xlabel(regmodel(basemodel).name)
        ylabel(regmodel(bestmodel).name)
        legend(group_names(3,:),group_names(2,:))
        axis square

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
end

function [B,d,dev,AIC,BIC,Briar,AUC,pihat,X,Y,stats] = modelcomparison(modelparams,category,n)
   % Logistic Regression
%     [B,dev,stats] = mnrfit(modelparams,category);
%     pihat = mnrval(B,modelparams);
%     pihat = pihat(:,2);
%     RSS = nansum(stats.resid(:,2).^2);
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