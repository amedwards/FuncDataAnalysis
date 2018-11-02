%  --------------- Load Data and Set Up Workspace  ----------------

% Add Ramsay's scripts to the path
addpath ('X:\Amanda\FuncDataAnalysis\Ramsay')
clear

% Choose which dataset to load
load('X:\Amanda\FuncDataAnalysis\UVA_CU_WUSTL_FPCA_Stats_day.mat')

% Load Pete IDs
load('X:\Amanda\HistogramsAndHeatmaps\peteids.mat')

% Load Karen IDs
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\KarenIDs.mat')
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\KarenPercentages.mat')

% Load Columbia Aim 3 Eligible
load('X:\Amanda\Columbia\Aim3Eligible.mat')

% Only Run Analysis on Variables of Interest
varofinterest = {'Mean HR','Mean SPO2-%'};
limofinterest = [100 180; 80 100];
[varlog,varnums] = ismember(vname,varofinterest);
nv = sum(varlog);

% Choose Parameters
weeks = 6;
maxday = 7*weeks; 
dayrange = [1,maxday];

dayrange28 = [1,28];
daytime = 1:maxday;
daytime28 = 1:28;

nbasis = ceil(maxday/2);
thresh = 0.5;
nharm  = 4;
grouping = 1; % set to 1 to split by site, set to 2 to split by gender, 3 to split by ELBW/VLBW, 4 to split by ega

% ------------------- Create Basis Functions -----------------------

Lbasis  = create_constant_basis(dayrange);  %  create a constant basis
Mbasis = create_monomial_basis(dayrange,nbasis); % create a monomial basis
Bbasis = create_bspline_basis(dayrange,nbasis,5); % create a bspline basis - needed to up the order to 5 for smoothing penalty requirement

% Choose which basis function to use
basis = Bbasis;

%  ----------  Set up the harmonic acceleration operator  --------------

Lcoef   = [0,(4/length(daytime))^2,0];    %  set up three coefficients
wfd     = fd(Lcoef,Lbasis);      % define an FD object for weight functions
wfdcell = fd2cell(wfd);          % convert the FD object to a cell object
harmaccelLfd = Lfd(3, wfdcell);  %  define the operator object


%  -------------  Identify which data we want to use  ---------------

% Keep only the 502 Pete babies and the Aim 3 included CU IDs
if exist('peteids','var')
    peterows = ones(n,1);
%     peterows(inst==1) = ismember(pid(inst==1),peteids); % ex. data for baby 1089 is stored in row 1081 out of 6755
    KarenIDs = KarenIDs(KarenPercentages>=0.25); % Keep data from babies who have at least 7 days of data in the first 28 days
    peterows(inst==1) = ismember(pid(inst==1),KarenIDs); % ex. data for baby 1089 is stored in row 1081 out of 6755
    goodCUids = EligibilityforAim3paper.DEIDMRN(EligibilityforAim3paper.Aim3Inclusion==1);
    peterows(inst==2) = ismember(pid(inst==2),goodCUids);
else
    peterows = ones(n,1);
end

% Initialize empty arrays and structs
percentavail = zeros(n,nv);
percentavail28 = zeros(n,nv);
vdata_interp = zeros(maxday,nv,size(dataday,3));
fdstruct = struct();
pcastruct = struct();

for v=1:nv
    column = find(varnums==v);
    % --------------------- Data Preparation ----------------------------
    % Pull the data for variable v
    vdata = squeeze(dataday(daytime,column,:));
    vdata28 = squeeze(dataday(daytime28,column,:));
    
    % Plot the number of babies with data available each day
    % figure(); plot(sum(~isnan(vdata),2))
    % xlabel('Day of age'); ylabel('Number of babies from which we have data that day'); title(vname{v})
    
    % How much data does each baby have?
    percentavail(:,v) = sum(~isnan(vdata),1)/length(daytime)';
    percentavail28(:,v) = sum(~isnan(vdata28),1)/28';
    
    % Interpolate/Extrapolate for missing data
    vdata_interp(:,v,:) = fillmissing(vdata,'linear','SamplePoints',daytime,'EndValues','nearest');
    
    % Find babies which have a certain percentage of data and are Pete babies
    good28 = percentavail28(:,v)>=0.25; % need at least 7 days of data in first 28 days
    gooddata = percentavail(:,v)>=thresh;
    goodindices = and(gooddata,peterows);
    goodindices = and(goodindices,good28);
    
    sprintf(['UVA   Babies Included: ' num2str(sum(goodindices(inst==1)))])
    sprintf(['CU    Babies Included: ' num2str(sum(goodindices(inst==2)))])
    sprintf(['WUSTL Babies Included: ' num2str(sum(goodindices(inst==3)))])
    
    gooddata = squeeze(vdata_interp(:,v,goodindices));
    
    % ------------------ Fit Basis Functions ---------------------------
    
    % Fit basis functions for each baby
    fdstruct(v).dayfd = smooth_basis(daytime, gooddata, basis);
    dayfd_fdnames = {'Day','ID',vname{v}};
    fdstruct(v).dayfd = putnames(fdstruct(v).dayfd, dayfd_fdnames);
    
    % Plot individual basis curves and values
    % figure(); plotfit_fd(gooddata, daytime, fdstruct(v).dayfd)
    
    % Plot all basis curves at once
    figure(); plot(fdstruct(v).dayfd)
    
    % ----------------- Determine Level of Smoothing ------------------
    % Choose level of smoothing using the generalized cross-validation 
    %          criterion with smoothing function smooth_basis.

    % Set up range of smoothing parameters in log_10 units
    loglam = (-5:1)';
    nlam   = length(loglam);
    dfsave  = zeros(nlam,1);
    gcvsave = zeros(nlam,1);

    % Loop through smoothing parameters
    for ilam=1:length(loglam)
        lambda = 10^loglam(ilam);
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
    
    % -------------------- Smooth our data -----------------------------

    %  Do final smooth with minimum GCV value
    lambda = 1;  %  minimum GCV estimate
    fdParobj = fdPar(basis, harmaccelLfd, lambda);
    fdstruct(v).dayfd = smooth_basis(daytime, gooddata, fdParobj);
    
    % Plot data and fit
    % subplot(1,1,1)
    % plotfit_fd(gooddata, daytime, fdstruct(v).dayfd, vname{v})
           
    % ------------------------- PCA ------------------------------------
    
    pcastruct(v).daypcastr = pca_fd(fdstruct(v).dayfd,nharm);
    figure(); plot_pca_fd(pcastruct(v).daypcastr,1)
    
    % -------------------- Functional ANOVA -------------------------
    % Names for groups
    switch grouping
        case 1
            group_names = ['All  '; 'UVA  ';'CU   ';'WUSTL'];
            % Site Indices
            category = inst(goodindices);
            category1 = find(category==1)'; % UVA
            category2 = find(category==2)'; % CU
            category3 = find(category==3)'; % WUSTL
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy contraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
            zmat(category3,4) = 1;
        case 2
            group_names = ['All   ';'Female'; 'Male  '];
            % Gender Indices
            category = pgen(goodindices);
            category1 = find(category==1)'; % Female
            category2 = find(category==2)'; % Male
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy contraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 3
            group_names = ['All ';'ELBW'; 'VLBW'];
            % Gender Indices
            category = pbw(goodindices);
            category1 = find(category<1000)'; % ELBW
            category2 = find(category>=1000)'; % VLBW non ELBW
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site/weight. Add a dummy contraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 4
            group_names = ['All        '; '<27 Weeks  ';'27-30 Weeks';'31-34 Weeks';'>34 Weeks  '];
            % Site Indices
            category = pega(goodindices);
            category1 = find(category<27)'; % <27 Weeks
            category2 = find(category>=27&category<30)'; % 27-30 Weeks
            category3 = find(category>=30&category<34)'; % 31-34 Weeks
            category4 = find(category>=34)'; % >34 Weeks
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy contraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
            zmat(category3,4) = 1;
            zmat(category4,5) = 1;
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

    % Set up the basis for the regression functions
    nbetabasis = 10;
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
    
    
    % -------------- Find Confidence Intervals ----------------------
    
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
    legend(group_names(2:end,:));
    ylim([limofinterest(v,1),limofinterest(v,2)]);
end
