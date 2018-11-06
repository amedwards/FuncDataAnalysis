%  --------------- Load Data and Set Up Workspace  ----------------

% Add Ramsay's scripts to the path
addpath ('X:\Amanda\FuncDataAnalysis\Ramsay')
clear
close all

% Choose which dataset to load
dataset = 4; % 1 = Sepsis/NEC, 2 = HERO Sepsis Cases, 3 = All HERO cultures, 4 = Control Hero Only, exclude bad pg's

% 1-Mean model, 2-Slope model, 3-Bspline model, 4-Last Hero Value
model = 3;

% Choose whether you want to subtract off the PRECEEDING mean
subtractoffmean = 0;

% Choose whether to use 1: Principal Components or 0:Raw output from basis functions for clustering
usePCs = 0;

% Subtract (NOT preceeding) mean from data, then add to array to help with prediction
colofmeans = 0;

% Add in birthweight
add_bwt = 0;

% Algorithm: 1 = logistic regression, 2 = kmeans clustering, 3 = Gaussian Mixture Model
alg = 2;

% Select Grouping: 1 - split by site, 2 - split by gender, 3 -split by ELBW/VLBW, 4 - split by ega, 5 - split by control vs Display for Hero, 6 - split by 30 day survival, 7 - split by pos/neg blood culture, 8 - negsep
grouping = 8; 

% Choose number of neighborhoods
neighborhoods = 6; % usually 6

% Choose our data window
daysbefore = -5;
daysafter = 0;

% Percent threshold for the amount of data that is required for inclusion
thresh = 0.85;

% Choose whether to weight risk score by distance
weightbydistance = 1;

% Turn on smoothing
smoothingon = 0;





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
    toinclude = pg>0 & control;
    
    % Only Run Analysis on Variables of Interest
    varofinterest = {'hero'};
    vname = varofinterest;
    limofinterest = [-1 7];
    varnums = 1;
    nv = 1;
    
    % Set basis parameters
    reductionfactor = 6;
    lambdabase = 0.1;
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
end

% ------------------- Create Basis Functions -----------------------

Lbasis  = create_constant_basis(dayrange);  %  create a constant basis
Mbasis = create_monomial_basis(dayrange,nbasis); % create a monomial basis
if model==3
    Bbasis = create_bspline_basis(dayrange,nbasis,5); % create a bspline basis - needed to up the order to 5 for smoothing penalty requirement
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
    goodindices = gooddata;
    
    % Subtract off the data mean
    datamean(v,:) = nanmean(squeeze(vdata_interp(:,v,:)));
    if colofmeans
        vdata_interp(:,v,:) = squeeze(vdata_interp(:,v,:))-datamean(v,:);
    end
    
    if dataset == 1
        fprintf(['UVA   Babies Included: ' num2str(sum(goodindices(inst==1))) '\n'])
        fprintf(['CU    Babies Included: ' num2str(sum(goodindices(inst==2))) '\n'])
    elseif dataset == 2
        fprintf(['Babies Included: ' num2str(sum(goodindices)) '\n'])
    end
    
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
    figure(); 
    [meanfd_fdmat,pc_fdmat] = plot_pca_fd(pcastruct(v).daypcastr,1); % store the principal components that we plot
    
    % -------------------- Functional ANOVA -------------------------
    % Names for groups
    switch grouping
        case 1
            group_names = ['All  '; 'UVA  ';'CU   '];
            % Site Indices
            category = inst(goodindices);
            category1 = find(category==1)'; % UVA
            category2 = find(category==2)'; % CU
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy constraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 2
            group_names = ['All   ';'Female'; 'Male  '];
            % Gender Indices
            category = pgen(goodindices);
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
            category = pbw(goodindices);
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
            category = pega(goodindices);
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
            % Gender Indices
            category = gg(goodindices);
            category1 = find(category==1)'; % Female
            category2 = find(category==2)'; % Male
            % Set up a design matrix having a column for the grand mean and a 
            % column for each gender/site. Add a dummy constraint observation.
            p = size(group_names,1);
            zmat = zeros(size(gooddata,2),p);
            zmat(:,1) = 1;
            zmat(category1,2) = 1;
            zmat(category2,3) = 1;
        case 6
            group_names = ['All         '; 'Death in 30d'; 'Survival    '];
            diedin30days = ddate(pnum)<pdate+30;
            category = double(diedin30days(goodindices));
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
            group_names = ['All     ';'Positive';'Negative'];
            category = pg(goodindices);
            %Switch the category labels
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
            category = double(negsep(goodindices));
            category1 = find(category==1); % Negsep
            category2 = find(category==0); % Not negsep (Clinsep + Sep)
            category(category==0) = 2;
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
        harmscr = [harmscr, datamean(goodindices)'];
    end
    
    if add_bwt
        harmscr = [harmscr, bwt(pnum(goodindices))];
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
            GMM = fitgmdist(harmscr,neighborhoods,'Options',options);
            clusterGMM = cluster(GMM,harmscr);
            if usePCs
                NearestNeighborPCs(GMM.mu,clusterGMM,pc_fdmat,meanfd_fdmat,daytime,vname{column},category,group_names(2:end,:),daysbefore,daysafter);
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
            if weightbydistance
                P = posterior(GMM,harmscr);
                prob_of_outcome = P*percent_in_cat(:,1)/100;
            else
                prob_of_outcome = percent_in_cat(clusterGMM);
            end
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
    fprintf(['Score One Method: Percent correct: ' num2str(round(percentcorrect*100)) '%% \n']);
    fprintf(['Score One Method: Improvement over incident rate: ' num2str(improvement_over_incident_rate) ' \n'])
    
    % -------------- Find Confidence Intervals ----------------------
    
    if model == 3 && dataset<3
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
        elseif dataset == 2
            labelordervector = unique(category,'stable')'; % NOTE: This only works if the indices are 1,2,3,etc.
            groupnamesnoall = group_names(2:end,:);
            legend(groupnamesnoall(labelordervector,:));
            ylim([limofinterest(v,1),limofinterest(v,2)]);
        end
    end
end
