load('load_data_for_loop_rng_value')
addpath('EnsembleClustering')
for e=1:100
    rng(e)
    randk = randi([mink maxk],niterk,1); % generate 100 random integers between mink and maxk
    allclusters = zeros(objects,niterk);
    for a=1:length(randk)
        [idx] = kmeans(harmscr,randk(a));
        allclusters(:,a) = idx;
    end
    [resultsLWEA,resultsLWGP] = LWEA_and_LWGP(allclusters,category,clusteriterations);
    C_LWEA = findclustercentroids(resultsLWEA,harmscr);
    C_LWGP = findclustercentroids(resultsLWGP,harmscr); 
    if ~multisignalclust
        if usePCs
            for clusts = 1:11
                figure()
                if exist('pnum')
                    percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),squeeze(pc_fdmat(v,:,:)),meanfd_fdmat(v,:),daytime,vname{column},category,c_struct,daysbefore,daysafter,limofinterest(v,:),fdstruct(v).dayfd,id(pnum(goodindices(v,:))));
                else
                    percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),squeeze(pc_fdmat(v,:,:)),meanfd_fdmat(v,:),daytime,vname{column},category,c_struct,daysbefore,daysafter,limofinterest(v,:),fdstruct(v).dayfd);
                end
            end
        else
            for clusts = 1:11
                figure()
                percent_in_cat = NearestNeighborBasic(C_LWEA(clusts).centroid,resultsLWEA(:,clusts),daytime,vname{column},category,c_struct,daysbefore,daysafter,basis);
            end
        end
    else
        colors = [1 0 0;0.5430 0 0;0 0 .7;0 0.5 0;0.25 0.25 0;0 0.75 0.1;0.1 0.4 0.2];
        if v==nv
            for clusts=1:11
                figure()
                for b=1:nv
                    column = find(varnums==b);
                    color = colors(b,:);
                    if exist('pnum')
                        percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid(:,b*4-3:b*4),resultsLWEA(:,clusts),squeeze(pc_fdmat(b,:,:)),meanfd_fdmat(b,:),daytime,{vname{varlog}}',category,c_struct,daysbefore,daysafter,limofinterest(b,:),fdstruct(b).dayfd,id(pnum(goodindices(v,:))),color);
                    else
                        percent_in_cat = NearestNeighborPCs(C_LWEA(clusts).centroid(:,b*4-3:b*4),resultsLWEA(:,clusts),squeeze(pc_fdmat(b,:,:)),meanfd_fdmat(b,:),daytime,{vname{varlog}}',category,c_struct,daysbefore,daysafter,limofinterest(b,:),fdstruct(b).dayfd,[],color);
                    end
                    hold on
                end
            end
        end
    end

    prob_of_outcome = percent_in_cat(resultsLWEA(:,clusts));
    clusters = resultsLWEA;
    
    
    
    
    
    
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

        % Repeat regression, this time outputting results for confidence intervals
        stderrStruct = fRegress_stderr(fRegressStruct, y2cMap, SigmaE);
        betastderrcell = stderrStruct.betastderr;

%         % Plot regression functions with confidence limits
%         for j=1:p
%             figure()
%             plotbeta(betaestcell{j}, betastderrcell{j}, daytime')
%             xlabel('\fontsize{19} Day')
%             ylabel(['\fontsize{19} ' vname{column}])
%             title(['\fontsize{16} ',group_names(j,:)])
%         end

%         % Plot predicted functions with shaded error bars
%         plotprederrbar(yhatfdobj, betastderrcell, daytime')
%         % legend([group_names(2:end,:);group_names(1,:)]) % The "All" group is at the end
%         if dataset == 1
%             legend(group_names(2:end,:));
%             if ~subtractoffmean
%                 ylim([limofinterest(v,1),limofinterest(v,2)]);
%             end
%         elseif dataset == 2 || dataset==5
%             labelordervector = unique(category,'stable')'; % NOTE: This only works if the indices are 1,2,3,etc.
%             groupnamesnoall = group_names(2:end,:);
%             legend(groupnamesnoall(labelordervector,:));
%             ylim([limofinterest(v,1),limofinterest(v,2)]);
%         end
    end
    
    % ------------ Create regression model comparison -------------------

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

    regmodel(7).params = [bwt(pnum(goodindices(v,:))),clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(7).name = 'BWT + 7 Neighborhoods';

    regmodel(8).params = [bwt(pnum(goodindices(v,:))),clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(8).name = 'BWT + 8 Neighborhoods';

    regmodel(9).params = [bwt(pnum(goodindices(v,:))),clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(9).name = 'BWT + 9 Neighborhoods';

    regmodel(10).params = [bwt(pnum(goodindices(v,:))), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(10).name = 'BWT + 10 Neighborhoods';

    regmodel(11).params = [bwt(pnum(goodindices(v,:))), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(11).name = 'BWT + 11 Neighborhoods';

    regmodel(12).params = [bwt(pnum(goodindices(v,:))), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(12).name = 'BWT + 12 Neighborhoods';

    regmodel(13).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v)];
    regmodel(13).name = 'BWT + Last Value';

    regmodel(14).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,1)==1];
    regmodel(14).name = 'BWT + Last Value + 2 Neighborhoods';

    regmodel(15).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,2)==1,clusters(:,2)==2];
    regmodel(15).name = 'BWT + Last Value + 3 Neighborhoods';

    regmodel(16).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
    regmodel(16).name = 'BWT + Last Value + 4 Neighborhoods';

    regmodel(17).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
    regmodel(17).name = 'BWT + Last Value + 5 Neighborhoods';

    regmodel(18).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(18).name = 'BWT + Last Value + 6 Neighborhoods';

    regmodel(19).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(19).name = 'BWT + Last Value + 7 Neighborhoods';

    regmodel(20).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(20).name = 'BWT + Last Value + 8 Neighborhoods';

    regmodel(21).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(21).name = 'BWT + Last Value + 9 Neighborhoods';

    regmodel(22).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(22).name = 'BWT + Last Value + 10 Neighborhoods';

    regmodel(23).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(23).name = 'BWT + Last Value + 11 Neighborhoods';

    regmodel(24).params = [bwt(pnum(goodindices(v,:))), lastvalues(goodindices(v,:),v), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(24).name = 'BWT + Last Value + 12 Neighborhoods';

    regmodel(25).params = priorprob; %[bwt(pnum(goodindices(v,:))), daysofage(goodindices(v,:))];
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

    regmodel(37).params = [priorprob, lastvalues(goodindices(v,:),v)]; %[bwt(pnum(goodindices(v,:))), daysofage(goodindices(v,:)), lastvalues(goodindices(v,:),v)];
    regmodel(37).name = 'BWT/DOA + Last Value';

    regmodel(38).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,1)==1];
    regmodel(38).name = 'BWT/DOA + Last Value + 2 Neighborhoods';

    regmodel(39).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,2)==1,clusters(:,2)==2];
    regmodel(39).name = 'BWT/DOA + Last Value + 3 Neighborhoods';

    regmodel(40).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,3)==1,clusters(:,3)==2,clusters(:,3)==3];
    regmodel(40).name = 'BWT/DOA + Last Value + 4 Neighborhoods';

    regmodel(41).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,4)==1,clusters(:,4)==2,clusters(:,4)==3,clusters(:,4)==4];
    regmodel(41).name = 'BWT/DOA + Last Value + 5 Neighborhoods';

    regmodel(42).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(42).name = 'BWT/DOA + Last Value + 6 Neighborhoods';

    regmodel(43).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,6)==1,clusters(:,6)==2,clusters(:,6)==3,clusters(:,6)==4,clusters(:,6)==5,clusters(:,6)==6];
    regmodel(43).name = 'BWT/DOA + Last Value + 7 Neighborhoods';

    regmodel(44).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,7)==1,clusters(:,7)==2,clusters(:,7)==3,clusters(:,7)==4,clusters(:,7)==5,clusters(:,7)==6,clusters(:,7)==7];
    regmodel(44).name = 'BWT/DOA + Last Value + 8 Neighborhoods';

    regmodel(45).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,8)==1,clusters(:,8)==2,clusters(:,8)==3,clusters(:,8)==4,clusters(:,8)==5,clusters(:,8)==6,clusters(:,8)==7,clusters(:,8)==8];
    regmodel(45).name = 'BWT/DOA + Last Value + 9 Neighborhoods';

    regmodel(46).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,9)==1,clusters(:,9)==2,clusters(:,9)==3,clusters(:,9)==4,clusters(:,9)==5,clusters(:,9)==6,clusters(:,9)==7,clusters(:,9)==8,clusters(:,9)==9];
    regmodel(46).name = 'BWT/DOA + Last Value + 10 Neighborhoods';

    regmodel(47).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,10)==1,clusters(:,10)==2,clusters(:,10)==3,clusters(:,10)==4,clusters(:,10)==5,clusters(:,10)==6,clusters(:,10)==7,clusters(:,10)==8,clusters(:,10)==9,clusters(:,10)==10];
    regmodel(47).name = 'BWT/DOA + Last Value + 11 Neighborhoods';

    regmodel(48).params = [priorprob, lastvalues(goodindices(v,:),v), clusters(:,11)==1,clusters(:,11)==2,clusters(:,11)==3,clusters(:,11)==4,clusters(:,11)==5,clusters(:,11)==6,clusters(:,11)==7,clusters(:,11)==8,clusters(:,11)==9,clusters(:,11)==10,clusters(:,11)==11];
    regmodel(48).name = 'BWT/DOA + Last Value + 12 Neighborhoods';

    regmodel(49).params = lastvalues(goodindices(v,:),v);
    regmodel(49).name = 'Last Value';

    regmodel(50).params = [lastvalues(goodindices(v,:),v),clusters(:,5)==1,clusters(:,5)==2,clusters(:,5)==3,clusters(:,5)==4,clusters(:,5)==5];
    regmodel(50).name = 'Last Value + 6 Neighborhoods';

    alld = zeros(length(regmodel),1);
    alldev = zeros(length(regmodel),1);
    allAIC = zeros(length(regmodel),1);
    allBIC = zeros(length(regmodel),1);
    allBriar = zeros(length(regmodel),1);
    allAUC = zeros(length(regmodel),1);
    allpihat = zeros(n,length(regmodel));
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
            allpihat(:,m) = pihat;
            allrownames{m} =  regmodel(m).name;
        end


%         alldev = round(alldev,1);
%         allAIC = round(allAIC,1);
%         allBIC = round(allBIC,1);
%         allBriar = round(allBriar,5);
%         allAUC = round(allAUC,2);

        disp(c_struct(c).group_names(2,:))
        Table = table(alld,alldev,allAIC,allBIC,allBriar,allAUC,'VariableNames',{'DOF','Dev','AIC','BIC','Briar','AUC'},'rowNames',allrownames);
        Table
        
        if e==1 && c==1
            savedalld = alld;
            savedalldev = alldev;
            savedallAIC = allAIC;
            savedallBIC = allBIC;
            savedallBriar = allBriar;
            savedallAUC = allAUC;
            savedrng = repmat(e,length(alld),1);
            savedc = repmat(c,length(alld),1);
        else
            savedalld = vertcat(savedalld,alld);
            savedalldev = vertcat(savedalldev,alldev);
            savedallAIC = vertcat(savedallAIC,allAIC);
            savedallBIC = vertcat(savedallBIC,allBIC);
            savedallBriar = vertcat(savedallBriar,allBriar);
            savedallAUC = vertcat(savedallAUC,allAUC);
            savedrng = vertcat(savedrng,repmat(e,length(alld),1));
            savedc = vertcat(savedc,repmat(c,length(alld),1));
        end
    end
    close all
end
save('rng_stats','savedalld','savedalldev','savedallAIC','savedallBIC','savedallBriar','savedallAUC','savedrng','savedc')

function [B,d,dev,AIC,BIC,Briar,AUC,pihat,X,Y,stats] = modelcomparison(modelparams,category,n)
   % Logistic Regression
%     [B,dev,stats] = mnrfit(modelparams,category);
%     pihat = mnrval(B,modelparams);
%     pihat = pihat(:,2);
%     RSS = nansum(stats.resid(:,2).^2);
    category(category==2) = 0;
    category = logical(category);
%     opts = statset('glmfit');
%     opts.MaxIter = 1000; % default value for glmfit is 100.
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