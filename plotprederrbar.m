function plotprederrbar(yhatfdobject, betastderrcell, argvals)
%  PLOTPREDERRBAR plots the predicted value along with the confidence
%  interval from individual fits

%  This has been modified from PLOTBETA by Amanda Zimmet

%  Arguments:
%  YHATFDOBJ      ... A cell object containing one or more functional
%                     parameter objects or functional data objects.
%  BETASTDERRCELL ... A cell object containing functional data objects
%                     for the standard error of the objects in
%                     BETAESTCELL.
% ------------------------------------------------------------------------

%  Check YHATFDOBJ

if strcmp(class(yhatfdobject), 'fdPar') || ...
   strcmp(class(yhatfdobject), 'fd')
    temp{1} = yhatfdobject;
    yhatfdobject = temp;
end

if ~iscell(yhatfdobject)
    error('YHATFDOBJECT is not a cell, fd, or fdpar object.');
end

%  Check BETASTDERRCELL

if nargin > 1

    if strcmp(class(betastderrcell), 'fd')
        temp{1} = betastderrcell;
        betastderrcell = temp;
    end

    if ~iscell(betastderrcell)
        error('BETASTDERRCELL is not a cell, or fd object.');
    end

end

%  Get range

if isa_fdPar(yhatfdobject{1})
    rangeval = getbasisrange(getbasis(getfd(yhatfdobject{1})));
elseif isa_fd(yhatfdobject{1})
    rangeval = getbasisrange(getbasis(yhatfdobject{1}));
else
    error(['A cell does not contain either a functional parameter ', ...
           'or a functional data object.']);
end

if nargin < 3
    argvals = linspace(rangeval(1),rangeval(2),51)';
end
n = length(argvals);
p = length(yhatfdobject);
figure()
for j=1:p
    if isa_fdPar(yhatfdobject{j})
        betavec = eval_fd(argvals, getfd(yhatfdobject{j}));
    elseif isa_fd(yhatfdobject{j})
        betavec = eval_fd(argvals, yhatfdobject{j});
    else
        error(['A cell of BETAESTCELL ', ...
               'does not contain either a functional parameter ', ...
               'or a functional data object.']);
    end
    plot(argvals, betavec, '-')
    line([rangeval(1),rangeval(2)],[0,0],'LineStyle',':','Color','r')
    if nargin > 1
        % Evaluate the standard error function at all time points
        betastderr = zeros(length(argvals),length(betastderrcell));
        for i=1:length(betastderrcell)
            betastderr(:,i) = eval_fd(argvals, betastderrcell{i});
        end
        
        % Keep only one copy of each group prediction
        betavec = unique(betavec','rows','stable')';
        
        % Compute the boundaries of the 95% confidence interval
        betavecp = betavec + 2.*betastderr;
        betavecm = betavec - 2.*betastderr;
        
        % Choose the color for each group
        Colors = {'r','g','b'};
        
        % Plot the shaded 95% Confidence Intervals
        for group = 1:size(betavec,2)-1
            h = fill([argvals;   argvals(n:-1:1,:)],[betavecp(:,group);betavecm(n:-1:1,group)],Colors{group}, 'LineStyle', ':');
            set(h,'facealpha',.1)
            hold on
        end
        
        % Add hash marks to one group
        % group = 1;
        % for i=1:n
        %     line([argvals(i),argvals(i)],[betavecm(i,group),betavecp(i,group)],'Color','red','LineWidth',0.25)
        % end
        
        % Plot the predicted lines
        for group=1:size(betavec,2)-1
            line(argvals, betavec(:,group),'Color',Colors{group},'LineWidth',4)
        end
        
        % Change x-axis to weeks of life
        if rangeval(1)==1 % only if you are starting at day of life == 1
            weekstartdays = rangeval(1):7:rangeval(2);
            weeks = (weekstartdays+6)/7;
            xticks(weekstartdays);
            xticklabels(string(weeks));
            xlim([rangeval(1), rangeval(2)])
            xlabel('Week of Life')
        end
        
    end
    title('\fontsize{16} Prediction with error bars')
    
    v = axis;
    axis([rangeval(1),rangeval(2),v(3),v(4)])
    if p > 1, pause; end
end
