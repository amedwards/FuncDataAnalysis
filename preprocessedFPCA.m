function preprocessedFPCA(institution,hrdaywk)
% preprocessedFPCA reads in stats files and converts them into the format
% required to run FPCA via preprocessedFPCA_run.m. 
% The user needs to select which institution from which to save data by 
% changing the institution value. If you would like to change which 
% variables are pulled from the file, change vname.


% Choose which dataset you would like to read in:
% 1 = UVA, 2 = CU, 3 = WUSTL
% institution = 
% Choose whether you would like to lump by
% 1 = HR, 2 = Day, 3 = Week
% hrdaywk = 

switch institution
    case 1
        % UVA Data
        folder = 'V:/NICU/pulseOxStats';
        savefile = 'UVA_FPCA_Stats';
        vname={'Mean HR','Mean SPO2-R','Mean SPO2-%','SD HR','SD SPO2-R','SD SPO2-%','Skewness HR','Skewness SPO2-R','Skewness SPO2-%','Kurtosis HR','Kurtosis SPO2-R','Kurtosis SPO2-%',...
            'Max XC HR SPO2-%','Lag Max XC HR SPO2-%','Min XC HR SPO2-%','Lag Min XC HR SPO2-%','Max XC SPO2-R SPO2-%','Lag Max XC SPO2-R SPO2-%','Min XC SPO2-R SPO2-%','Lag Min XC SPO2-R SPO2-%'}';
    case 2
        % CU Data
%         folder = 'V:/CU/pulseOxStats';
        folder = 'X:/Amanda/Columbia/ColumbiaThirdBatchJob/forAmanda_01092019/pulseOxStats';
        savefile = 'CU_FPCA_Stats';
        vname = {'Mean vt_hr','Mean vt_spo2r','Mean vt_spo2','SD vt_hr','SD vt_spo2r','SD vt_spo2','Skewness vt_hr','Skewness vt_spo2r','Skewness vt_spo2','Kurtosis vt_hr','Kurtosis vt_spo2r','Kurtosis vt_spo2',...
            'Max XC vt_hr vt_spo2','Lag Max XC vt_hr vt_spo2','Min XC vt_hr vt_spo2','Lag Min XC vt_hr vt_spo2', 'Max XC vt_spo2r vt_spo2','Lag Max XC vt_spo2r vt_spo2','Min XC vt_spo2r vt_spo2','Lag Min XC vt_spo2r vt_spo2'}';
    case 3
        % WUSTL Data
        folder = 'V:/WUSTL/pulseOxStats'; 
        savefile = 'WUSTL_FPCA_Stats';
        vname={'Mean HR','Mean SPO2-R','Mean SPO2-%','SD HR','SD SPO2-R','SD SPO2-%','Skewness HR','Skewness SPO2-R','Skewness SPO2-%','Kurtosis HR','Kurtosis SPO2-R','Kurtosis SPO2-%',...
            'Max XC HR SPO2-%','Lag Max XC HR SPO2-%','Min XC HR SPO2-%','Lag Min XC HR SPO2-%','Max XC SPO2-R SPO2-%','Lag Max XC SPO2-R SPO2-%','Min XC SPO2-R SPO2-%','Lag Min XC SPO2-R SPO2-%'}';
end



filelist = dir([folder '/*.mat']);
n = length(filelist);

% Choose which variables you would like to store
% vname={'Mean HR','Mean SPO2-%'}';

nv=length(vname);
vcols = zeros(nv,1);

% Initialize variables
dh=60*60; % number of seconds in an hour
dy=86400; % number of seconds in a day
maxthr = 12*7*24; % maximum time since birth (hrs) (12 weeks chronologic age)
maxtday = 12*7; % maximum time since birth (days) (12 weeks chronologic age)
maxtweek = 12; % maximum time since birth (weeks) (12 weeks chronologic age)

% Name the save file and initialize empty matrices
switch hrdaywk
    case 1
        savefile = [savefile '_hr'];
        datahr = ones(maxthr,nv,n)*nan;
        numsamps = zeros(maxthr,nv,n);
    case 2
        savefile = [savefile '_day'];
        dataday = ones(maxtday,nv,n)*nan;
        numsamps = zeros(maxtday,nv,n);
    case 3
        savefile = [savefile '_wk'];
        dataweek = ones(maxtweek,nv,n)*nan;
        numsamps = zeros(maxtweek,nv,n);
end

pid = zeros(n,1);
pbd = zeros(n,1);
pbw = zeros(n,1);
pega = zeros(n,1);
pgen = zeros(n,1);

f = waitbar(0,'Percentage of Files Completed','Name','Progress');

tic
for i=1:n
    % Load data from patient file
    filename = filelist(i).name;
    fullfilename = fullfile(folder,filename);
    load(fullfilename)
    
    if exist('pname','var')
        
        % Find patient ID
        idindex = find(contains(pname,'PatientID'));
        
        if institution == 2
            id = str2double(extractBetween(filename,'_','_'));
            pdatarow = find(pdata(:,idindex)==id);
        end

        % Find birth date
        bdayindex = find(contains(pname,'BirthDate'));
        
        % Find birth weight
        bwindex = find(contains(pname,'BirthWeight'));
        
        % Find Estimated Gestational Age
        egaindex = find(contains(pname,'GestAge'));
        if ~isempty(egaindex) % Need to handle possibility of 'GestAgeDays'
            egaindex = ismember(pname,'GestAge');
        end
        
        % Find Gender
        genindex = find(contains(pname,'Sex'));
        
        if institution==2
            pid(i) = pdata(pdatarow,idindex);
            pbd(i) = pdata(pdatarow,bdayindex);
            pbw(i) = pdata(pdatarow,bwindex);
            pega(i) = pdata(pdatarow,egaindex);
            pgen(i) = pdata(pdatarow,genindex);
        else
            pid(i) = pdata(idindex);
            pbd(i) = pdata(bdayindex); % seconds of age since some arbitrary time point
            pbw(i) = pdata(bwindex);
            pega(i) = pdata(egaindex);
            pgen(i) = pdata(genindex);
        end
    elseif institution==3
        pid(i) = str2double(filename(4:7));
    end
    
    % Find chronologic age time stamps
    agesec = xt; % chronologic age in seconds
    agehours = ceil(agesec/60/60); % chronologic age in hours
    agedays = ceil(agesec/60/60/24); % chronologic age in days
    ageweeks = ceil(agesec/60/60/24/7); %chronologic age in weeks

    % Find the data columns that correspond to the variables of interest
    for v=1:nv
        vcols(v) = find(ismember(xname,vname(v)));
    end
    
    switch hrdaywk
        case 1
            % Take hourly data average
            for h=1:maxthr % iterate through hours since birth
                indices=(agehours==h); % all samples taken in that hour
                numsamps(h,:,i) = sum(~isnan(x(indices,vcols)));
                datahr(h,:,i)=nanmean(x(indices,vcols),1); % take the average over the hour of data
                matrixtosave = 'datahr';
            end
        case 2
            % Take daily data average
            for d=1:maxtday % iterate through days since birth
                indices=(agedays==d); % all samples taken in that day
                numsamps(d,:,i) = sum(~isnan(x(indices,vcols)));
                dataday(d,:,i)=nanmean(x(indices,vcols),1); % take the average over the day of data
                matrixtosave = 'dataday';
            end
        case 3
            % Take weekly data average
            for w=1:maxtweek % iterate through days since birth
                indices=(ageweeks==w); % all samples taken in that day
                numsamps(w,:,i) = sum(~isnan(x(indices,vcols)));
                dataweek(w,:,i) = nanmean(x(indices,vcols),1); % take the average over the week of data
                matrixtosave = 'dataweek';
            end
    end
    
    % Update the waitbar
    waitbar(i/n,f,[num2str(round(i/n*100,1)),'% Complete'])
    clear pname
end
inst = institution*ones(n,1);
close(f)
save(savefile,'n','nv',matrixtosave,'vname','pgen','pid','pbd','pbw','pega','inst','numsamps')