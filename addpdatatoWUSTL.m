folder = 'V:/WUSTL/pulseOxStats'; 

filelist = dir([folder '/*.mat']);
n = length(filelist);

f = waitbar(0,'Percentage of Files Completed','Name','Progress');

load('X:/Data/wupdata.mat')
pdatafull = pdata;
pnamefull = pname;
clear pdata pname

tic
for i=1:n
    % Load data from patient file
    filename = filelist(i).name;
    fullfilename = fullfile(folder,filename);
    load(fullfilename)
    
    if ~exist('pname','var')
        currentid = str2double(filename(4:7));
        if sum(id==currentid)
            index = find(id==currentid);
            pdata = pdatafull(index,:);
            pname = pnamefull;
        end
    end
    
    fullfilename = fullfile([folder '/AmandaUpdated'],filename);
    if ~exist('pdata','var')
        disp(fullfilename)
        save(fullfilename,'x','xname','xt');
    else
        save(fullfilename,'pdata','pname','x','xname','xt');
    end
    clear pdata pname x xname xt 
    % Update the waitbar
    waitbar(i/n,f,[num2str(round(i/n*100,1)),'% Complete'])
end
        