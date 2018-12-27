function mergedatasets(dataset1,dataset2,mergename)

% mergedatasets('UVA_FPCA_Stats_wk.mat','CU_FPCA_Stats_wk.mat','UVA_CU_FPCA_Stats_wk.mat')
% mergedatasets('UVA_CU_FPCA_Stats_wk.mat','WUSTL_FPCA_Stats_wk.mat','UVA_CU_WUSTL_FPCA_Stats_wk.mat')

% mergedatasets('UVA_FPCA_Stats_day.mat','CU_FPCA_Stats_day.mat','UVA_CU_FPCA_Stats_day.mat')
% mergedatasets('UVA_CU_FPCA_Stats_day.mat','WUSTL_FPCA_Stats_day.mat','UVA_CU_WUSTL_FPCA_Stats_day.mat')

% mergedatasets('UVA_FPCA_Stats_hr.mat','CU_FPCA_Stats_hr.mat','UVA_CU_FPCA_Stats_hr.mat')
% mergedatasets('UVA_CU_FPCA_Stats_hr.mat','WUSTL_FPCA_Stats_hr.mat','UVA_CU_WUSTL_FPCA_Stats_hr.mat')

load(dataset1,'data*','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
if exist('datahr','var')
    datamerge1 = datahr;
elseif exist('dataday','var')
    datamerge1 = dataday;
elseif exist('dataweek','var')
    datamerge1 = dataweek;
else
    datamerge1 = datamerge;
end

nmerge = n;
nvmerge = nv;
pbdmerge = pbd;
pbwmerge = pbw;
pegamerge = pega;
pgenmerge = pgen;
pidmerge = pid;
institution1 = inst;
numsampsmerge = numsamps;

load(dataset2,'data*','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
if exist('datahr','var')
    datamerge2 = datahr;
elseif exist('dataday','var')
    datamerge2 = dataday;
elseif exist('dataweek','var')
    datamerge2 = dataweek;
else
    datamerge2 = datamerge;
end

datamerge = cat(3,datamerge1,datamerge2);
n = nmerge + n;
if nv~=nvmerge
    error('The number of variables is not the same between datasets.')
end
pbd = [pbdmerge;pbd];
pbw = [pbwmerge;pbw];
pega = [pegamerge;pega];
pgen = [pgenmerge;pgen];
pid = [pidmerge;pid];
inst = [institution1; inst];
numsamps = cat(3,numsampsmerge,numsamps);

save(mergename,'datamerge','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')