function mergedatasets(dataset1,dataset2,mergename)
% mergedatasets('UVA_FPCA_Stats.mat','WUSTL_FPCA_Stats.mat','UVA_WUSTL_FPCA_Stats')
% mergedatasets('UVA_FPCA_Stats_wk.mat','CU_FPCA_Stats_wk.mat','UVA_CU_FPCA_Stats_wk.mat')
% mergedatasets('UVA_CU_FPCA_Stats_wk.mat','WUSTL_FPCA_Stats_wk.mat','UVA_CU_WUSTL_FPCA_Stats_wk.mat')
% mergedatasets('UVA_FPCA_Stats_day.mat','CU_FPCA_Stats_day.mat','UVA_CU_FPCA_Stats_day.mat')
% mergedatasets('UVA_CU_FPCA_Stats_day.mat','WUSTL_FPCA_Stats_day.mat','UVA_CU_WUSTL_FPCA_Stats_day.mat')

% load(dataset1,'dataweek','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
load(dataset1,'dataday','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
% dataweekmerge = dataweek;
datadaymerge = dataday;
nmerge = n;
nvmerge = nv;
pbdmerge = pbd;
pbwmerge = pbw;
pegamerge = pega;
pgenmerge = pgen;
pidmerge = pid;
institution1 = inst;
numsampsmerge = numsamps;

% load(dataset2,'dataweek','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
load(dataset2,'dataday','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
% dataweek = cat(3,dataweekmerge,dataweek);
dataday = cat(3,datadaymerge,dataday);
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

% save(mergename,'dataweek','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')
save(mergename,'dataday','n','nv','pbd','pbw','pega','pgen','pid','vname','inst','numsamps')