%% Create FPCA Stats from pulseOxStats
disp('Processing UVA data...')
preprocessedFPCA(1,1) % create UVA_FPCA_Stats_hr
preprocessedFPCA(1,2) % create UVA_FPCA_Stats_day
preprocessedFPCA(1,3) % create UVA_FPCA_Stats_wk

disp('Processing CU data...')
preprocessedFPCA(2,1) % create CU_FPCA_Stats_hr
preprocessedFPCA(2,2) % create CU_FPCA_Stats_day
preprocessedFPCA(2,3) % create CU_FPCA_Stats_wk

disp('Processing WUSTL data...')
preprocessedFPCA(3,1) % create WUSTL_FPCA_Stats_hr
preprocessedFPCA(3,2) % create WUSTL_FPCA_Stats_day
preprocessedFPCA(3,3) % create WUSTL_FPCA_Stats_wk

%% Merge Datasets
disp('Merging hourly FPCA stats files...')
mergedatasets('UVA_FPCA_Stats_hr.mat','CU_FPCA_Stats_hr.mat','UVA_CU_FPCA_Stats_hr.mat')
mergedatasets('UVA_CU_FPCA_Stats_hr.mat','WUSTL_FPCA_Stats_hr.mat','UVA_CU_WUSTL_FPCA_Stats_hr.mat')

disp('Merging daily FPCA stats files...')
mergedatasets('UVA_FPCA_Stats_day.mat','CU_FPCA_Stats_day.mat','UVA_CU_FPCA_Stats_day.mat')
mergedatasets('UVA_CU_FPCA_Stats_day.mat','WUSTL_FPCA_Stats_day.mat','UVA_CU_WUSTL_FPCA_Stats_day.mat')

disp('Merging weekly FPCA stats files...')
mergedatasets('UVA_FPCA_Stats_wk.mat','CU_FPCA_Stats_wk.mat','UVA_CU_FPCA_Stats_wk.mat')
mergedatasets('UVA_CU_FPCA_Stats_wk.mat','WUSTL_FPCA_Stats_wk.mat','UVA_CU_WUSTL_FPCA_Stats_wk.mat')

%% Find babies with 7 out of first 28 days with data
addpath('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets');
disp('Finding babies with 7/28 days of data...')
AmountofDataAvailable; % Creates SevenIn28_UVA, SevenIn28_CU, SevenIn28_WUSTL

%% Create the GoodIDs variable
disp('Getting GoodIDs')
getGoodIDs