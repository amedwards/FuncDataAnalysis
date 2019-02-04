load('X:\Amanda\FuncDataAnalysis\UVA_CU_WUSTL_FPCA_Stats_day.mat')
GoodIDs = zeros(length(pid),1);
goodidids = pid;
goodidinst = inst;

% UVA
% load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_ega28_UVA.mat') % loads goodids, which contains binary representation of whether or not there are 7 days of data in first 28 days
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_UVA.mat') % loads goodids, which contains binary representation of whether or not there are 7 days of data in first 28 days
[A,~] = ismember(pid(inst==1),goodpids(goodids));
GoodIDs(inst==1) = A;

% CU
% load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_ega28_CU.mat')
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_CU.mat')
[A,~] = ismember(pid(inst==2),goodpids(goodids));
GoodIDs(inst==2) = A;

% WUSTL
% load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_ega28_WUSTL.mat')
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_WUSTL.mat')
[A,~] = ismember(pid(inst==3),goodpids(goodids)); % figure out which indices of the pid list meet all the criteria
GoodIDs(inst==3) = A; % save which babies at WUSTL meet all criteria in the GoodIDs list

GoodIDs = logical(GoodIDs); % convert 1's and 0's to a logical vector
save('GoodIDs','GoodIDs','goodidids','goodidinst')