load('X:\Amanda\FuncDataAnalysis\UVA_CU_WUSTL_FPCA_Stats_day.mat')
GoodIDs = zeros(length(pid),1);
goodidids = pid;
goodidinst = inst;

% UVA
% load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\KarenIDs.mat')
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_ega28_UVA.mat') % loads goodids, which contains binary representation of whether or not there are 7 days of data in first 28 days
% GoodIDsUVA = KarenIDs(goodids);
% [A,~] = ismember(pid(inst==1),GoodIDsUVA);
[A,~] = ismember(pid(inst==1),goodpids(goodids));
GoodIDs(inst==1) = A;

% CU
% load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\CUAim3Eligible.mat')
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_ega28_CU.mat')
% GoodIDsCU = CUAim3Eligible.DeidentifiedMRN(CUAim3Eligible.Aim3Inclusion==0); % 0 = to be included
% [A,~] = ismember(GoodIDsCU,pid); % Only keep the list of GoodIDs which have data files
% GoodIDsCU = GoodIDsCU(A);
% GoodIDsCU = GoodIDsCU(goodids);
% [A,~] = ismember(pid(inst==2),GoodIDsCU);
[A,~] = ismember(pid(inst==2),goodpids(goodids));
GoodIDs(inst==2) = A;

% WUSTL
load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\SevenIn28_ega28_WUSTL.mat')
% load('X:\Amanda\FuncDataAnalysis\TablesAndSpreadsheets\WUSTLAim3Demographics181218.mat')
% GoodIDsWUSTL = pid(inst==3); % ids from WUSTL with data
% GoodIDsWUSTL = GoodIDsWUSTL(goodids); % ids which have at least 7/first 28 days and are eligible
% [A,~] = ismember(pid(inst==3),GoodIDsWUSTL); % figure out which indices of the pid list meet all the criteria
[A,~] = ismember(pid(inst==3),goodpids(goodids)); % figure out which indices of the pid list meet all the criteria
GoodIDs(inst==3) = A; % save which babies at WUSTL meet all criteria in the GoodIDs list

GoodIDs = logical(GoodIDs); % convert 1's and 0's to a logical vector
save('GoodIDs_ega28','GoodIDs','goodidids','goodidinst')