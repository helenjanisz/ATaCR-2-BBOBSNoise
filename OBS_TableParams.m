% OBS_TableParams
% Script to keep naming and plotting organized

figoutpath = '~/FIGURES/OBS_Noise/';
if ~exist(figoutpath)
    mkdir(figoutpath);
end

exp_name = {'AACSE','LAU','PAPUA','ALBACORE','SEGMENT','MOANA','PLUME','ENAM',...
    'NOMELT','SCOOBA','GREECE','SHATSKY','MARIANA','HOBITSS','CHILE',...
    'BLANCO','TAIGER','CASCADIA_INITIATIVE','GORDA','NOOTKA','ALST','CASCADIA_KECK',...
    'GLIMPSE','MARIANA_O','YOUNG_ORCA'};
exp_abbr = {'Aa','La','Pa','Al','Se','Ma','Pl','En',...
    'Nm','Sc','Gr','Sh','Ms','Ho','Ch',...
    'Bl','Ta','Ci','Go','No','A','Ck',...
    'Gl','Mo','Yo'};

indir = '~/DATA_OBS_NOISE_4Paper/'; % locations where the output mat files from the ATaCR package live
specdir = 'AVG_STA_GOOD'; % station average spectra files 
TFdir = 'TRANSFUN_GOOD'; % station average transfer functions

extable = 'OBS_Working4Paper.xlsx'; %location of organizational table, excel format, akin to TableS2, this needs to be compiled by the user on their own
mattable = 'OBS_Working4Paper.mat'; %location of organizational table, matlab format

flo_vec = 1./[10 20 100 200]; % main frequency bands for average trends
fhi_vec = 1./[2 10 20 100];

% Colors for station types, consistent throughout
scg(1,:) = [166,206,227]; %SIO AB
scg(2,:) = [255,127,0]; % WHOI ARRA
scg(7,:) = [178,223,138]; % WHOI KECK
scg(4,:) = [202,178,214]; % LDEO APG
scg(5,:) = [106,61,154]; % LDEO DPG
scg(3,:) = [227,26,28]; % SIO BB
scg(8,:) = [31,120,180]; %TRM
scg(6,:) = [51,160,44]; % WHOI BB
scg= scg./255;

% symbols for seismometer types, consistent throughout
seissym{1} = 'o'; % guralp
seissym{2} = 'x'; % 240
seissym{3} = '^'; % compact
