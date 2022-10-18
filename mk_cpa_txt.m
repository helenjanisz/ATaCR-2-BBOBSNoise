% mk_cpa_txt
% quick script to make the coherence, phase, and admittance files
% not octave averaged

clear;
close all

OBS_TableParams;

outpath1 = './Coh_Supp_Files';
if ~exist(outpath1,'dir')
    mkdir(outpath1);
end

outpath2 = './Phs_Supp_Files';
if ~exist(outpath2,'dir')
    mkdir(outpath2);
end

outpath3 = './Adm_Supp_Files';
if ~exist(outpath3,'dir')
    mkdir(outpath3);
end

flo = 1/1000;
fhi = 2.5;

load(mattable);
OBS_table_orig = OBS_table;

stations = OBS_table.Station;
networks = OBS_table.Network;
expfldr = OBS_table.ExperimentFolder;
expabrv = OBS_table.ExperimentAbbreviation;
waterdepth = OBS_table.WaterDepth;
% isgood = OBS_table.IsGood;
Zisgood = OBS_table.ZIsGood;
H1isgood = OBS_table.H1IsGood;
H2isgood = OBS_table.H2IsGood;
Pisgood = OBS_table.PIsGood;
statype = OBS_table.InstrumentDesign;
seismometer = OBS_table.Seismometer;
prestypes = OBS_table.PressureGauge;

experiments_all = unique(expfldr);
pressure_all = unique(prestypes);

ocav = 16;

ie = 1;

ii = 1;
Ts = 1/fhi;
while Ts<1/flo
    Tl = Ts*2;
    Tc(ii) = sqrt(Ts*Tl);
    Ts=Ts*2^(1/ocav);
    ii = ii+1;
end

fc = 1./Tc;
for iexp = 1:length(experiments_all)
    %     for iexp = 14
    exp = char(experiments_all(iexp));
    inpath_dir = sprintf('%s%s/%s',indir,exp,specdir);
    idx = strcmp(exp,expfldr);
    idxsta = find(idx==1);
    for id = 1:length(idxsta)
        ista = idxsta(id);
        sta = char(stations(ista));
        net = char(networks(ista));
        station_file = sprintf('%s/%s%s_spectraavg.mat',inpath_dir,net,sta);
        if exist(station_file,'file') == 0
            continue
        end
        load(station_file);
        
        f = avgprop.params.f;
        clear Tcsta
        ii = 1;
        Tssta = 1/max(f);
        while Tssta<1/f(2)
            Tl = Tssta*2;
            Tcsta(ii) = sqrt(Tssta*Tl);
            Tssta=Tssta*2^(1/ocav);
            ii = ii+1;
        end
        fc_sta =1./Tcsta;
        PI = 4*atan(1);
        w = 2*PI*avgprop.params.f;
        wc = 2*PI*fc_sta;
        
        % coherences
        coh1z = abs(staavg.cross.c1z_mean).^2./(staavg.power.c11_mean.*staavg.power.czz_mean);
        coh2z = abs(staavg.cross.c2z_mean).^2./(staavg.power.c22_mean.*staavg.power.czz_mean);
        cohpz = abs(staavg.cross.cpz_mean).^2./(staavg.power.cpp_mean.*staavg.power.czz_mean);
        coh1p = abs(staavg.cross.c1p_mean).^2./(staavg.power.c11_mean.*staavg.power.cpp_mean);
        coh2p = abs(staavg.cross.c2p_mean).^2./(staavg.power.c22_mean.*staavg.power.cpp_mean);
        coh12 = abs(staavg.cross.c12_mean).^2./(staavg.power.c22_mean.*staavg.power.c11_mean);
        
        %phases
        ph1z = 180/pi.*atan2(imag(staavg.cross.c1z_mean),real(staavg.cross.c1z_mean));
        ph2z = 180/pi.*atan2(imag(staavg.cross.c2z_mean),real(staavg.cross.c2z_mean));
        phpz = 180/pi.*atan2(imag(staavg.cross.cpz_mean),real(staavg.cross.cpz_mean));
        ph1p = 180/pi.*atan2(imag(staavg.cross.c1p_mean),real(staavg.cross.c1p_mean));
        ph2p = 180/pi.*atan2(imag(staavg.cross.c2p_mean),real(staavg.cross.c2p_mean));
        ph12 = 180/pi.*atan2(imag(staavg.cross.c12_mean),real(staavg.cross.c12_mean));
        
        % admittances
        ad1z = abs(staavg.cross.c1z_mean)./staavg.power.c11_mean;
        ad2z = abs(staavg.cross.c2z_mean)./staavg.power.c22_mean;
        adpz = abs(staavg.cross.cpz_mean)./staavg.power.cpp_mean;
        ad1p = abs(staavg.cross.c1p_mean)./staavg.power.cpp_mean;
        ad2p = abs(staavg.cross.c2p_mean)./staavg.power.cpp_mean;
        ad12 = abs(staavg.cross.c12_mean)./staavg.power.c11_mean;
        
        
        % coherence
        filename1 = sprintf('%s/%s_%s_coh_Z1.txt',outpath1,net,sta);
        filename2 = sprintf('%s/%s_%s_coh_Z2.txt',outpath1,net,sta);
        filename3 = sprintf('%s/%s_%s_coh_ZP.txt',outpath1,net,sta);
        filename4 = sprintf('%s/%s_%s_coh_1P.txt',outpath1,net,sta);
        filename5 = sprintf('%s/%s_%s_coh_2P.txt',outpath1,net,sta);
        filename6 = sprintf('%s/%s_%s_coh_12.txt',outpath1,net,sta);
        
        Z1mat = [f;coh1z];
        Z2mat = [f;coh2z];
        ZPmat = [f;cohpz];
        H1Pmat = [f;coh1p];
        H2Pmat = [f;coh2p];
        H12mat = [f;coh12];
       
        % write files
        writematrix(Z1mat',filename1,'Delimiter','tab');
        writematrix(Z2mat',filename2,'Delimiter','tab');
        writematrix(ZPmat',filename3,'Delimiter','tab');
        writematrix(H1Pmat',filename4,'Delimiter','tab');
        writematrix(H2Pmat',filename5,'Delimiter','tab');
        writematrix(H12mat',filename6,'Delimiter','tab');
        
        %phases
        filename1 = sprintf('%s/%s_%s_phs_Z1.txt',outpath2,net,sta);
        filename2 = sprintf('%s/%s_%s_phs_Z2.txt',outpath2,net,sta);
        filename3 = sprintf('%s/%s_%s_phs_ZP.txt',outpath2,net,sta);
        filename4 = sprintf('%s/%s_%s_phs_1P.txt',outpath2,net,sta);
        filename5 = sprintf('%s/%s_%s_phs_2P.txt',outpath2,net,sta);
        filename6 = sprintf('%s/%s_%s_phs_12.txt',outpath2,net,sta);
        
        Z1mat = [f;ph1z];
        Z2mat = [f;ph2z];
        ZPmat = [f;phpz];
        H1Pmat = [f;ph1p];
        H2Pmat = [f;ph2p];
        H12mat = [f;ph12];
       
        % write files
        writematrix(Z1mat',filename1,'Delimiter','tab');
        writematrix(Z2mat',filename2,'Delimiter','tab');
        writematrix(ZPmat',filename3,'Delimiter','tab');
        writematrix(H1Pmat',filename4,'Delimiter','tab');
        writematrix(H2Pmat',filename5,'Delimiter','tab');
        writematrix(H12mat',filename6,'Delimiter','tab');
        
        % admittances
        
        filename1 = sprintf('%s/%s_%s_adm_Z1.txt',outpath3,net,sta);
        filename2 = sprintf('%s/%s_%s_adm_Z2.txt',outpath3,net,sta);
        filename3 = sprintf('%s/%s_%s_adm_ZP.txt',outpath3,net,sta);
        filename4 = sprintf('%s/%s_%s_adm_1P.txt',outpath3,net,sta);
        filename5 = sprintf('%s/%s_%s_adm_2P.txt',outpath3,net,sta);
        filename6 = sprintf('%s/%s_%s_adm_12.txt',outpath3,net,sta);
        
        Z1mat = [f;ad1z];
        Z2mat = [f;ad2z];
        ZPmat = [f;adpz];
        H1Pmat = [f;ad1p];
        H2Pmat = [f;ad2p];
        H12mat = [f;ad12];
       
        % write files
        writematrix(Z1mat',filename1,'Delimiter','tab');
        writematrix(Z2mat',filename2,'Delimiter','tab');
        writematrix(ZPmat',filename3,'Delimiter','tab');
        writematrix(H1Pmat',filename4,'Delimiter','tab');
        writematrix(H2Pmat',filename5,'Delimiter','tab');
        writematrix(H12mat',filename6,'Delimiter','tab');
        
    end
    
end

