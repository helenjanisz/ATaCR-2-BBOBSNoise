% mk_spec_corr_txt
% quick script to make the text files for the tilt and compliance corrected
% spectra
% octave averaged

clear; close all

OBS_TableParams;

outpath = './Corr_Spectra_Supp_Files';
if ~exist(outpath,'dir')
    mkdir(outpath);
end

flo = 1/1000;
fhi = 2.5;

load(mattable);
OBS_table_orig = OBS_table;

%DO NOT EDIT ORDER
TF_label_list = {'Z2-1','ZP','ZP-21'};


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

ocav = 8;

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
    exp = char(experiments_all(iexp));
    inpath_dir = sprintf('%s%s/%s',indir,exp,TFdir);
    idx = strcmp(exp,expfldr);
    idxsta = find(idx==1);
    for id = 1:length(idxsta)
        ista = idxsta(id);
        sta = char(stations(ista));
        net = char(networks(ista));
        station_file = sprintf('%s/%s%s/%s%s_AVERAGE_transfun.mat',inpath_dir,net,sta,net,sta);
        if exist(station_file,'file') == 0
            continue
        end
        load(station_file);
        
        f = transprop.params.f;
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
        w = 2*PI*transprop.params.f;
        wc = 2*PI*fc_sta;
        
        for it = 1:length(TF_label_list)
            
            TF_label = char(TF_label_list(it));
            for ic = 1:length(corr)
                if strcmp(corr(ic).label,TF_label)==1;
                    idxx(it) = ic;
                    break
                end
            end
        end
        
        
        % for compliance containing - need to edit these such that only the
        % bad conflicting component is excluded.
        spec = 10*log10(smoothSpectrum_octave(corr(idxx(it)).Zspec,f',fc_sta,ocav))+40*log10(wc);
        spc_stack(:,it) = interp1(fc_sta,spec,fc);
        
        filenamez = sprintf('%s/%s_%s_spectra_Z_tiltcorrected.txt',outpath,net,sta);
        filename1 = sprintf('%s/%s_%s_spectra_Z_compcorrected.txt',outpath,net,sta);
        filename2 = sprintf('%s/%s_%s_spectra_Z_allcorrected.txt',outpath,net,sta);
        
        
        Zmat = [fc',spc_stack(:,1)];
        H1mat = [fc',spc_stack(:,2)];
        H2mat = [fc',spc_stack(:,3)];
        
        % write files
        writematrix(Zmat,filenamez,'Delimiter','tab');
        writematrix(H1mat,filename1,'Delimiter','tab');
        writematrix(H2mat,filename2,'Delimiter','tab');
        
        
    end
 
end


