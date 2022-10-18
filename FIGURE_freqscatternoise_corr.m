% FIGURE_freqscatternoise_corr

clear; close all

OBS_TableParams;
load(mattable);


flo = 1/1000;
fhi = 2.5;

TF_label_list = {'ZP-21'}; % minus tilt then compliance - will be a different figure for comparing the o

stations = OBS_table.Station;
networks = OBS_table.Network;
expfldr = OBS_table.ExperimentFolder;
experiments = OBS_table.Experiment;
expabrv = OBS_table.ExperimentAbbreviation;
waterdepth = OBS_table.WaterDepth;
sedthk = OBS_table.SedimentThickness;
statype = OBS_table.InstrumentDesign;
seismometer = OBS_table.Seismometer;
prestypes = OBS_table.PressureGauge;
deplystart = OBS_table.Start;
landdist = OBS_table.DistanceToLandCoarse;
starttimes = OBS_table.Start;
Zisgood = OBS_table.ZIsGood;
H1isgood = OBS_table.H1IsGood;
H2isgood = OBS_table.H2IsGood;
Pisgood = OBS_table.PIsGood;

[experiments_all,ia,ic] = unique(expfldr);
experimentnames_all = experiments(ia);
pressure_all = unique(prestypes);
statypes = unique(statype);
seistypes = unique(seismometer);

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
        
        
        for it =1:length(TF_label_list)
            isbad=0;
        if Zisgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            spc_stack_orig(:,ie,it) = NaN(size(fc));
            spec_diff(:,ie,it) = NaN(size(fc));
            isbad=1;
        end
        % for tilt containing
        if contains(TF_label,'1') == 1
            if H1isgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            spc_stack_orig(:,ie,it) = NaN(size(fc));
            spec_diff(:,ie,it) = NaN(size(fc));
            isbad=1;
            end
        end
        
        if contains(TF_label,'2') == 1
            if H2isgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            spc_stack_orig(:,ie,it) = NaN(size(fc));
            spec_diff(:,ie,it) = NaN(size(fc));
            isbad=1;
            end
        end
        
        if contains(TF_label,'P') == 1
            if Pisgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            spc_stack_orig(:,ie,it) = NaN(size(fc));
            spec_diff(:,ie,it) = NaN(size(fc));
            isbad=1;
            end
        end
        
        % bad conflicting component is excluded.
        if isbad==0
            spec = 10*log10(smoothSpectrum_octave(corr(idxx(it)).Zspec,f',fc_sta,ocav))+40*log10(wc);
            spc_stack(:,ie,it) = interp1(fc_sta,spec,fc);
        end
        end
        
        elevvec(ie) = waterdepth(ista);
        seiscol(ie,:) = scg(find(strcmp(statypes,statype(ista))==1),:);
        stavec(ie) = statype(ista);
        smtvec(ie) = seismometer(ista);
   
        ie = ie+1;
  
        
    end
end

nsta = length(elevvec);
%%
for ifreq = 1:length(flo_vec) % switch to scatter if possible...
    flo = flo_vec(ifreq);
    fhi = fhi_vec(ifreq);
    idxf = find(fc<=fhi & fc>=flo);
    fnew = fc(idxf);
    [LL,HH,FF] = noise_models(100);
    idxfm = find(FF<=fhi & FF>=flo);
    LAVG = mean(LL(idxfm));
    HAVG = mean(HH(idxfm));
    
    
    ii = 1;
    for iseis = 1:length(seistypes)
            idxseis = find(strcmp(seistypes(iseis),smtvec));
            clear meannoise watvec secol
            ie =1;
    for ista = 1:nsta
        if isempty(intersect(ista,idxseis))==1
            continue
        end
    noise = spc_stack(idxf,ista,ii);
    meannoise(ie) = mean(noise);
    watvec(ie)= elevvec(ista);
    secol(ie,:) = seiscol(ista,:);
    ie = ie+1;
    end
    
    figure(1)
    subplot(length(flo_vec),1,ifreq);
    scatter(watvec,meannoise,10,secol,char(seissym(iseis))); hold on
    end
    plot([0 6500], [LAVG LAVG],'-k','linewidth',1);
    plot([0 6500], [HAVG HAVG],'-k','linewidth',1);
    ylim([-190 -60])
    xlim([0 6500])
    set(gca,'FontSize',8)
    box on
xlabel('Water Depth (m)')
ylabel('Power (m^2/s^4/Hz, dB)');
    end
    

figure(1)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[8.5 11]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 3 6]);
filename=sprintf('%s/scatterZ_TPcorr_all.pdf',figoutpath);
print(gcf,'-dpdf',filename)
