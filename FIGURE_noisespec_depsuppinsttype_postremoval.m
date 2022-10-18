% FIGURE_noisespec_depsupp_posttiltremoval
% Figure to make supplementarty figures that show the noise spectra for
% each deployment individually

clear; close all

OBS_TableParams;

figoutpath = [figoutpath,'/Supplement/'];

flo = 1/1000;
fhi = 2.5;

load(mattable);
OBS_table_orig = OBS_table;

TF_label_list = {'ZP-21'}; % minus tilt then compliance - will be a different figure for comparing the o

stations = OBS_table.Station;
networks = OBS_table.Network;
expfldr = OBS_table.ExperimentFolder;
expname = OBS_table.Experiment;
expabrv = OBS_table.ExperimentAbbreviation;
waterdepth = OBS_table.WaterDepth;
Zisgood = OBS_table.ZIsGood;
H1isgood = OBS_table.H1IsGood;
H2isgood = OBS_table.H2IsGood;
Pisgood = OBS_table.PIsGood;
statype = OBS_table.InstrumentDesign;
seismometer = OBS_table.Seismometer;
prestypes = OBS_table.PressureGauge;

experiments_all = unique(expfldr);
pressure_all = unique(prestypes);
statypes = unique(statype);

ocav = 8;

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
    ie = 1;
    clear spc_stack ndays elev_vec isAPG
    exp = char(experiments_all(iexp));
    expn = char(experiments_all(iexp));
    inpath_dir = sprintf('%s%s/%s',indir,exp,TFdir);
    idx = strcmp(exp,expfldr);
    idxsta = find(idx==1);
    for id = 1:length(idxsta)
        ista = idxsta(id);
        expn = char(expname(ista));
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
            isbad=1;
        end
        % for tilt containing
        if contains(TF_label,'1') == 1
            if H1isgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            isbad=1;
            end
        end
        
        if contains(TF_label,'2') == 1
            if H2isgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            isbad=1;
            end
        end
        
        if contains(TF_label,'P') == 1
            if Pisgood(ista) ~=1
            spc_stack(:,ie,it) = NaN(size(fc));
            isbad=1;
            end
        end
        
        % bad conflicting component is excluded.
        if isbad==0
            spec = 10*log10(smoothSpectrum_octave(corr(idxx(it)).Zspec,f',fc_sta,ocav))+40*log10(wc);
            spc_stack(:,ie,it) = interp1(fc_sta,spec,fc);
        end
        end
        
        elev_vec(ie) = waterdepth(ista);
        stavec(ie) = statype(ista);
        ie = ie+1;
    end

     fnew = fc;
    tt = 0;
    nsta = length(elev_vec);
    
    for isttyp = 1:length(statypes)
        idxsttyp = find(strcmp(statypes(isttyp),stavec));
        if isempty(idxsttyp)==1
            continue
        end
        
    % corr spectra
        ii = 1;
        suff = expn;
        noise = spc_stack(:,:,ii);
        figure(1);
        subplot(6,3,iexp)
        for ista = 1:nsta
            if isempty(find(ista==idxsttyp))==1
                continue
            end
              semilogx(fnew,noise(:,ista),'-','Color',scg(isttyp,:),'linewidth',.5);hold on
        end
[LL,HH,FF] = noise_models(100);
        semilogx(FF,HH,'-k', 'linewidth',2); hold on;
        semilogx(FF,LL,'-k', 'linewidth',2); hold on;
        xlim([flo,1]);
        ylim([-200 -50])
        
        text(1/1000,-40,suff,'FontSize',8)
        if iexp~= 1 & iexp ~= 4 & iexp ~= 7 & iexp ~= 10 & iexp ~= 13 & iexp ~= 16
            set(gca,'ytick',[])
        end
        if iexp<16
            set(gca,'xtick',[])
        end
        if iexp>=16
            xlabel('Frequency (Hz)')
        end
        if iexp== 1 | iexp == 4 | iexp == 7 | iexp == 10 | iexp == 13 | iexp == 16
            ylabel('m^2/s^4/Hz, dB')
        end
        xticks([1/1000,1/100,1/10,1])
        set(gca,'FontSize',8)



    end
    
end
figure(1)

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 11]);
filename=sprintf('%s/NoiseSpec_StaType_DepZ_tiltcomprm.pdf',figoutpath);
print(gcf,'-dpdf',filename)

