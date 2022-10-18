% FIGURE_noisespec_depsupp
% Figure to make supplementarty figures that show the noise spectra for
% each deployment individually
% colored by instrument type

clear; close all

OBS_TableParams;

figoutpath = [figoutpath,'/Supplement/'];

flo = 1/1000;
fhi = 2.5;

load(mattable);
OBS_table_orig = OBS_table;

stations = OBS_table.Station;
networks = OBS_table.Network;
expfldr = OBS_table.ExperimentFolder;
expname = OBS_table.Experiment;
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
    inpath_dir = sprintf('%s%s/%s',indir,exp,specdir);
    idx = strcmp(exp,expfldr);
    idxsta = find(idx==1);
    for id = 1:length(idxsta)
        ista = idxsta(id);
        expn = char(expname(ista));
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
        
        % noise, in db
        if Zisgood(ista) ==1
            spec = 10*log10(smoothSpectrum_octave(staavg.power.czz_mean,avgprop.params.f,fc_sta,ocav))+40*log10(wc);
            spc_stack(:,ie,1) = interp1(fc_sta,spec,fc);
        else
            spc_stack(:,ie,1) = NaN(size(fc));
        end
        if H1isgood(ista) ==1
            spec = 10*log10(smoothSpectrum_octave(staavg.power.c11_mean,avgprop.params.f,fc_sta,ocav))+40*log10(wc);
            spc_stack(:,ie,2) = interp1(fc_sta,spec,fc);
        else
            spc_stack(:,ie,2) = NaN(size(fc));
        end
        if H2isgood(ista) ==1
            spec = 10*log10(smoothSpectrum_octave(staavg.power.c22_mean,avgprop.params.f,fc_sta,ocav))+40*log10(wc);
            spc_stack(:,ie,3) = interp1(fc_sta,spec,fc);
        else
            spc_stack(:,ie,3) = NaN(size(fc));
        end
        if Pisgood(ista) ==1
            spec = 10*log10(smoothSpectrum_octave(staavg.power.cpp_mean,avgprop.params.f,fc_sta,ocav));
            spc_stack(:,ie,4) = interp1(fc_sta,spec,fc);
        else
            spc_stack(:,ie,4) = NaN(size(fc));
        end
        
        elev_vec(ie) = waterdepth(ista);
        stavec(ie) = statype(ista);
        % APG v DPG Check
        if strcmp(prestypes(ista),'APG')==1
            isAPG(ie) = 1;
        else
            isAPG(ie) = 0;
        end
        ndays(ie) = avgprop.params.ndays;
        ie = ie+1;
    end
    
    fnew = fc;
    nsta = length(elev_vec);
    
    for isttyp = 1:length(statypes)
        idxsttyp = find(strcmp(statypes(isttyp),stavec));
        if isempty(idxsttyp)==1
            continue
        end
        
        % vertical spectra
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
%         if isttyp == 1
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
%         end

        
        % horizontal spectra
        figure(2)
        ii = 2;
        noise = spc_stack(:,:,ii);
        
        subplot(6,3,iexp)
        for ista = 1:nsta
            if isempty(find(ista==idxsttyp))==1
                continue
            end
            semilogx(fnew,noise(:,ista),'-','Color',scg(isttyp,:),'linewidth',.5);   hold on
        end
        ii = 3;
        noise = spc_stack(:,:,ii);
        for ista = 1:nsta
            if isempty(find(ista==idxsttyp))==1
                continue
            end
            semilogx(fnew,noise(:,ista),'-','Color',scg(isttyp,:),'linewidth',.5);        hold on
        end
%         if isttyp == 1
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
%         end
        
        
        % APG
        figure(3)
        ii = 4;
        noise = spc_stack(:,:,ii);
        
        subplot(6,3,iexp)
        for ista = 1:nsta
            noiseplot = noise(:,ista);
            if isempty(find(ista==idxsttyp))==1
            noiseplot = NaN(size(noise(:,ista)));
            end
            if isAPG(ista)==0
            noiseplot = NaN(size(noise(:,ista)));
            end
            semilogx(fnew,noiseplot,'-','Color',scg(isttyp,:),'linewidth',.5);   hold on
        end
        xlim([flo,1]);
        ylim([-150 130])
        text(1/1000,150,suff,'FontSize',8)
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
            ylabel('Pa^2/Hz, dB')
        end
        xticks([1/1000,1/100,1/10,1])
        set(gca,'FontSize',8)
        
        % DPG
        figure(4)
        ii = 4;
        noise = spc_stack(:,:,ii);
        
        subplot(6,3,iexp)
        for ista = 1:nsta
            noiseplot = noise(:,ista);
            if isempty(find(ista==idxsttyp))==1
            noiseplot = NaN(size(noise(:,ista)));
            end
            if isAPG(ista)==1
            noiseplot = NaN(size(noise(:,ista)));
            end
            semilogx(fnew,noiseplot,'-','Color',scg(isttyp,:),'linewidth',.5);   hold on
        end
        xlim([flo,1]);
        ylim([-150 130])
        text(1/1000,150,suff,'FontSize',8)
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
            ylabel('Pa^2/Hz, dB')
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
filename=sprintf('%s/NoiseSpec_StaType_DepZ.pdf',figoutpath);
print(gcf,'-dpdf',filename)

figure(2)

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 11]);
filename=sprintf('%s/NoiseSpec_StaType_DepH.pdf',figoutpath);
print(gcf,'-painters','-dpdf',filename)


figure(3)

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 11]);
filename=sprintf('%s/NoiseSpec_StaType_DepAPG.pdf',figoutpath);
print(gcf,'-painters','-dpdf',filename)


figure(4)

set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 11]);
filename=sprintf('%s/NoiseSpec_StaType_DepDPG.pdf',figoutpath);
print(gcf,'-painters','-dpdf',filename)

