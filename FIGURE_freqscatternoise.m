% FIGURE_freqscatternoise

clear;
close all

OBS_TableParams;

flo = 1/1000;
fhi = 2.5;

load(mattable);
OBS_table_orig = OBS_table;

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
%     for iexp = [1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18];
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
        
        elevvec(ie) = waterdepth(ista);
        % APG v DPG Check
        if strcmp(prestypes(ista),'APG')==1
            isAPG(ie) = 1;
        else
            isAPG(ie) = 0;
        end
        ndays(ie) = avgprop.params.ndays;
        seiscol(ie,:) = scg(find(strcmp(statypes,statype(ista))==1),:);
        stavec(ie) = statype(ista);
        smtvec(ie) = seismometer(ista);
  
        ie = ie+1;
  
    end
end

nsta = length(elevvec);

for ifreq = 1:length(flo_vec) % switch to scatter if possible...
    flo = flo_vec(ifreq);
    fhi = fhi_vec(ifreq);
    idxf = find(fc<=fhi & fc>=flo);
    fnew = fc(idxf);
    [LL,HH,FF] = noise_models(100);
    idxfm = find(FF<=fhi & FF>=flo);
    LAVG = mean(LL(idxfm));
    HAVG = mean(HH(idxfm));
    n = 1;
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
    if meannoise(ie)<HAVG;
        n=n+1;
    end
    ie = ie+1;
    end
        
    figure(1)
    subplot(length(flo_vec),1,ifreq);
    scatter((watvec),meannoise,10,secol,char(seissym(iseis))); hold on
    idxp = isnan(meannoise);
    end
%     plot([0 6500], [LAVG LAVG],'-k','linewidth',1);
%     plot([0 6500], [HAVG HAVG],'-k','linewidth',1);
    ylim([-190 -60])
%     xlim([0 6500])
    set(gca,'FontSize',8)
    box on
xlabel('Water Depth (m)')
ylabel('Power (m^2/s^4/Hz, dB)');

    
    ii = 2;
    for iseis = 1:length(seistypes)
            idxseis = find(strcmp(seistypes(iseis),smtvec));
            clear meannoise watvec secol
            ie =1;
    for ista = 1:nsta
        if isempty(intersect(ista,idxseis))==1
            continue
        end
    noise = smooth(spc_stack(idxf,ista,ii),40);
    meannoise(ie) = mean(noise);
    watvec(ie)= elevvec(ista);
    secol(ie,:) = seiscol(ista,:);
    ie = ie+1;
    end
    
    figure(2)
    subplot(length(flo_vec),1,ifreq);
    scatter((watvec),meannoise,10,secol,char(seissym(iseis))); hold on
    end
    
    ii = 3;
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
    
    figure(2)
    subplot(length(flo_vec),1,ifreq);
    scatter((watvec),meannoise,10,secol,char(seissym(iseis))); hold on
    end
%     plot([0 6500], [LAVG LAVG],'-k','linewidth',1);
%     plot([0 6500], [HAVG HAVG],'-k','linewidth',1);
    ylim([-190 -60])
%     xlim([0 6500])
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
filename=sprintf('%s/scatterucorr_Z.pdf',figoutpath);
print(gcf,'-dpdf',filename)

figure(2)
set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[8.5 11]);
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 3 6]);
filename=sprintf('%s/scatterucorr_H.pdf',figoutpath);
print(gcf,'-dpdf',filename)

