% FIGURE_coherences

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
        
        if Zisgood(ista)==0
            coh_stack(:,ie,1) = NaN(size(fc));
            coh_stack(:,ie,2) = NaN(size(fc));
            coh_stack(:,ie,3) = NaN(size(fc));
        else
            % CAREFUL THIS IS CURRENTLY SET UP FOR THE PRESSURES (E.G. SUPP
            % FIGURE, NOT MAIN TEXT). HARDWIRED, NEED TO MANUALLY ADJUST
            coh1z = abs(staavg.cross.c1p_mean).^2./(staavg.power.c11_mean.*staavg.power.cpp_mean);
            coh1zsmooth = smoothSpectrum_octave(coh1z,avgprop.params.f,fc_sta,ocav);
            coh_stack(:,ie,1) = interp1(fc_sta,coh1zsmooth,fc);
            
            coh2z = abs(staavg.cross.c2p_mean).^2./(staavg.power.c22_mean.*staavg.power.cpp_mean);
            coh2zsmooth = smoothSpectrum_octave(coh2z,avgprop.params.f,fc_sta,ocav);
            coh_stack(:,ie,2) = interp1(fc_sta,coh2zsmooth,fc);
            
            cohpz = abs(staavg.cross.c12_mean).^2./(staavg.power.c11_mean.*staavg.power.c22_mean);
            cohpzsmooth = smoothSpectrum_octave(cohpz,avgprop.params.f,fc_sta,ocav);
            coh_stack(:,ie,3) = interp1(fc_sta,cohpzsmooth,fc);
        end
        if H1isgood(ista)==0
            coh_stack(:,ie,1) = NaN(size(fc));
        end
        if H2isgood(ista)==0
            coh_stack(:,ie,2) = NaN(size(fc));             
        end
        if Pisgood(ista)==0
            coh_stack(:,ie,3) = NaN(size(fc));           
        end
      
        
        elev_vec(ie) = waterdepth(ista);
        % APG v DPG Check
        if strcmp(prestypes(ista),'APG')==1
            isAPG(ie) = 1;
        else
            isAPG(ie) = 0;
        end
        ndays(ie) = avgprop.params.ndays;
        
        ie = ie+1;
    end
    
end

fnew = fc;
nsta = length(elev_vec);

% make figures
for ii = 1:3
    % set up figure naming properties
    if ii == 1;
        suff = 'Z1';
    elseif ii == 2;
        suff = 'Z2';
    else
        suff = 'ZP';
    end
    
    cohere = coh_stack(:,:,ii);
    flimvec1 = sqrt(9.8./(2*pi*[0:1:6500]));
    flimvec2 = sqrt(9.8./(1.6*pi*[0:1:6500]));
    flimvec3 = 1.5./sqrt([0:1:6500]);
    
    figure(ii); clf
    for ista = 1:nsta
        scatter(fnew,elev_vec(ista).*ones(size(fnew)),7,(1-(cohere(:,ista))).*ones(length(fnew),3),'filled'); hold on
    end
    set(gca,'XScale','log');
    xlim([1/1000,1]);
    plot([1/7,1/7],[0,7000],'-r');
    plot([1/14,1/14],[0,7000],'-b');
    plot([1/10,1/10],[0,7000],'-k');
    plot(flimvec1,[0:6500],'-g');
    plot(flimvec2,[0:6500],'-c');
    plot(flimvec3,[0:6500],'-y');
    ylim([0,6500]);
    xlabel('Frequency (Hz)');
    ylabel('Elevation (m)');
    figure(ii)
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPosition',[.05 .05 4 4.5]);
    
        filename=sprintf('%s/CoherencesALT_%s.pdf',figoutpath,suff);
exportgraphics(gcf,filename,'ContentType','vector')
end

