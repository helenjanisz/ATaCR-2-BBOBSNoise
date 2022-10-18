% FIGURE_coherences_postremoval
% FIGURE_coherences

clear; 
% close all

OBS_TableParams;

flo = 1/1000;
fhi = 2.5;

load(mattable);
OBS_table_orig = OBS_table;

%DO NOT EDIT ORDER
TF_label_list = {'Z2-1','ZP'}; % check coherence after removing tilt and removing compliance

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
            coh_stack(:,ie,1) = NaN(size(fc)); %ZP coh - tilt
            coh_stack(:,ie,2) = NaN(size(fc)); %Z1 coh - comp
            coh_stack(:,ie,3) = NaN(size(fc)); %Z2 coh - comp
            isbad=1;
        end
        
        % indcies are hardwired do not edit
        % ZP coh - tilt
        
        if it ==1
            cohpz = corr(idxx(it)).cohere;
            cohpzsmooth = smoothSpectrum_octave(cohpz,f',fc_sta,ocav);
            coh_stack(:,ie,1) = interp1(fc_sta,cohpzsmooth,fc);
            if H1isgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end
            if H2isgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end
            if Pisgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end

        % ZH coh - compliance
        elseif it ==2
            coh1z = corr(idxx(it)).cohere1;
            coh1zsmooth = smoothSpectrum_octave(coh1z,f',fc_sta,ocav);
            coh_stack(:,ie,2) = interp1(fc_sta,coh1zsmooth,fc);
            if H1isgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end
            if Pisgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end
            
            coh2z = corr(idxx(it)).cohere2;
            coh2zsmooth = smoothSpectrum_octave(coh2z,f',fc_sta,ocav);
            coh_stack(:,ie,3) = interp1(fc_sta,coh2zsmooth,fc);
            if H2isgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end
            if Pisgood(ista) ~=1
                coh_stack(:,ie,1) = NaN(size(fc));
            end
        end

        end
        elev_vec(ie) = waterdepth(ista);
        
        ie = ie+1;
    end
    
end

fnew = fc;
nsta = length(elev_vec);

% make figures
for ii = 1:3
    % set up figure naming properties
    if ii == 1;
        suff = 'ZP_Tilt';
    elseif ii == 2;
        suff = 'Z1_Compliance';
    else
        suff = 'Z2_Compliance';
    end
    flimvec1 = sqrt(9.8./(2*pi*[0:1:6500]));
    flimvec2 = sqrt(9.8./(1.6*pi*[0:1:6500]));
    flimvec3 = 1.5./sqrt([0:1:6500]);
    cohere = coh_stack(:,:,ii);
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
    colormap(gray);
caxis([0 1]);
colorbar
    set(gcf,'PaperPositionMode','manual');
    set(gcf,'PaperUnits','inches');
    set(gcf,'PaperOrientation','portrait');
    set(gcf,'PaperPosition',[.05 .05 4 4.5]);
    
        filename=sprintf('%s/Coherences_Corr_%s.pdf',figoutpath,suff);
exportgraphics(gcf,filename,'ContentType','vector')
end

