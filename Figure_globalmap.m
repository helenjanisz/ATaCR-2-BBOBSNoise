%FIGURE_globalmap

clear; close all

OBS_TableParams;
load(mattable);

stations = OBS_table.Station;
networks = OBS_table.Network;
expfldr = OBS_table.ExperimentFolder;
expabrv = OBS_table.ExperimentAbbreviation;
waterdepth = OBS_table.WaterDepth;
isgood = OBS_table.IsGood;
statype = OBS_table.InstrumentDesign;
seismometer = OBS_table.Seismometer;
prestypes = OBS_table.PressureGauge;
lats = OBS_table.Latitude;
lons = OBS_table.Longitude;

experiments_all = unique(expfldr);
pressure_all = unique(prestypes);

ii = 1;
for ista = 1:length(stations);
    if isgood(ista) == 0
        continue
    end
    stalat(ii) = lats(ista);
    if lons(ista) <0
        stalon(ii) = 360+lons(ista);
    else
        stalon(ii) = lons(ista);
    end
    ii = ii+1;
end

% make global map
figure(1)
clf;
hold on
set(gcf,'color','white')
m_proj('robinson','lon',[0 360]);
% ax1=axes;
colormap([m_colmap('blues')]);
caxis([-10000 0]);
[CS,CH]=m_etopo2('contourf','edgecolor','none');
m_coast('patch',[.7 .7 .7],'edgecolor','none');
m_grid('linestyle','none','tickdir','out');
m_plot(stalon,stalat,'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',2)
m_text(290,35,'ENAM','color','r','fontsize',8);
m_text(240,37,'ALBACORE','color','r','fontsize',8);
m_text(245,20,'SCOOBA','color','r','fontsize',8);
m_text(238,45,'CK, CI, GORDA, BLANCO','color','r','fontsize',8);
m_text(212,55,'AACSE','color','r','fontsize',8);
m_text(212,20,'PLUME','color','r','fontsize',8);
m_text(210,3,'NOMELT','color','r','fontsize',8);
m_text(225,-12,'YOUNG ORCA','color','r','fontsize',8);
m_text(155,33,'PLATE','color','r','fontsize',8);
m_text(150,20,'MARIANA','color','r','fontsize',8);
m_text(152,-10,'PAPUA','color','r','fontsize',8);
m_text(188,-20,'LAU','color','r','fontsize',8);
m_text(182,-40,'HOBITSS','color','r','fontsize',8);
m_text(175,-50,'MOANA','color','r','fontsize',8);
m_text(37,-12,'SEGMENT','color','r','fontsize',8);


set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperOrientation','portrait');
set(gcf,'PaperPosition',[.05 .05 8 8]);
filename=sprintf('%s/StationMap.pdf',figoutpath);
print(gcf,'-dpdf',filename)