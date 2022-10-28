load /Users/Seth/Desktop/2022_diagnostic_Ni_model/allJs.mat
load /Users/Seth/Desktop/2022_diagnostic_Ni_model/allSols.mat
load /Users/Seth/Desktop/2022_diagnostic_Ni_model/allPres.mat
% load /Users/Seth/Desktop/2022_diagnostic_Ni_model/WOA18Si.mat
% load /Users/Seth/Desktop/2022_diagnostic_Ni_model/WOA18P.mat
load Niclimatology.mat
load gridding/GT_IDP21_Ni.mat
load ../data/ao.mat
addpath '../../../Matlab_Utils_Data/utils/cptcmap'

% make a new colormap based on balance, but with offset zero
newmap=cmocean('balance',1000);
% newmap=cptcmap('GMT_polar','ncol',1000);
% newcmap=colormap(brewermap('rdbu'),1000);
iremove=(500:2:1000);
newmap(iremove,:)=[];
iremove=(500:2:748);
newmap(iremove,:)=[];
iremove=(500:2:624);
newmap(iremove,:)=[];

% grab the obs
NIOBS = GT_IDP21_Ni;
niobs = NIOBS(ao.iocn);

% Calculate the transport convergences and errors for plotting
J_Ni = repmat(ao.nanOCN,1,1,1,size(allJs,2));
for i=1:size(allJs,2);
    J=ao.nanOCN; J(ao.iocn)=allJs(:,i,1); J(82:91,:,:)=0;
    J_Ni(:,:,:,i)=J;
end
J_Ni_surf = nanmean(J_Ni(:,:,1:3,:),3);
nitransects = squeeze(nanmean(J_Ni_surf,2));
nitransect = -nanmean(nitransects,2); nitransect(isnan(nitransect))=0;
err_nitrans = std(nitransects,0,2); err_nitrans(isnan(err_nitrans))=0;
niprofs = squeeze(nanmean(nanmean(J_Ni,1),2));
niprof = nanmean(niprofs,2); niprof=movmean(niprof,3,1);
err_niprof = std(niprofs,0,2); err_niprof=movmean(err_niprof,3,1);
cumni = cumsum(niprof(4:24)'.*ao.height(4:24))/max(cumsum(niprof(4:24)'.*ao.height(4:24)));

J_P = repmat(ao.nanOCN,1,1,1,size(allJs,2));
for i=1:size(allJs,2);
    J=ao.nanOCN; J(ao.iocn)=allJs(:,i,2); J(82:91,:,:)=0;
    J_P(:,:,:,i)=J*1000;
end
j_P_surf = squeeze(nanmean(J_P(:,:,1:3,:),3));
ptransects = squeeze(nanmean(j_P_surf,2));
ptransect = -nanmean(ptransects,2); ptransect(isnan(ptransect))=0;
err_ptrans = std(ptransects,0,2); err_ptrans(isnan(err_ptrans))=0;
pprofs = squeeze(nanmean(nanmean(J_P,1),2));
pprof = nanmean(pprofs,2); pprof=movmean(pprof,3,1);
err_pprof = std(pprofs,0,2); err_pprof=movmean(err_pprof,3,1);
cump = cumsum(pprof(4:24)'.*ao.height(4:24))/max(cumsum(pprof(4:24)'.*ao.height(4:24)));

J_Si = repmat(ao.nanOCN,1,1,1,size(allJs,2));
for i=1:size(allJs,2);
    J=ao.nanOCN; J(ao.iocn)=allJs(:,i,3); J(82:91,:,:)=0;
    J_Si(:,:,:,i)=J*1000;
end
J_Si_surf = nanmean(J_Si(:,:,1:3,:),3);
sitransects = squeeze(nanmean(J_Si_surf,2));
sitransect = -nanmean(sitransects,2); sitransect(isnan(sitransect))=0;
err_sitrans = std(sitransects,0,2); err_sitrans(isnan(err_sitrans))=0;
siprofs = squeeze(nanmean(nanmean(J_Si,1),2));
siprof = nanmean(siprofs,2); siprof=movmean(siprof,3,1);
err_siprof = std(siprofs,0,2); err_siprof=movmean(err_siprof,3,1);

%% calculate cumulative sums 
cumni = cumsum(niprof(4:24)'.*ao.height(4:24))/max(cumsum(niprof(4:24)'.*ao.height(4:24)));
cump = cumsum(pprof(4:24)'.*ao.height(4:24))/max(cumsum(pprof(4:24)'.*ao.height(4:24)));
cumsi = cumsum(siprof(4:24)'.*ao.height(4:24))/max(cumsum(siprof(4:24)'.*ao.height(4:24)));

cumnidepths = interp1(ao.depth(4:24),cumni,[1000 2000 3000 4000]);
cumpdepths = interp1(ao.depth(4:24),cump,[1000 2000 3000 4000]);
cumsidepths = interp1(ao.depth(4:24),cumsi,[1000 2000 3000 4000]);

%% calculate ratios
% niptrans=nitrans./ptrans; niptrans(isnan(niptrans))=0;
% err_niptrans=niptrans.*sqrt(((err_nitrans./nitrans).^2)+((err_ptrans./ptrans).^2)); err_niptrans(isnan(err_niptrans))=0;
% nisitrans=nitrans./sitrans; nisitrans(isnan(nisitrans))=0;
% err_nisitrans=nisitrans.*sqrt(((err_nitrans./nitrans).^2)+((err_sitrans./sitrans).^2)); err_nisitrans(isnan(err_nisitrans))=0;
% 
% nipprof=niprof./pprof; nipprof(isnan(nipprof))=0;
% err_nipprof=nipprof.*sqrt(((err_niprof./niprof).^2)+((err_pprof./pprof).^2)); err_nipprof(isnan(err_nipprof))=0;
% nisiprof=niprof./siprof; nisiprof(isnan(nisiprof))=0;
% err_nisiprof=nisiprof.*sqrt(((err_niprof./niprof).^2)+((err_siprof./siprof).^2)); err_nisiprof(isnan(err_nisiprof))=0;

%% process the restoring-smoothed distributions
Ni = mean(allSols(:,:,1),2);
NI = ao.nanOCN;
NI(ao.iocn)=Ni;
NI(82:91,:,:)=0;

p = mean(allSols(:,:,2),2)*1000;
P = ao.nanOCN;
P(ao.iocn)=p;
P(82:91,:,:)=0;

Si = mean(allSols(:,:,3),2)*1000;
SI = ao.nanOCN;
SI(ao.iocn)=Si;
SI(82:91,:,:)=0;

%% load the original neural net Ni, and WOAobs

mlrni=Niclimatology;
MLRNI = ao.nanOCN;
MLRNI(ao.iocn)=mlrni;

% woaSi = WOA18Si*1000;
% WOASI = ao.nanOCN;
% WOASI(ao.iocn)=woaSi;
% WOASI(82:91,:,:)=0;
% 
% woaP = WOA18P*1000;
% WOAP = ao.nanOCN;
% WOAP(ao.iocn)=woaP;
% WOAP(82:91,:,:)=0;

%% create the obs and model distributions with depth
% mlrni0100 = nanmean(MLRNI(:,:,1:3),3);
% mlrni100500 = nanmean(MLRNI(:,:,4:8),3);
% mlrni5001500 = nanmean(MLRNI(:,:,9:14),3);
% mlrni15004500 = nanmean(MLRNI(:,:,14:21),3);
% obsni0100 = nanmean(NIOBS(:,:,1:3),3); i0100=find(~isnan(obsni0100));
% obsni100500 = nanmean(NIOBS(:,:,4:8),3); i100500=find(~isnan(obsni100500));
% obsni5001500 = nanmean(NIOBS(:,:,9:14),3); i5001500=find(~isnan(obsni5001500));
% obsni15004500 = nanmean(NIOBS(:,:,14:21),3); i15004500=find(~isnan(obsni15004500));

%% find the low-latitude Ni flux profiles and the high-latitude profiles,
% low lat is only below 40S, mid lat is 40S to 40N

J_Ni_LOLAT=J_Ni; J_Ni_LOLAT(26:91,:,:)=NaN;
niprof_lolat = squeeze(nanmean(nanmean(J_Ni_LOLAT,1),2));
niprof_lolat = movmean(niprof_lolat,3);
J_Ni_MIDLAT=J_Ni; J_Ni_MIDLAT(1:25,:,:)=NaN; J_Ni_MIDLAT(67:91,:,:)=NaN;
niprof_midlat = squeeze(nanmean(nanmean(J_Ni_MIDLAT,1),2));
niprof_midlat = movmean(niprof_midlat,3);

J_P_LOLAT=J_P; J_P_LOLAT(26:91,:,:)=NaN;
pprof_lolat = squeeze(nanmean(nanmean(J_P_LOLAT,1),2));
pprof_lolat = movmean(pprof_lolat,3);
J_P_MIDLAT=J_P; J_P_MIDLAT(1:25,:,:)=NaN; J_P_MIDLAT(67:91,:,:)=NaN;
pprof_midlat = squeeze(nanmean(nanmean(J_P_MIDLAT,1),2));
pprof_midlat = movmean(pprof_midlat,3);

J_Si_LOLAT=J_Si; J_Si_LOLAT(26:91,:,:)=NaN;
siprof_lolat = squeeze(nanmean(nanmean(J_Si_LOLAT,1),2));
siprof_lolat = movmean(siprof_lolat,3);
J_Si_MIDLAT=J_Si; J_Si_MIDLAT(1:25,:,:)=NaN; J_Si_MIDLAT(67:91,:,:)=NaN;
siprof_midlat = squeeze(nanmean(nanmean(J_Si_MIDLAT,1),2));
siprof_midlat = movmean(siprof_midlat,3);

%% load and process the preformed nutrients

NOARCTIC=ao.nanOCN; NOARCTIC(82:91,:,:)=NaN;

ni_pre = nanmean(allPres(:,:,1),2);
NI_PRE=ao.nanOCN; NI_PRE(ao.iocn)=ni_pre; NI_PRE=NI_PRE.*NOARCTIC;
nipre_prof = squeeze(nansum(nansum(NI_PRE.*ao.VOL.*NOARCTIC))./nansum(nansum(ao.VOL.*NOARCTIC)));
ni_post = nanmean(allSols(:,:,1),2);
NI_POST=ao.nanOCN; NI_POST(ao.iocn)=ni_post; NI_POST=NI_POST.*NOARCTIC;
nipost_prof = squeeze(nansum(nansum(NI_POST.*ao.VOL.*NOARCTIC))./nansum(nansum(ao.VOL.*NOARCTIC)));
nireg_prof = nipost_prof-nipre_prof;

p_pre = nanmean(allPres(:,:,2),2)*1000;
P_PRE=ao.nanOCN; P_PRE(ao.iocn)=p_pre; P_PRE=P_PRE.*NOARCTIC;
ppre_prof = squeeze(nansum(nansum(P_PRE.*ao.VOL.*NOARCTIC))./nansum(nansum(ao.VOL.*NOARCTIC)));
p_post = nanmean(allSols(:,:,2),2)*1000;
P_POST=ao.nanOCN; P_POST(ao.iocn)=p_post; P_POST=P_POST.*NOARCTIC;
ppost_prof = squeeze(nansum(nansum(P_POST.*ao.VOL.*NOARCTIC))./nansum(nansum(ao.VOL.*NOARCTIC)));
preg_prof = ppost_prof-ppre_prof;

si_pre = nanmean(allPres(:,:,3),2)*1000;
SI_PRE=ao.nanOCN; SI_PRE(ao.iocn)=si_pre; SI_PRE=SI_PRE.*NOARCTIC;
sipre_prof = squeeze(nansum(nansum(SI_PRE.*ao.VOL.*NOARCTIC))./nansum(nansum(ao.VOL.*NOARCTIC)));
si_post = nanmean(allSols(:,:,3),2)*1000;
SI_POST=ao.nanOCN; SI_POST(ao.iocn)=si_post; SI_POST=SI_POST.*NOARCTIC;
sipost_prof = squeeze(nansum(nansum(SI_POST.*ao.VOL.*NOARCTIC))./nansum(nansum(ao.VOL.*NOARCTIC)));
sireg_prof = sipost_prof-sipre_prof;
%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 7alt, heatmap of relative uptake of Ni compared to P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surfupni=(J_Ni(1:81,:,1)); surfupni=surfupni(~isnan(surfupni));
surfupp=(J_P(1:81,:,1)); surfupp=surfupp(~isnan(surfupp));
% surfupsi=(SURFUPSI(:,:,1)); surfupsi=surfupsi(~isnan(surfupsi));
surfni=squeeze(MLRNI(1:81,:,1)); surfni=surfni(~isnan(surfni));
ni2p = surfupni./surfupp; ni2p(ni2p>15)=NaN; ni2p(ni2p<0)=NaN;

xbins = 0:0.1:8;
ybins = 0:0.1:12;
data = hist3([ni2p,surfni],'edges',{ybins,xbins}); data = data/(max(max(data)))*100;
% data=smoothdata(data,1,'gaussian',2);
% data=smoothdata(data,2,'gaussian',2);


figure(7); clf; papersize = [12 12];
set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
colormap(flipud(hot));
plot(surfni,ni2p,'.k','markersize',1); hold on;
% contourf(xbins,ybins,data,50,'linestyle','none'); hold on;
surf(xbins,ybins,data,'linestyle','none','facealpha',0.5); hold on;
view([0 90])
% set(gca,'ylim',[0 5],'xlim',[0 12.5],'xtick',[0.5 2.5 4.5 6.5 8.5 10.5 12.5],'xticklabel',{'2 nM','3 nM','4 nM','5 nM','6 nM','7 nM','8 nM'},'fontsize',10,'position',[0.2 0.2 0.7 0.7])
xlabel('Surface Ni concentration','fontsize',12,'position',[6 -1])
ylabel('Surface Ni uptake/P uptake (nM/\muM)','fontsize',12,'position',[-2 2.5])

pdfname=['/Users/Seth/Desktop/temp/fig7alt.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)


% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Figure 8 preformed and regenerated nutrients
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure(8); clf; papersize = [30 15];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% subplot(2,4,1)
% h1=fill([ppost_prof;zeros(24,1)],[ao.depth'/1000;flipud(ao.depth'/1000)],[1 .6 .7],'linestyle','none'); hold on;
% h2=plot(ppost_prof,ao.depth/1000,'color',[.2 0 0],'Linewidth',1); hold on;
% h3=fill([ppre_prof;zeros(24,1)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.8 .2 .3],'linestyle','none'); hold on;
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 3],'XAxisLocation','top','Position',[0.05 .1 .18 .7])
% xlabel('Phosphate (\muM)')
% ylabel('Depth (km)')
% legend([h3 h1 h2],'Preformed P','Regenerated P','Total P','Location','SouthOutside')
% 
% subplot(2,4,2)
% h1=fill([sipost_prof;zeros(24,1)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.7 1 .6],'linestyle','none'); hold on;
% h2=plot(sipost_prof,ao.depth/1000,'color',[0 .2 0],'Linewidth',1); hold on;
% h3=fill([sipre_prof;zeros(24,1)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.2 .8 .1],'linestyle','none'); hold on;
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 150],'XAxisLocation','top','Position',[.28 .1 .18 .7])
% xlabel('Silicate (\muM)')
% legend([h3 h1 h2],'Preformed Si','Regenerated Si','Total Si','Location','SouthOutside')
% 
% subplot(2,4,3)
% h1=fill([nipost_prof;zeros(24,1)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.6 .7 1],'linestyle','none'); hold on;
% h2=plot(nipost_prof,ao.depth/1000,'color',[0 0 .2],'Linewidth',1); hold on;
% h3=fill([nipre_prof;zeros(24,1)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.2 .3 .8],'linestyle','none'); hold on;
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 10],'XAxisLocation','top','Position',[.51 .1 .18 .7])
% xlabel('Nickel (nM)')
% legend([h3 h1 h2],'Preformed Ni','Regenerated Ni','Total Ni','Location','SouthOutside')
% 
% subplot(2,4,4)
% plot(preg_prof/max(preg_prof),ao.depth/1000,'color',[.8 0 0],'Linewidth',1); hold on;
% plot(sireg_prof/max(sireg_prof),ao.depth/1000,'color',[0 .8 0],'Linewidth',1); hold on;
% plot(nireg_prof/max(nireg_prof),ao.depth/1000,'color',[0 0 .8],'Linewidth',1); hold on;
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 1.1],'XAxisLocation','top','Position',[.74 .1 .18 .7])
% xlabel('Normalized concentration')
% legend('Regenerated P','Regenerated Si','Regenerated Ni','Location','SouthOutside')
% 
% pdfname=['/Users/Seth/Desktop/temp/fig8.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 7 histogram of relative uptake of Ni compared to P
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% surfupni=(J_Ni(1:81,:,1)); surfupni=surfupni(~isnan(surfupni));
% surfupp=(J_P(1:81,:,1)); surfupp=surfupp(~isnan(surfupp));
% % surfupsi=(SURFUPSI(:,:,1)); surfupsi=surfupsi(~isnan(surfupsi));
% surfni=squeeze(MLRNI(1:81,:,1)); surfni=surfni(~isnan(surfni));
% 
% [n,edges,i]=histcounts(surfni,(2:.5:8));
% ni2p = surfupni./surfupp; ni2p(ni2p>15)=NaN; ni2p(ni2p<0)=NaN;
% % ni2p = ni2p*(nansum(nansum(nansum(WOA18P)))/nansum(nansum(nansum(NI))));
% 
% figure(7); clf; papersize = [12 12];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% boxplot(ni2p,i,'Symbol',''); hold on;
% % set(gca,'Labels',{'2-2.5 nM','2.5-3 nM','3-3.5 nM','3.5-4 nM','4-4.5 nM','5-5.5 nM','5.5-6 nM','6-6.5 nM','6.5-7 nM','7-7.5 nM','7.5-8 nM','2-3 nM'})
% set(gca,'ylim',[0 5],'xlim',[0 12.5],'xtick',[0.5 2.5 4.5 6.5 8.5 10.5 12.5],'xticklabel',{'2 nM','3 nM','4 nM','5 nM','6 nM','7 nM','8 nM'},'fontsize',10,'position',[0.2 0.2 0.7 0.7])
% xtickangle(60);
% xlabel('Surface Ni concentration','fontsize',12,'position',[6 -1])
% ylabel('Surface Ni uptake/P uptake (nM/\muM)','fontsize',12,'position',[-2 2.5])
% 
% pdfname=['/Users/Seth/Desktop/temp/fig7.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 6 Nickel model compared to obs at different depth surfaces
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(6); clf; papersize = [45 30];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% surfmap = ao.OCN(:,:,1);
% [LON,LAT] = meshgrid(ao.lon,ao.lat);
% 
% subplot(2,2,1);
% worldmap([-80 70],[0 360])
% setm(gca,'MapProjection','miller','MLabelParallel','south','FontSize',12)
% contourfm(ao.LAT(:,:,1),ao.LON(:,:,1),mlrni0100,50,'linestyle','none'); hold on;
% scatterm(LAT(i0100),LON(i0100),20,obsni0100(i0100),'filled','k'); hold on;
% scatterm(LAT(i0100),LON(i0100),14,obsni0100(i0100),'filled'); hold on;
% geoshow('landareas.shp', 'FaceColor', [0.7 0.7 0.7]); hold on;
% set(gca,'clim',[0 10],'ytick',[-60 -30 0 30 60])
% cmocean('deep')
% c=colorbar; title(c,'Ni (nM)','fontsize',12);
% title('0-100 m','fontsize',18)
% 
% subplot(2,2,2);
% worldmap([-80 70],[0 360])
% setm(gca,'MapProjection','miller','MLabelParallel','south','FontSize',12)
% contourfm(ao.LAT(:,:,1),ao.LON(:,:,1),mlrni100500,50,'linestyle','none'); hold on;
% scatterm(LAT(i100500),LON(i100500),20,obsni100500(i100500),'filled','k'); hold on;
% scatterm(LAT(i100500),LON(i100500),14,obsni100500(i100500),'filled'); hold on;
% geoshow('landareas.shp', 'FaceColor', [0.7 0.7 0.7]); hold on;
% set(gca,'clim',[0 10],'ytick',[-60 -30 0 30 60])
% cmocean('deep')
% c=colorbar; title(c,'Ni (nM)','fontsize',12);
% title('100-500 m','fontsize',18)
% 
% subplot(2,2,3);
% worldmap([-80 70],[0 360])
% setm(gca,'MapProjection','miller','MLabelParallel','south','FontSize',12)
% contourfm(ao.LAT(:,:,1),ao.LON(:,:,1),mlrni5001500,50,'linestyle','none'); hold on;
% scatterm(LAT(i5001500),LON(i5001500),20,obsni5001500(i5001500),'filled','k'); hold on;
% scatterm(LAT(i5001500),LON(i5001500),14,obsni5001500(i5001500),'filled'); hold on;
% geoshow('landareas.shp', 'FaceColor', [0.7 0.7 0.7]); hold on;
% set(gca,'clim',[0 10],'ytick',[-60 -30 0 30 60])
% cmocean('deep')
% c=colorbar; title(c,'Ni (nM)','fontsize',12);
% title('500-1500 m','fontsize',18)
% 
% subplot(2,2,4);
% worldmap([-80 70],[0 360])
% setm(gca,'MapProjection','miller','MLabelParallel','south','FontSize',12)
% contourfm(ao.LAT(:,:,1),ao.LON(:,:,1),mlrni15004500,50,'linestyle','none'); hold on;
% scatterm(LAT(i15004500),LON(i15004500),20,obsni15004500(i15004500),'filled','k'); hold on;
% scatterm(LAT(i15004500),LON(i15004500),14,obsni15004500(i15004500),'filled'); hold on;
% geoshow('landareas.shp', 'FaceColor', [0.7 0.7 0.7]); hold on;
% set(gca,'clim',[0 10],'ytick',[-60 -30 0 30 60])
% cmocean('deep')
% c=colorbar; title(c,'Ni (nM)','fontsize',12);
% title('1500-4500 m','fontsize',18)
% 
% pdfname=['/Users/Seth/Desktop/temp/fig6.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 5 correlation between model and obs
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(5); clf; papersize = [13 11];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% subplot(1,1,1)
% % scatter (niobs,mlrni,20,ao.Depth)
% % set(gca,'ColorScale','log','clim',[10 10000])
% % cb = colorbar('xtick',[10 100 1000],'xTickLabel',{'10' '100' '1000'});
% 
% 
% bins = 0:0.2:12;
% data = hist3([mlrni,niobs],'edges',{bins,bins}); data = data/(max(max(data)))*100;
% data=smoothdata(data,1,'gaussian',2);
% data=smoothdata(data,2,'gaussian',2);
% contourf(bins,bins,data,50,'linestyle','none'); hold on;
% % scatter(mlrni,niobs,'.k')
% line([0 12],[0 12],'color','black','linestyle','--'); hold on;
% crameri('-oslo')
% % contour(bins,bins,data,(0:10:100),'-k'); hold on;
% set(gca,'clim',[0 100])
% cb = colorbar;
% % text(4.5,13, ['\it R^2 = ' num2str(r2,2)],'fontweight','bold','fontsize',14);
% xlabel('observed [Ni]')
% ylabel ('model [Ni]');
% 
% pdfname=['/Users/Seth/Desktop/temp/fig5.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 4 original versus restoring-smoothed sections at 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(4); clf; papersize = [35 25];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% subplot(3,3,1)
% pcolor(ao.lat,ao.depth/1000,squeeze(WOAP(:,100,:))'); shading flat;
% set(gca,'clim',[0 3.2],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.05 .7 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('deep')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('World Ocean Atlas 2018 PO_4^- (\muM)','FontWeight','Normal')
% 
% subplot(3,3,4)
% pcolor(ao.lat,ao.depth/1000,squeeze(P(:,100,:))'); shading flat;
% set(gca,'clim',[0 3.2],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.05 .4 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('deep')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('Restoring-smoothed PO_4^- (\muM)','FontWeight','Normal')
% 
% subplot(3,3,2)
% pcolor(ao.lat,ao.depth/1000,squeeze(WOASI(:,100,:))'); shading flat;
% set(gca,'clim',[0 180],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.35 .7 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('deep')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('World Ocean Atlas 2018 Si (\muM)','FontWeight','Normal')
% 
% subplot(3,3,5)
% pcolor(ao.lat,ao.depth/1000,squeeze(SI(:,100,:))'); shading flat;
% set(gca,'clim',[0 180],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.35 .4 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('deep')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('Restoring-smoothed Si (\muM)','FontWeight','Normal')
% 
% subplot(3,3,3)
% pcolor(ao.lat,ao.depth/1000,squeeze(MLRNI(:,100,:))'); shading flat;
% set(gca,'clim',[0 10],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.65 .7 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('deep')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('MLR predicted Ni (nM)','FontWeight','Normal')
% 
% subplot(3,3,6)
% pcolor(ao.lat,ao.depth/1000,squeeze(NI(:,100,:))'); shading flat;
% set(gca,'clim',[0 10],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.65 .4 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('deep')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('Restoring-smoothed Ni (nM)','FontWeight','Normal')
% 
% subplot(3,3,7)
% pcolor(ao.lat,ao.depth/1000,(squeeze(WOAP(:,100,:))'-squeeze(P(:,100,:))')); shading flat;
% set(gca,'clim',[-.5 .5],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.05 .1 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('diff')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('WOA 2018 PO_4^- - smoothed PO_4^-','FontWeight','Normal')
% 
% subplot(3,3,8)
% pcolor(ao.lat,ao.depth/1000,(squeeze(WOASI(:,100,:))'-squeeze(SI(:,100,:))')); shading flat;
% set(gca,'clim',[-20 20],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.35 .1 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('diff')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('WOA 2018 PO_4^- - smoothed PO_4^-','FontWeight','Normal')
% 
% subplot(3,3,9)
% pcolor(ao.lat,ao.depth/1000,(squeeze(MLRNI(:,100,:))'-squeeze(NI(:,100,:))')); shading flat;
% set(gca,'clim',[-1 1],'xlim',[-80 60],'ydir','reverse','ylim',[0 5.6],'color',[0.6 0.6 0.6],'position',[.65 .1 .25 .2])
% % cptcmap('GMT_haxby')
% cmocean('diff')
% colorbar
% xlabel('Latitude')
% ylabel('Depth (km)')
% title('NN Ni - smoothed Ni','FontWeight','Normal')
% 
% pdfname=['/Users/Seth/Desktop/temp/fig4.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Alternative figure 3 flux profiles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3); clf; papersize = [15 20];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% subplot(3,3,1)
% plot(pprof,ao.depth/1000,'color',[.8 .2 .2],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 .05])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('P regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,2)
% plot(siprof,ao.depth/1000,'color',[.2 .8 .2],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 1])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('Si regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,3)
% plot(niprof,ao.depth/1000,'color',[.2 .2 .8],'linewidth',2); hold on;
% plot(cumni,ao.depth(4:24)/1000,'color',[.6 .6 .8],'linewidth',2); hold on;
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 .05])
% xlabel('(nM/yr)')
% ylabel('Depth (km)')
% title('Ni regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,4)
% plot(pprof_lolat,ao.depth/1000,'color',[.8 .2 .2],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 .05])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('P regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,5)
% plot(siprof_lolat,ao.depth/1000,'color',[.2 .8 .2],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 1])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('Si regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,6)
% plot(niprof_lolat,ao.depth/1000,'color',[.2 .2 .8],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 .05])
% xlabel('(nM/yr)')
% ylabel('Depth (km)')
% title('Ni regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,7)
% plot(pprof_midlat,ao.depth/1000,'color',[.8 .2 .2],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 .05])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('P regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,8)
% plot(siprof_midlat,ao.depth/1000,'color',[.2 .8 .2],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 1])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('Si regeneration flux','FontWeight','Normal')
% 
% subplot(3,3,9)
% plot(niprof_midlat,ao.depth/1000,'color',[.2 .2 .8],'linewidth',2)
% set(gca,'ylim',[0 3],'ydir','reverse','xlim',[0 .05])
% xlabel('(nM/yr)')
% ylabel('Depth (km)')
% title('Ni regeneration flux','FontWeight','Normal')
% 
% pdfname=['/Users/Seth/Desktop/temp/fig3alt.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 3 flux profiles
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(3); clf; papersize = [15 20];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% subplot(3,2,1)
% fill([pprof-2*err_pprof;flipud(pprof+2*err_pprof)],[ao.depth'/1000;flipud(ao.depth'/1000)],[1 .8 .8],'linestyle','none'); hold on;
% plot(pprof,ao.depth/1000,'color',[.8 .2 .2],'linewidth',1)
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 .04],'position',[.06 .75 .35 .2])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('P regeneration flux','FontWeight','Normal')
% 
% subplot(3,2,3)
% fill([siprof-2*err_siprof;flipud(siprof+2*err_siprof)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.8 1 .8],'linestyle','none'); hold on;
% plot(siprof,ao.depth/1000,'color',[.2 .8 .2],'linewidth',1)
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 .7],'position',[.06 .425 .35 .2])
% xlabel('(\muM/yr)')
% ylabel('Depth (km)')
% title('Si regeneration flux','FontWeight','Normal')
% 
% subplot(3,2,5)
% fill([niprof-2*err_niprof;flipud(niprof+2*err_niprof)],[ao.depth'/1000;flipud(ao.depth'/1000)],[.8 .8 1],'linestyle','none'); hold on;
% plot(niprof,ao.depth/1000,'color',[.2 .2 .8],'linewidth',1)
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 .05],'position',[.06 .1 .35 .2])
% xlabel('(nM/yr)')
% ylabel('Depth (km)')
% title('Ni regeneration flux','FontWeight','Normal')
% 
% subplot(3,2,2)
% fill([cump';zeros(21,1)],[ao.depth(4:24)'/1000;flipud(ao.depth(4:24)'/1000)],[.92 .9 .9],'linestyle','none'); hold on;
% plot(cump,ao.depth(4:24)/1000,'-','color',[.8 .2 .2],'linewidth',1); hold on;
% plot([0 1],[1 1],'--k','linewidth',.5); hold on; plot([0 1],[2 2],'--k','linewidth',.5); hold on; plot([0 1],[3 3],'--k','linewidth',.5); hold on; plot([0 1],[4 4],'--k','linewidth',.5); hold on; 
% text(0.05,.7,'77%','fontsize',8); text(0.05,1.7,'90%','fontsize',8); text(0.05,2.7,'93%','fontsize',8); text(0.05,3.7,'96%','fontsize',8); 
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 1],'xtick',(0:0.2:1),'position',[.56 .75 .35 .2])
% title('Cumulative P regeneration','FontWeight','Normal')
% 
% subplot(3,2,4)
% fill([cumsi';zeros(21,1)],[ao.depth(4:24)'/1000;flipud(ao.depth(4:24)'/1000)],[.9 .92 .9],'linestyle','none'); hold on;
% plot(cumsi,ao.depth(4:24)/1000,'-','color',[.2 .8 .2],'linewidth',1); hold on;
% plot([0 1],[1 1],'--k','linewidth',.5); hold on; plot([0 1],[2 2],'--k','linewidth',.5); hold on; plot([0 1],[3 3],'--k','linewidth',.5); hold on; plot([0 1],[4 4],'--k','linewidth',.5); hold on; 
% text(0.05,.7,'36 %','fontsize',8); text(0.05,1.7,'57%','fontsize',8); text(0.05,2.7,'68%','fontsize',8); text(0.05,3.7,'82%','fontsize',8); 
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 1],'xtick',(0:0.2:1),'position',[.56 .425 .35 .2])
% title('Cumulative Si regeneration','FontWeight','Normal')
% 
% subplot(3,2,6)
% fill([cumni';zeros(21,1)],[ao.depth(4:24)'/1000;flipud(ao.depth(4:24)'/1000)],[.9 .9 .95],'linestyle','none'); hold on;
% plot(cumni,ao.depth(4:24)/1000,'-','color',[.2 .2 .8],'linewidth',2); hold on;
% plot([0 1],[1 1],'--k','linewidth',.5); hold on; plot([0 1],[2 2],'--k','linewidth',.5); hold on; plot([0 1],[3 3],'--k','linewidth',.5); hold on; plot([0 1],[4 4],'--k','linewidth',.5); hold on; 
% text(0.05,.7,'58%','fontsize',8); text(0.05,1.7,'79%','fontsize',8); text(0.05,2.7,'86%','fontsize',8); text(0.05,3.7,'93%','fontsize',8); 
% set(gca,'ylim',[0 5],'ydir','reverse','xlim',[0 1],'xtick',(0:0.2:1),'position',[.56 .1 .35 .2])
% xlabel('Fraction regenerated')
% title('Cumulative Ni regeneration','FontWeight','Normal')
% 
% 
% pdfname=['/Users/Seth/Desktop/temp/fig3.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 2 flux transects across the surface
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(2); clf; papersize = [25 15];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% subplot(3,2,1)
% fill([ao.lat';flipud(ao.lat')],[ptransect-2*err_ptrans;flipud(ptransect+2*err_ptrans)],[1 .9 .9],'linestyle','none'); hold on;
% plot(ao.lat,ptransect,'color',[.8 .2 .2],'linewidth',1); hold on;
% set(gca,'xlim',[-60 60],'ylim',[0 0.6],'position',[.05 .75 .35 .2])
% ylabel('P uptake (\muM/yr)')
% title('Surface ocean P uptake','FontWeight','Normal')
% 
% subplot(3,2,3)
% fill([ao.lat';flipud(ao.lat')],[sitransect-2*err_sitrans;flipud(sitransect+2*err_sitrans)],[.9 1 .9],'linestyle','none'); hold on;
% plot(ao.lat,trans_Si_merid,'color',[.2 .8 .2],'linewidth',1);
% set(gca,'xlim',[-60 60],'ylim',[0 25],'position',[.05 .42 .35 .2])
% ylabel('Si uptake (\muM/yr)')
% title('Surface ocean Si uptake','FontWeight','Normal')
% 
% subplot(3,2,5)
% fill([ao.lat';flipud(ao.lat')],[nitrans-2*err_nitrans;flipud(nitrans+2*err_nitrans)],[.9 .9 1],'linestyle','none'); hold on;
% plot(ao.lat,nitrans,'color',[.2 .2 .8],'linewidth',1)
% set(gca,'xlim',[-60 60],'ylim',[0 .8],'position',[.05 .1 .35 .2])
% xlabel('Latitude')
% ylabel('Ni uptake (nM/yr)')
% title('Surface ocean Ni uptake','FontWeight','Normal')
% 
% subplot(3,2,4)
% fill([ao.lat';flipud(ao.lat')],[niptrans-2*err_niptrans;flipud(niptrans+2*err_niptrans)],[.9 .9 .9],'linestyle','none'); hold on;
% plot(ao.lat,niptrans,'-k','linewidth',1); hold on;
% set(gca,'xlim',[-60 60],'ylim',[0 4],'position',[.55 .6 .3 .2])
% xlabel('Latitude')
% ylabel('(\muM/nM)')
% title('Ni uptake/P uptake','FontWeight','Normal')
% 
% subplot(3,2,6)
% fill([ao.lat';flipud(ao.lat')],[nisitrans-2*err_nisitrans;flipud(nisitrans+2*err_nisitrans)],[.9 .9 .9],'linestyle','none'); hold on;
% plot(ao.lat,nisitrans,'-k','color',[.2 .2 .2],'linewidth',1)
% set(gca,'xlim',[-60 60],'ylim',[0 .5],'position',[.55 .2 .3 .2])
% xlabel('Latitude')
% ylabel('(\muM/nM)')
% title('Ni uptake/Si uptake','FontWeight','Normal')
% 
% 
% pdfname=['/Users/Seth/Desktop/temp/fig2.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Figure 1 flux maps
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1); clf; papersize = [30 20];
% set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);
% 
% % map for contouring
% surfmap = ao.OCN(:,:,1);
% meridmap = squeeze((max(ao.OCN,[],2)))'; meridmap(24,:)=meridmap(23,:)*0.5;
% 
% 
% subplot(3,2,1);
% pcolor(ao.lon,ao.lat,J_P_surf); shading flat; hold on;
% contour(ao.lon,ao.lat,surfmap,[1],'-k','linewidth',1);
% set(gca,'clim',[-1 .125],'color',[0.6 0.6 0.6],'ylim',[-80 70]);
% colormap(newmap)
% colorbar;
% ylabel('Latitude')
% title('Surface ocean P flux (\muM/yr)','FontWeight','Normal')
% 
% subplot(3,2,2);
% pcolor(ao.lat,ao.depth/1000,J_P_merid'); shading flat; hold on;
% contour(ao.lat,ao.depth/1000,meridmap,[1],'-k','linewidth',1);
% set(gca,'ydir','reverse','clim',[-1 .125],'color',[0.6 0.6 0.6],'xlim',[-80 70],'ylim',[0 5.6],'ytick',(0:1:5));
% colormap(newmap)
% colorbar;
% ylabel('Depth (km)')
% title('Meridional mean P flux (\muM/yr)','FontWeight','Normal')
% 
% subplot(3,2,3);
% pcolor(ao.lon,ao.lat,J_Si_surf); shading flat; hold on;
% contour(ao.lon,ao.lat,ao.OCN(:,:,1),'-k','linewidth',.3);
% set(gca,'clim',[-50 6.25],'color',[0.6 0.6 0.6],'ylim',[-80 70]);
% colormap(newmap)
% colorbar;
% ylabel('Latitude')
% title('Surface ocean Si flux (\muM/yr)','FontWeight','Normal')
% 
% subplot(3,2,4);
% pcolor(ao.lat,ao.depth/1000,J_Si_merid'); shading flat; hold on;
% contour(ao.lat,ao.depth/1000,meridmap,[1],'-k','linewidth',1);
% set(gca,'ydir','reverse','clim',[-50 6.25],'color',[0.6 0.6 0.6],'xlim',[-80 70],'ylim',[0 5.6],'ytick',(0:1:5));
% colormap(newmap)
% colorbar;
% ylabel('Depth (m)')
% title('Meridional mean Si flux (\muM/yr)','FontWeight','Normal')
% 
% subplot(3,2,5);
% pcolor(ao.lon,ao.lat,J_Ni_surf); shading flat; hold on;
% contour(ao.lon,ao.lat,ao.OCN(:,:,1),'-k','linewidth',.3);
% set(gca,'clim',[-2 .25],'color',[0.6 0.6 0.6],'ylim',[-80 70]);
% colormap(newmap)
% colorbar;
% xlabel('Longitude')
% ylabel('Latitude')
% title('Surface ocean Ni flux (nM/yr)','FontWeight','Normal')
% 
% subplot(3,2,6);
% pcolor(ao.lat,ao.depth/1000,J_Ni_merid'); shading flat; hold on;
% contour(ao.lat,ao.depth/1000,meridmap,[1],'-k','linewidth',1);
% set(gca,'ydir','reverse','clim',[-2 .25],'color',[0.6 0.6 0.6],'xlim',[-80 70],'ylim',[0 5.6],'ytick',(0:1:5));
% colormap(newmap)
% colorbar;
% xlabel('Latitude')
% ylabel('Depth (m)')
% title('Meridional mean Ni flux (nM/yr)','FontWeight','Normal')
% 
% pdfname=['/Users/Seth/Desktop/temp/fig1.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)


