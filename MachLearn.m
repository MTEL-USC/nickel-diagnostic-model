% load ../seth_Ni/GTGP15TARANi.mat
load gridding/GT_IDP21_Ni.mat
load ../data/WOA09/WOANO3.mat
load ../data/WOA09/WOAO2.mat
load ../data/WOA09/WOAPO4.mat
load ../data/WOA09/WOASAL.mat
load ../data/WOA09/WOASI.mat
load ../data/WOA09/WOATEMP.mat
load ../data/ao.mat
load Niclimatology.mat

% niobs = GTGP15TARANi(ao.iocn);
niobs = GT_IDP21_Ni(ao.iocn);
no3obs = WOANO3(ao.iocn);
o2obs = WOAO2(ao.iocn);
pobs = WOAPO4(ao.iocn);
salobs = WOASAL(ao.iocn);
siobs = WOASI(ao.iocn);
tempobs = WOATEMP(ao.iocn);
depth = ao.DEPTH(ao.iocn);
SINLAT = sin(ao.LAT/180*pi()); sinlat = SINLAT(ao.iocn);
COSLAT = cos(ao.LAT/180*pi()); coslat = COSLAT(ao.iocn);

% NaN out the data from the Arctic
niobs(ao.Lat>75)=NaN;

rerun=1;
if rerun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict the distributions using a machine learning technique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nruns = 100;

runANN = 0;
runTREE = 0;
runMLR = 1;

allclimatologies=zeros(200160,nruns)*NaN;
allR2=zeros(nruns,1)*NaN;
allR2in=zeros(nruns,1)*NaN;
allR2out=zeros(nruns,1)*NaN;
allSLOPE=zeros(nruns,1)*NaN;
allSLOPEin=zeros(nruns,1)*NaN;
allSLOPEout=zeros(nruns,1)*NaN;

if runMLR
    allts=zeros(nruns,10)*NaN;
end
if runTREE
    allimps=zeros(nruns,9)*NaN;
end

tic
for i=1:nruns

    %hold out a tenth of the data
    itenth=round(rand(200160,1)*(5/9));
    tenthobs = niobs; tenthobs(find(~itenth))=NaN;

    if runANN
    model = fitrnet([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],tenthobs,"Standardize",true);
    climatology = predict(model,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
    elseif runTREE
    model = fitrtree([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],tenthobs);
    climatology = predict(model,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
    imp = predictorImportance(model);
    allimps(i,:)=imp;
    elseif runMLR
    model = fitlm([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],tenthobs);
    climatology = predict(model,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
    t = abs(model.Coefficients.tStat);
    allts(i,:)=t;
    end

    allclimatologies(:,i)=climatology;

    jall=find(~isnan(niobs));
    R=corrcoef(niobs(jall),climatology(jall)); R2 = R(2)^2;
    allR2(i) = R2;
    SLOPE=niobs(jall)\climatology(jall);
    allSLOPE(i) = SLOPE;

    jin=find(~isnan(tenthobs));
    Rin=corrcoef(niobs(jin),climatology(jin)); R2in = Rin(2)^2;
    allR2in(i) = R2in;
    SLOPEin=niobs(jin)\climatology(jin);
    allSLOPEin(i) = SLOPE;

    kout=~ismember(jall,jin);
    jout=jall(kout);
    Rout=corrcoef(niobs(jout),climatology(jout)); R2out = Rout(2)^2;
    allR2out(i) = R2out;
    SLOPEout=niobs(jout)\climatology(jout);
    allSLOPEout(i) = SLOPE;

end

if runANN
model_allin = fitrnet([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],niobs,"Standardize",true);
climatology_allin = predict(model,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
elseif runTREE
model_allin = fitrtree([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],niobs);
climatology_allin = predict(model,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
imp = predictorImportance(model_allin)
elseif runMLR
model_allin = fitlm([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],niobs);
climatology_allin = predict(model,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
model_allin
end

end
toc

aveSLOPE = mean(allSLOPE)
errSLOPE = std(allSLOPE)
aveSLOPEin = mean(allSLOPEin)
errSLOPEin = std(allSLOPEin)
aveSLOPEout = mean(allSLOPEout)
errSLOPEout = std(allSLOPEout)

aveR2 = mean(allR2)
errR2 = std(allR2)
aveR2in = mean(allR2in)
errR2in = std(allR2in)
aveR2out = mean(allR2out)
errR2out = std(allR2out)

% plot the runs
figure(1); clf; papersize = [10 20];
set(gcf,'InvertHardCopy','Off','Renderer','opengl', 'color', 'w','PaperUnits','centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [100 100 500 500]);

bins = 0:0.2:12;
data = hist3([allclimatologies(:,1),niobs],'edges',{bins,bins}); data = data/(max(max(data)))*100;
data=smoothdata(data,1,'gaussian',2);
data=smoothdata(data,2,'gaussian',2);

subplot(2,1,1)
for i=1:nruns;
    scatter(niobs,allclimatologies(:,i),.5,'ok','filled'); hold on;
end
set(gca,'xlim',[0 12],'ylim',[0 12],'box','on','position',[.15 .5 .65 .3])
xlabel('Observed [Ni] (nM)');
ylabel ('Predicted [Ni] (nM)');


subplot(2,1,2)
contourf(bins,bins,data,50,'linestyle','none'); hold on;
% scatter(netni,niobs,'.k')
line([0 12],[0 12],'color','black','linestyle','--'); hold on;
crameri('-oslo')
% contour(bins,bins,data,(0:10:100),'-k'); hold on;
cb = colorbar;
% text(4.5,13, ['\it R^2 = ' num2str(r2,2)],'fontweight','bold','fontsize',14);
set(gca,'xlim',[0 12],'ylim',[0 12],'clim',[0 100],'position',[.15 .1 .65 .3]);
xlabel('Observed [Ni] (nM)');
ylabel ('Predicted [Ni] (nM)');

pdfname=['/Users/Seth/Desktop/temp/fig1.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)

Niclimatology = climatology_allin;

save('Niclimatology.mat','Niclimatology')






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % predict the distributions using Multiple Linear Regression
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nruns = 0;
% allMLRs=zeros(200160,nruns)*NaN;
% allR2=zeros(nruns,1)*NaN;
% allR2in=zeros(nruns,1)*NaN;
% allR2out=zeros(nruns,1)*NaN;
% 
% tic
% for i=1:nruns
% 
%     %hold out tenth the data
%     itenth=round(rand(200160,1)*(5/9));
%     tenthobs = niobs; tenthobs(find(~itenth))=NaN;
% 
%     MLR = fitlm([no3obs o2obs pobs salobs siobs tempobs depth sinlat coslat],tenthobs);
%     MLRni = predict(MLR,[no3obs o2obs pobs salobs siobs tempobs depth sinlat coslat]);
%     allMLRs(:,i)=MLRni;
% 
%     jni=find(~isnan(niobs));
%     R=corrcoef(niobs(jni),MLRni(jni));
%     R2 = R(2)^2;
%     allR2(i) = R2;
% 
%     jin=find(~isnan(tenthobs));
%     Rin=corrcoef(niobs(jin),MLRni(jin));
%     R2in = Rin(2)^2;
%     allR2in(i) = R2in;
% 
%     ktenthout=~ismember(jni,jin);
%     jtenthout=jni(ktenthout);
%     Rout=corrcoef(niobs(jtenthout),MLRni(jtenthout));
%     R2out = Rout(2)^2;
%     allR2out(i) = R2out;
% 
% end
% toc
% 
% % plot the runs
% figure(1); clf;
% for i=1:nruns;
%     scatter(niobs,allMLRs(:,i)); hold on;
% end
% 
% % remove failed attempts
% ifail = find(allR2<0.8);
% allMLRs(:,ifail)=[];
% allR2(ifail)=[];
% 
% % plot the runs
% figure(2); clf;
% for i=1:nruns-length(ifail);
%     scatter(niobs,allMLRs(:,i)); hold on;
% end
% 
% save('allMLRs.mat','allMLRs')
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % predict the distributions using regression trees
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nruns = 10;
% allTREEs=zeros(200160,nruns)*NaN;
% allR2=zeros(nruns,1)*NaN;
% allR2in=zeros(nruns,1)*NaN;
% allR2out=zeros(nruns,1)*NaN;
% 
% tic
% for i=1:nruns
% 
%     %hold out tenth the data
%     itenth=round(rand(200160,1)*(5/9));
%     tenthobs = niobs; tenthobs(find(~itenth))=NaN;
% 
%     TREE = fitrtree([no3obs o2obs pobs salobs siobs tempobs depth sinlat coslat],tenthobs);
%     TREEni = predict(TREE,[no3obs o2obs pobs salobs siobs tempobs depth sinlat coslat]);
%     allTREEs(:,i)=TREEni;
% 
%     jni=find(~isnan(niobs));
%     R=corrcoef(niobs(jni),TREEni(jni));
%     R2 = R(2)^2;
%     allR2(i) = R2;
% 
%     jin=find(~isnan(tenthobs));
%     Rin=corrcoef(niobs(jin),TREEni(jin));
%     R2in = Rin(2)^2;
%     allR2in(i) = R2in;
% 
%     ktenthout=~ismember(jni,jin);
%     jtenthout=jni(ktenthout);
%     Rout=corrcoef(niobs(jtenthout),TREEni(jtenthout));
%     R2out = Rout(2)^2;
%     allR2out(i) = R2out;
% 
% end
% toc
% 
% % plot the runs
% figure(1); clf;
% for i=1:nruns;
%     scatter(niobs,allTREEs(:,i)); hold on;
% end
% 
% % remove failed attempts
% ifail = find(allR2<0.8);
% allMLRs(:,ifail)=[];
% allR2(ifail)=[];
% 
% % plot the runs
% figure(2); clf;
% for i=1:nruns-length(ifail);
%     scatter(niobs,allTREEs(:,i)); hold on;
% end
% 
% save('allMLRs.mat','allMLRs')
% 
% 






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % predict the distributions using ANNs
% % old method using random numbers of layers
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % general approach, 2 to 4 hidden layers, 10, 20, 30 nodes per
% % layer, withhold 50% of data each time
% nruns = 0;
% allnets=zeros(200160,nruns)*NaN;
% allR2=zeros(nruns,1)*NaN;
% allR2in=zeros(nruns,1)*NaN;
% allR2out=zeros(nruns,1)*NaN;
% 
% tic
% for i=1:nruns
% 
%     %hold out tenth the data
%     itenth=round(rand(200160,1)*(5/9));
%     tenthobs = niobs; tenthobs(find(~itenth))=NaN;
% 
%     nlayers = randi([2 4]);
%     nodes = [10 20 30];
%     nnodes = [nodes(randi([1 3])) nodes(randi([1 3])) nodes(randi([1 3])) nodes(randi([1 3]))];
%     if nlayers==2;
%         nnet = fitrnet([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],tenthobs,"LayerSizes",[nnodes(1) nnodes(2)],"Standardize",true);
%     elseif nlayers==3;
%         nnet = fitrnet([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],tenthobs,"LayerSizes",[nnodes(1) nnodes(2) nnodes(3)],"Standardize",true);
%     else nlayers==4;
%         nnet = fitrnet([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],tenthobs,"LayerSizes",[nnodes(1) nnodes(2) nnodes(3) nnodes(4)],"Standardize",true);
%     end
%     
%     nnetni = predict(nnet,[no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat]);
%     allnets(:,i)=nnetni;
% 
%     jni=find(~isnan(niobs));
%     R=corrcoef(niobs(jni),nnetni(jni));
%     R2 = R(2)^2;
%     allR2(i) = R2;
% 
%     jtenth=find(~isnan(tenthobs));
%     Rin=corrcoef(niobs(jtenth),nnetni(jtenth));
%     R2in = Rin(2)^2;
%     allR2in(i) = R2in;
% 
%     ktenthout=~ismember(jni,jtenth);
%     jtenthout=jni(ktenthout);
%     Rout=corrcoef(niobs(jtenthout),nnetni(jtenthout));
%     R2out = Rout(2)^2;
%     allR2out(i) = R2out;
% 
% end
% toc
% 
% allR2;
% 
% % plot the runs
% figure(1); clf;
% for i=1:nruns;
%     scatter(niobs,allnets(:,i)); hold on;
% end
% 
% % remove failed attempts
% ifail = find(allR2<0.8);
% allnets(:,ifail)=[];
% allR2(ifail)=[];
% 
% % plot the runs
% figure(2); clf;
% for i=1:nruns-length(ifail);
%     scatter(niobs,allnets(:,i)); hold on;
% end
% 
% save('allnets.mat','allnets')

% test= fitrnet([no3obs o2obs pobs salobs siobs tempobs depth,sinlat,coslat],niobs);
% test2=predict(test,[no3obs o2obs pobs salobs siobs tempobs depth]);

% tree = fitrtree([no3obs o2obs pobs salobs siobs tempobs depth],niobs);
% 
% % imp = predictorImportance(tree)
% % view(tree,'Mode','graph')
% 
% MLni = predict(tree,[no3obs o2obs pobs salobs siobs tempobs depth]);
% 
% i=find(~isnan(niobs));
% R = corr(niobs(i),MLni(i));
% R2 = R^2;
% 
% MLNI = ao.OCN;
% MLNI(:) = MLni;
% MLNI = MLNI.*ao.nanOCN;
% 
% % save ('MLNI.mat','MLNI')
% % 
% % figure(1); clf; papersize = [15 15];
% % set(gcf,'Renderer','opengl', 'color', 'w','PaperUnits', 'centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [0 0 500 300]);
% % 
% % scatter(ni,MLni,2,'k')
% % xlabel('GEOTRACES Ni')
% % ylabel('predicted Ni')
% 
% % pdfname=['/Users/Seth/Desktop/temp/fig.pdf' ]; print('-dpdf','-r300',pdfname); open (pdfname)
% 
% 
% mlni = MLNI(ao.iocn);
% woap = WOAPO4(ao.iocn);
% woasi = WOASI(ao.iocn);
% 
% % restoring timescales of 5 months in the surface ocean (top 3 layers) and
% % 30 months in the deep ocean, after Roshan and DeVries
% TAU = ao.OCN; TAU(:,:,1:3)=.5; TAU(:,:,4:24)=25; TAU=TAU.*ao.OCN;
% % TAU=ones(size(ao.OCN));
% tau=TAU(ao.iocn);
% tau=tau*2;
% Atau = speye(200160); Atau(Atau==1)=(1./tau);
% 
% % fprintf('\nSetting to work on po4\n')
% % p=(water_transport-Atau)\(-woap./tau);
% % fprintf('\nSetting to work on Si\n')
% % si=(water_transport-Atau)\(-woasi./tau);
% % fprintf('\nSetting to work on Ni\n');
% % ni=(water_transport-Atau)\(-mlni./tau);
% 
% P=ao.nanOCN; P(P==1)=p;
% SI=ao.nanOCN; SI(SI==1)=si;
% NI=ao.nanOCN; NI(NI==1)=ni;
% 
% jp = -water_transport*p;
% jsi = -water_transport*si;
% jni = -water_transport*ni;
% 
% % jp = (woap-p)./tau;
% % jsi = (woasi-si)./tau;
% % jni = (mlni-ni)./tau;
% 
% NOARCTIC = ao.nanOCN; NOARCTIC(80:90,:,:)=NaN;
% 
% JP=ao.nanOCN; JP(JP==1)=jp; JP=JP.*NOARCTIC;
% JSI=ao.nanOCN; JSI(JSI==1)=jsi; JSI=JSI.*NOARCTIC;
% JNI=ao.nanOCN; JNI(JNI==1)=jni; JNI=JNI.*NOARCTIC;
% % JP=JP.*ao.PAC;
% % JSI=JSI.*ao.PAC;
% % JNI=JNI.*ao.PAC;
% 
% UP=-JP; UP(UP<0)=0;
% USI=-JSI; USI(USI<0)=0;
% UNI=-JNI; UNI(UNI<0)=0;
% 
% RP=JP; RP(RP<0)=0;
% RSI=JSI; RSI(RSI<0)=0;
% RNI=JNI; RNI(RNI<0)=0;
% 
% pprof = squeeze(nansum(squeeze(nansum(RP,1)),1))./squeeze(nansum(squeeze(nansum(ao.VOL,1)),1));
% siprof = squeeze(nansum(squeeze(nansum(RSI,1)),1))./squeeze(nansum(squeeze(nansum(ao.VOL,1)),1));
% niprof = squeeze(nansum(squeeze(nansum(RNI,1)),1))./squeeze(nansum(squeeze(nansum(ao.VOL,1)),1));
% 
% pprof(1:4)=NaN; siprof(1:4)=NaN; niprof(1:4)=NaN;
% 
% plat = nansum(nansum(UP,3),2)./nansum(nansum(ao.VOL(:,:,1:1),3),2);
% silat = nansum(nansum(USI,3),2)./nansum(nansum(ao.VOL(:,:,1:1),3),2);
% nilat = nansum(nansum(UNI,3),2)./nansum(nansum(ao.VOL(:,:,1:1),3),2);
% 
% figure(2); clf; papersize = [12 5];
% set(gcf,'Renderer','opengl', 'color', 'w','PaperUnits', 'centimeters','PaperSize', [papersize(1) papersize(2)],'PaperPosition', [1 1 papersize(1)-1 papersize(2)-1],'position', [0 0 500 1000]);
% subplot(4,3,1)
% plot(pprof,ao.depth,'g','linewidth',2)
% set(gca,'ydir','reverse','xlim',[0 inf])
% title('P regeneration')
% 
% subplot(4,3,2)
% plot(siprof,ao.depth,'b','linewidth',2)
% set(gca,'ydir','reverse','xlim',[0 inf])
% title('Si regeneration')
% 
% subplot(4,3,3)
% plot(niprof,ao.depth,'r','linewidth',2)
% set(gca,'ydir','reverse','xlim',[0 inf])
% title('Ni regeneration')
% 
% subplot(4,3,4)
% plot(pprof./pprof,ao.depth,'g','linewidth',2)
% set(gca,'ydir','reverse')
% title('P regen/P regen')
% 
% 
% subplot(4,3,5)
% plot((siprof./pprof)*(nansum(pprof)/nansum(siprof)),ao.depth,'b','linewidth',2)
% set(gca,'ydir','reverse')
% title('Si regen/P regen')
% 
% subplot(4,3,6)
% plot(niprof./pprof*(nansum(pprof)/nansum(niprof)),ao.depth,'r','linewidth',2)
% set(gca,'ydir','reverse')
% title('Ni regen/P regen')
% 
% subplot(4,3,7)
% plot(ao.lat,plat,'g','linewidth',2)
% set(gca,'xlim',[-80 70],'ylim',[0 inf])
% title('P uptake')
% 
% 
% subplot(4,3,8)
% plot(ao.lat,silat,'b','linewidth',2)
% set(gca,'xlim',[-80 70],'ylim',[0 inf])
% title('Si uptake')
% 
% subplot(4,3,9)
% plot(ao.lat,nilat,'r','linewidth',2)
% set(gca,'xlim',[-80 70],'ylim',[0 inf])
% title('Ni uptake')
% 
% subplot(4,3,10)
% plot(ao.lat,plat./plat,'g','linewidth',2)
% set(gca,'xlim',[-80 70])
% title('P uptake/P uptake')
% 
% subplot(4,3,11)
% plot(ao.lat,silat./plat*(nansum(plat)/nansum(silat)),'b','linewidth',2)
% set(gca,'xlim',[-80 70])
% title('Si uptake/P uptake')
% 
% subplot(4,3,12)
% plot(ao.lat,nilat./plat*(nansum(plat)/nansum(nilat)),'r','linewidth',2)
% set(gca,'xlim',[-80 70])
% title('Ni uptake/P uptake')
% 
% 
% % figure(1); clf;
% % 
% % subplot(2,3,1)
% % pcolor(ao.lat,ao.depth,squeeze(WOAPO4(:,100,:))');
% % set(gca,'ydir','reverse','clim',[0 3]); colorbar;
% % title('WOA P at 160W')
% % 
% % subplot(2,3,2)
% % pcolor(ao.lat,ao.depth,squeeze(WOASI(:,100,:))');
% % set(gca,'ydir','reverse','clim',[0 150]); colorbar;
% % title('WOA SI at 160W')
% % 
% % subplot(2,3,3)
% % pcolor(ao.lat,ao.depth,squeeze(MLNI(:,100,:))');
% % set(gca,'ydir','reverse','clim',[0 10]); colorbar;
% % title('Machine learned Ni at 160W')
% % 
% % subplot(2,3,4)
% % pcolor(ao.lat,ao.depth,squeeze(P(:,100,:))');
% % set(gca,'ydir','reverse','clim',[0 3]); colorbar;
% % title('P at 160W smoothed by restoring')
% % 
% % subplot(2,3,5)
% % pcolor(ao.lat,ao.depth,squeeze(SI(:,100,:))');
% % set(gca,'ydir','reverse','clim',[0 150]); colorbar;
% % title('Si at 160W smoothed')
% % 
% % subplot(2,3,6)
% % pcolor(ao.lat,ao.depth,squeeze(NI(:,100,:))');
% % set(gca,'ydir','reverse','clim',[0 10]); colorbar;
% % title('Ni at 160W smoothed')
