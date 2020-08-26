%% plotting ERP data

addpath(genpath('/Users/Romy/Beruflich/EEGfunctions'))
PATH='/Users/Romy/Dropbox (Brown)/Confidence/';
%% get information on single trials (file made in R from logfiles)

load(sprintf('%sExport/behav_n.mat',PATH), 'a1')
%%
load(sprintf('%sExport/plotes.mat', PATH), 'Emat')

load(sprintf('%sExport/plotts.mat', PATH), 'Tmat')

load(sprintf('%sExport/plotes2.mat', PATH), 'Emat2')

load(sprintf('%sExport/plotts2.mat', PATH), 'Tmat2')
% a1 = 
% 
%              vpn: {8250x1 cell}
%            Trial: [8250x1 double]
%               RT: [8250x1 double]
%       absE: [8250x1 double]     % deviation from target interval
%      sabsE: [8250x1 double]     % Feedback
%        Direction: [8250x1 double]     % estimated deviation
%       Confidence: [8250x1 double]     % confidence rating
%            block: [8250x1 double]
%          block25: [8250x1 double]
%             absE: [8250x1 double]     % absolute error (unsigned, too
%                                       % slow and too fast reactions collapsed)
%              acE: [8250x1 double]     % Feedback - estimated deviation
%     asabsE: [8250x1 double]     % absolute feedback (no direction)
%       aDirection: [8250x1 double]     % absolute estimated deviation
%              raw: [8500x1 double]
%             fraw: {8500x1 cell}
%% get unique subject names for later processing and preselection of partially existing data
 vps=unique(a1.vpn);

 %% make inset for FRN figure
Rs= [linspace(0,153/255, 50) linspace(153/255,0.8, 50) linspace(0.8,153/255, 50) linspace(153/255,0, 50)];
Gs= [linspace(0,153/255, 50) linspace(153/255,0, 50) linspace(0,153/255, 50) linspace(153/255,0, 50)];
Bs= [linspace(0.8,153/255, 50) linspace(153/255,0, 50) linspace(0,153/255, 50) linspace(153/255,0.8, 50)];
bgrgb= [Rs(:), Gs(:), Bs(:)];

x1 = [-1:.01:1];
norm1 = normpdf(x1,0,0.1);

norm2 = normpdf(x1,-0.7,0.1);

norm3 = normpdf(x1,0.3,0.15);
 figure()
 set(gcf,'Color', [1 1 1])

 subplot(2,1,1)
  hold on
f0 = fill([-1 -1 0.2 0.2],[-2 7 7 -2],[0.8 0 0]);
cdata=get(f0,'xdata');
cdata=(cdata-min(cdata))/(max(cdata)-min(cdata)); 
set(f0,'CData',cdata,'FaceColor','interp')
f0.EdgeColor = 'none';
f0b = fill([-0.2 -0.2 1 1],[-2 7 7 -2],[0.8 0 0]);
cdata=get(f0b,'xdata');
cdata=(cdata-min(cdata))/(max(cdata)-min(cdata)); 
set(f0b,'CData',cdata,'FaceColor','interp')
f0b.EdgeColor = 'none';
f0c = fill([-0.4 -0.4 0.4 0.4],[-2 7 7 -2],[0.8 0 0]);
f0c.EdgeColor = 'none';
colormap(bgrgb)
 xlabel('Deviation from Goal','fontsize', 14)
 ylabel('a.u.','fontsize', 14)
 plot(x1,norm2, 'k','LineWidth',2)
 k1=arrow([-0.5, 5],[-0.5 0.2],'Width', 1,'Color',[153/255, 153/255, 153/255]);
 plot(zeros(length([-2:1:7])), [-2:1:7], 'k','Linestyle',':','LineWidth',2)
 t=gridxy(0, 'Linestyle',':','LineWidth',2);
 t=gridxy(-0.5, 'Linestyle',':','LineWidth',1);
 t=gridxy(0.5, 'Linestyle',':','LineWidth',1);
 ylim([0,5])
  xlim([-1,1])
 set(gca,'TickDir','out');
 subplot(2,1,2)
  hold on
 f1 = fill([-0.6 -0.6 0.6 0.6],[-2 7 7 -2],[0.8 0 0]);
cdata=[-0.7 -0.7 0.7 0.7];
set(f1,'CData',cdata,'FaceColor','interp')
f1.EdgeColor = 'none';
f2 = fill([-1 -1 -0.6 -0.6],[-2 7 7 -2],[0 0 0.8]);
f2.EdgeColor = 'none';
f3 = fill([0.6 0.6 1 1],[-2 7 7 -2],[0 0 0.8]);
f3.EdgeColor = 'none';
  k1=arrow([-0.5, 5],[-0.5 0.2],'Width', 1,'Color',[153/255, 153/255, 153/255]);
colormap(bgrgb)
 xlabel('Deviation from Goal','fontsize', 14)
 ylabel('a.u.','fontsize', 14)
 plot(x1,norm3, 'k','LineWidth',2)
 plot(zeros(length([-2:1:7])), [-2:1:7], 'k','Linestyle',':','LineWidth',2)
 t=gridxy(0, 'Linestyle',':','LineWidth',2);
 t=gridxy(-0.5, 'Linestyle',':','LineWidth',1);
 t=gridxy(0.5, 'Linestyle',':','LineWidth',1);
 ylim([0,5])
 set(gca,'TickDir','out');
%%
set(gcf,'units','centimeters','position',[0 0 12 12])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/FRN_Predictions','-pdf','-painters'); 

 %% make inset for P3a figure;
 
 x1 = [-1:.01:1];
norm1 = normpdf(x1,0,0.1);

norm2 = normpdf(x1,0,0.2);

norm3 = normpdf(x1,0,0.5);


mean(4-norm1)
mean(4-norm2)
mean(4-norm3)
trapz(x1,norm1)
trapz(x1,norm2)
trapz(x1,norm3)
 figure()
 set(gcf,'Color', [1 1 1])
 title('Internal representation of Prediction','fontsize', 14)
 hold on
 L1=plot(x1,norm1, 'linewidth',2, 'color',[86/255, 180/255, 233/255]);
 L2=plot(x1,norm2, 'linewidth',2, 'color',[230/255, 159/255, 0/255]);
 L3=plot(x1,norm3,'linewidth',2, 'color',[153/255, 153/255, 153/255]);
 
 ylim([0,4.5])
 set(gca,'TickDir','out');
 t=gridxy(0, 'Linestyle',':','LineWidth',1);
 xlabel('Deviation from Prediction','fontsize', 14)
 ylabel('a.u.','fontsize', 14)
 h =legend([L3, L2, L1], 'low','medium','high',  'Location', 'southoutside', 'Orientation','horizontal');
 h.Title.String = 'Confidence';
 lh=findall(gcf,'tag','legend');
 set(lh,  'fontsize', 12)
set(lh,'box','off')
 set(gcf,'units','centimeters','position',[0 0 12 12])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_Predictions','-pdf','-painters'); 


%% compute shannon information for values in x under different confidence distributions
cvals = [0.001:0.033:1 ];
SurpriseMat= zeros(length(cvals), length(x1));
figure
hold on
for ii = 1:length(cvals)
    
    pdist = normpdf(x1,0,cvals(ii));
    pdist = pdist/sum(pdist);
    plot(x1,pdist)
    for iii = 1:length(x1)
        if pdist(x1==x1(iii))< 1.2178e-05
            
            pdist(x1==x1(iii))=  1.2178e-05;
        end
        SurpriseMat(ii, iii)= -log(pdist(x1==x1(iii)));
    end
end




%%
figure()
 set(gcf,'Color', [1 1 1])
h = heatmap(SurpriseMat,'ColorLimits',[5 10],'GridVisible','off' )
%% compute mean surprise given OPE:

Confsurprise = mean(SurpriseMat, 2);
OPEsurprise = mean(SurpriseMat, 1);
%%
figure()
 set(gcf,'Color', [1 1 1])
 subplot(2,1,1)
h = heatmap(Confsurprise,'ColorLimits',[5 10],'GridVisible','off' )

 subplot(2,1,2)
h = heatmap(OPEsurprise,'ColorLimits',[5 10],'GridVisible','off' )
%%
Rs= [ linspace(153/255,0.8, 50) ];
Gs= [ linspace(153/255,0, 50) ];
Bs= [ linspace(153/255,0, 50) ];
gr= [Rs(:), Gs(:), Bs(:)];

S= ' ';
xsublab = cell(1,length(x1));
[xsublab{1,:}]= deal(S);
xsublab{1,1} = '-1';
xsublab{1,find(round(x1,4)==-0.5)}= '-0.5';
xsublab{1,find(round(x1,4)==0)}= '0';
xsublab{1,find(round(x1,4)==0.5)}= '0.5';
xsublab{1,end}= '1';
figure()
set(gcf,'Color', [1 1 1])
subplot(9,9, [2:9 9+2:2*9 2*9+2:3*9 3*9+2:4*9 4*9+2:5*9 5*9+2:6*9 6*9+2:7*9 7*9+2:8*9 ])
h = heatmap(SurpriseMat,'ColorLimits',[5 9],'GridVisible','off','Colormap',gr )
colormap(gr)
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];
subplot(9,9, [1 9+1 2*9+1 3*9+1 4*9+1 5*9+1 6*9+1 7*9+1])
h = heatmap(Confsurprise,'ColorLimits',[5 9],'GridVisible','off', 'YLabel','Confidence','FontSize',14  ,'ColorbarVisible','off','Colormap',gr ,'CellLabelColor', 'none')
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];
subplot(9, 9, [74:81])
h = heatmap(OPEsurprise,'ColorLimits',[5 9],'GridVisible','off','ColorbarVisible','off','XLabel','Deviation from Predition','FontSize',14 ,'Colormap',gr )
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];

%%
figure()
set(gcf,'Color', [1 1 1])
subplot(9,9, [2:9 9+2:2*9 2*9+2:3*9 3*9+2:4*9 4*9+2:5*9 5*9+2:6*9 6*9+2:7*9 7*9+2:8*9 ])
h = heatmap(SurpriseMat(:, 101:end),'ColorLimits',[5 9],'GridVisible','off','Colormap',gr )
colormap(gr)
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];
subplot(9,9, [1 9+1 2*9+1 3*9+1 4*9+1 5*9+1 6*9+1 7*9+1])
h = heatmap(Confsurprise,'ColorLimits',[5 9],'GridVisible','off', 'YLabel','Marginal Confidence','FontSize',14  ,'ColorbarVisible','off','Colormap',gr ,'CellLabelColor', 'none')
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];
subplot(9, 9, [74:81])
h = heatmap(OPEsurprise(:, 101:end),'ColorLimits',[5 9],'GridVisible','off','ColorbarVisible','off','XLabel','Marginal OPE','FontSize',14 ,'Colormap',gr )
old_warning_state = warning('off', 'MATLAB:structOnObject');
hs = struct(h);
warning(old_warning_state);
hs.XAxis.TickValues = [];
hs.YAxis.TickValues = [];

export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_Predictions_Shannon_w_Marginals','-pdf','-painters'); 

%% load EEG data

load(sprintf('%sExport/Time_FRN.mat', PATH), 'Time_FRN');

%% Define electrode numbers:
FP1=1;  FPz=2;  FP2=3;  AF3=4;  AFz=5;  AF4=6;   F7=7;   F5=8;   F3=9;   F1=10; 
Fz=11;  F2=12;  F4=13;  F6=14;  F8=15;
FT9=16;  FT7=17; FC5=18; FC3=19; FC1=20;
FCz=21; FC2=22; FC4=23; FC6=24; FT8=25; FT10=26; T7=27;  C5=28;  C3=29;  C1=30; 
C2=31;  C4=32;  C6=33;  T8=34;  TP9=35; TP7=36;  CP5=37; CP3=38; CP1=39; CPz=40;
CP2=41; CP4=42; CP6=43; TP8=44; TP10=45;P7=46;   P5=47;  P3=48;  P1=49;  Pz=50;
P2=51;  P4=52;  P6=53;  P8=54;  PO3=55; POz=56;  PO4=57; O1=58;  Oz=59;  O2=60;
LO1=61; IO1=62; IO2=63; LO2=64; Cz=65; 

%%

load(sprintf('%sExport/chanlocs.mat', PATH), 'chanlocs');
%% 
TIME    = -200:2:798;
XMIN = -200;
XMAX =  800;

%% Plot effects

%% Confidence
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('Confidence effect')
topoplot(Emat.oConfidence,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat.oConfidence,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1}, 'conv', 'on' );set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_conf','-eps','-nocrop');%     

%% OPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('OPE effect')
topoplot(Emat.ssacE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat.ssacE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_OPE','-eps','-nocrop');
%% Block
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('Block effect')
topoplot(Emat.block,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat.block,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
%%
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('Error magnitude effect')
topoplot(Emat.sassDifference,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat.sassDifference,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_EM','-eps','-nocrop');
%% RPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('RPE effect')
topoplot(Emat.ssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat.ssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_RPE','-eps','-nocrop');
%%
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('Error magnitude effect')
topoplot(Emat2.sassDifference,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat2.sassDifference,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3b_EM','-eps','-nocrop');
%% RPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('RPE effect')
topoplot(Emat2.ssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat2.ssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3b_RPE','-eps','-nocrop');
%% OPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('OPE effect')
topoplot(Emat2.ssacE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat2.ssacE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3b_OPE','-eps','-nocrop');

%% Block
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('Block effect')
topoplot(Emat2.block,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat2.block,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
 set(gcf,'units','centimeters','position',[0 0 20 8])

%% blockbyoConfidencebyssPE

PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(1,2,1)
title('RPE by Confidence by block effect')
topoplot(Emat2.blockbyoConfidencebyssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');
subplot(1,2,2)
title('t-values')
MAPLIM = [-3 3];
topoplot(Tmat2.blockbyoConfidencebyssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1}, 'conv', 'on');set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
colorbar('EastOutside');

%% Big fucking figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Confidence
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(2,2,1)
title('Confidence', 'fontsize', 12)
topoplot(Emat.oConfidence,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

% OPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
subplot(2,2,2)
title('OPE', 'fontsize', 12)
topoplot(Emat.ssacE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

subplot(2,2,3)
title('RPE', 'fontsize', 12)
topoplot(Emat.ssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
hp4 = get(subplot(2,2,4),'Position')

colorbar('Eastoutside','Position', [hp4(1)+hp4(3)+hp4(3)/25  hp4(2)+ hp4(4)/6  hp4(3)/8  hp4(4)*2.5-hp4(4)/4]);

PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
subplot(2,2,4)
title('Error magnitude', 'fontsize', 12)
topoplot(Emat.sassDifference,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ F1, Fz, F2,FC1, FCz, FC2, C1, Cz, C2],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

set(gcf,'units','centimeters','position',[0 0 10 10])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3a_alleff_new','-pdf','-painters'); 

%% P3b 


% Confidence
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
figure();
subplot(2,2,4)
title(sprintf('Confidence\n x RPE x Block'), 'fontsize', 12)
topoplot(Emat2.blockbyoConfidencebyssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

% OPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
subplot(2,2,3)
title('OPE', 'fontsize', 12)
topoplot(Emat2.ssacE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

% RPE
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
subplot(2,2,1)
title('RPE', 'fontsize', 12)
topoplot(Emat2.ssPE,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])
hp4 = get(subplot(2,2,4),'Position')

colorbar('Eastoutside','Position', [hp4(1)+hp4(3)+hp4(3)/25  hp4(2)+ hp4(4)/6  hp4(3)/8  hp4(4)*2.5-hp4(4)/4]);
%
PLOTRAD = .70;
HEADRAD = .70;
INTRAD  =  0.9;
MAPLIM = [-5 5];
subplot(2,2,2)
title('Error magnitude', 'fontsize', 12)
topoplot(Emat2.sassDifference,chanlocs,'electrodes','off','maplimits',MAPLIM,'plotrad',PLOTRAD,'headrad',HEADRAD,'intrad',INTRAD, 'emarker2', {[ CP1, CPz, CP2, P1, Pz, P2, PO3, POz, PO4],'o','k', 3,1});set(gca,'Color','none');
xlim([-0.6 0.6])
ylim([-0.6 0.6])

 set(gcf,'units','centimeters','position',[0 0 10 10])
 export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/P3b_alleff_new','-pdf','-painters'); 

%% Peaks P3a, P3b:


GrM = nanmean(Time_FRN(:,:,:),3);

%%
YMIN=-2;
YMAX=7;
figure
subplot(3,1,1)
hold on;
f0 = fill([150 150 300 300],[-2 7 7 -2], [.96 .96 .96]);
f0.EdgeColor = 'none';
f0.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
title('FCz')
set(gcf,'Color', [1 1 1])
plot(TIME,GrM(FCz,:),'k','linewidth',1.5)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabh = get(gca,'XLabel');
ylabel('amplitude [µV]');%, 'fontsize', 14
set(xlabh,'Position',get(xlabh,'Position') - [-700 0 0]);
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
set(gca,'TickDir','out');
gridxy(0, 'Linestyle',':','LineWidth',1);
subplot(3,1,2)
hold on;
f1 = fill([330 330 430 430],[-2 7 7 -2], [.89 .89 .89]);
f1.EdgeColor = 'none';
f1.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
title('FCz')
set(gcf,'Color', [1 1 1])
plot(TIME,GrM(FCz,:),'k','linewidth',1.5)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabh = get(gca,'XLabel');
ylabel('amplitude [µV]');%, 'fontsize', 14
set(xlabh,'Position',get(xlabh,'Position') - [-700 0 0]);
set(gca,'TickDir','out');
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
subplot(3,1,3)
set(gcf,'Color', [1 1 1])
title('Pz')
set(gcf,'Color', [1 1 1])
hold on;
f2 = fill([416 416 516 516],[-2 7 7 -2], [.75 .75 .75]);
f2.EdgeColor = 'none';
f2.EdgeColor = 'none';
set(gcf,'Color', [1 1 1])
A1=plot(TIME,GrM(Pz,:),'k','linewidth',1.5)
xlim([XMIN XMAX])
ylim([YMIN YMAX])
xlabh = get(gca,'XLabel');
ylabel('amplitude [µV]', 'fontsize', 14);
set(gca,'TickDir','out');
xlabel('time [ms]', 'fontsize', 14);
%set(xlabh,'Position',get(xlabh,'Position') - [-600 -0.8 0]);
h = legend([f0, f1, f2],'FRN time window','P3a time window', 'P3b time window'); 
     lh=findall(gcf,'tag','legend');
     lp=get(lh,'position')
     set(lh,'position',[.50,.0,lp(3:4)])%,'fontsize',16
     set(lh,'box','off')
hline = refline([0 0]);
set(hline,'Color','black','Linestyle',':','LineWidth',1);
gridxy(0, 'Linestyle',':','LineWidth',1);
set(gcf,'units','centimeters','position',[0 0 8 40])
export_fig('/Users/Romy/Beruflich/Arbeit/Freiberufliche Projekte/Seminar_Confidence/Plots/FRN_P3a_b_time_intervals','-pdf','-painters'); 


