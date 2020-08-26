%% export data for analyses in R 


%% export PEAKs FCz 

P2PFCZ = p2p(Time_FRN, 200, 250, FCz, 50);
PATH = '/Volumes/daten/romy/Confidence/';
save(sprintf('%sExport/P2PFCZ.mat', PATH), 'P2PFCZ', '-v6');



%% export mean P3a

P3=nanmean(Time_FRN(:,265:315,:),2);

P3=reshape(P3, [65, 10000]);

P3=P3';

PATH = '/Volumes/daten-2/romy/Confidence/';
save(sprintf('%sExport/P3.mat', PATH), 'P3', '-v6');

%% export mean P3b

P3b=nanmean(Time_FRN(:,308: 358,:),2);

P3b=reshape(P3b, [65, 10000]);

P3b=P3b';

PATH = '/Volumes/daten-2/romy/Confidence/';
save(sprintf('%sExport/P3b.mat', PATH), 'P3b', '-v6');
