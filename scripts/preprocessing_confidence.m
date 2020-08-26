

%% preparation
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

PATH = '/Volumes/daten/romy/'; %edit to your needs 

SOURCEFILES = dir(strcat(PATH, 'Confidence/raw/*.vhdr')); %all brain vision header files in the folder
SUBJECTS = 1:numel(SOURCEFILES);
%%
for s =1:numel(SOURCEFILES)%SUBJECTS
%% load Data
    
    fn   = SOURCEFILES(s).name(1:2);
    s_id = str2double(SOURCEFILES(s).name(1:2)) % reads in subject ID and converts it to double (does not work for '01')
    if isnan(s_id) % if operation returns NaN (because of zero)
    s_id = str2double(SOURCEFILES(s).name(2)) % then convert only second part of the string (2-9 in our Data)  
    end
    
    
    %%
    EEG = pop_loadbv(sprintf('%sConfidence/raw/',PATH),sprintf('%s',SOURCEFILES(s).name));


%% do OC 

%add empty channel Cz
    EEG.data(end+1,:) = 0;
    EEG.nbchan = size(EEG.data,1);
    EEG.chanlocs(end+1).labels = 'Cz';
%%
  
    %load chanloc information
    %EEG = pop_chanedit(EEG, 'lookup',sprintf('%sOccular_Correction/werfen_BESA.elp',PATH));    
    EEG=pop_chanedit(EEG, 'lookup',sprintf('%seeglab13_5_4b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp', PATH));
    %re-reference to average
    EEG = pop_reref( EEG, []);

%%

 msec_matrix = dlmread(sprintf('%sConfidence/OC/%s_Cali.matrix', PATH, fn), '\t', 1, 1); % Namen einfügen
        
 EEG.data = msec_matrix*EEG.data;
 EEG = pop_select(EEG,'channel',1:65);
  
 %% Filtering 
 %data, lowcut, hicut
    %high-cut
    EEG = pop_eegfiltnew(EEG, [], 40);
    %low-cut
    EEG = pop_eegfiltnew(EEG, 0.5, []); 
    
    
 %%   save clean data
  
    %save set
    pop_saveset(EEG, 'filename',sprintf('%sConfidence/EEG_new/%s_clean',PATH,fn));

  
end