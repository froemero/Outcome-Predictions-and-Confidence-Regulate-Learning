%% File to extract single trials:

% note that this was written in 2015/2016. If you want a better and evolved
% version to rely on see: https://osf.io/hdxvb/


%% load all_trials

%clear; clc;
%addpath(genpath('/Volumes/daten/romy/eeglab13_5_4b'));

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

PATH = '/Volumes/daten/romy/';



SOURCEFILES = dir(strcat(PATH, 'Confidence/EEG_new/*.set')); %all brain vision header files in the folder
SUBJECTS = 1:numel(SOURCEFILES);
%%
Time_ERN=nan(65,500,40*250); % Matrix: electrode * time*participants*trials
Time_FRN=nan(65,500,40*250);
Time_SPN=nan(65,800,40*250);

contmat=zeros(40,3);

Trials=(1:250);
%%
for s =1:numel(SOURCEFILES)
    %% 
    
    vpn=SOURCEFILES(s).name(1:2); % read out of dataset name
    pcnt=(s*250)-250; % counter line ntrials
    
    EEG = pop_loadset('filename',sprintf('%s_clean.set',vpn),'filepath',sprintf('%s/Confidence/EEG_new/',PATH));

    
    EEGC=EEG; % save copy in temp memory
    
    %% special function to select only marker when another marker is given at a specific position

    % here: get button presses, only following tones
    [EEG, indices]=bool_epoch(EEG,[103] , [99],-1, [-0.2 .8]);
    
    
    contmat(s,1)=EEG.trials;
 
    %% Select epochs
    %  Trialnumber is always exactly 2 positions before the button press in the urevent structure 
    
    trigs=[{EEG.event(:).type}]';
    trig_db=zeros(size(trigs));
    bt=zeros(EEG.trials,3);
    for n=1:length(trigs)
        f=trigs{n};
        trig_db(n) = str2double(f(2:end));
    end
    err=0;
     % 
    for n=1:EEG.trials  %go trough all epochs  
        if EEG.epoch(n).eventtype{1} == 'S103' % find button presses
            k = EEG.epoch(n).eventurevent{1}; % get position in urevent structure
            if strcmp(EEG.urevent(k-1).type,'S 99') % check whether the previous trigger was a tone 
                f=EEG.urevent(k-2).type; % get the trigger 2 prior 
                tr=str2double(f(2:end)); % transform it into numeric  
            bt(n,1)=tr(1); % remember trial numbers for all variables
            bt(n,2)=n; % remember epoch number
            end 
            
        end  
        
     end
    
    %%
    % Artifact rejection: Tirals are not deleted but identified and not read out later. 
    % (Irej: index of rejected trials)
    
    [EEGn, Irej] = pop_eegthresh(EEG,1,[1:65] ,-150,150, -0.2,.798,1,0);
   
    %% %   >> [rej rejE] = rejtrend( signal, winsize, maxslope, minR, step);

    [rej, rejE] = rejtrend(EEG.data,500,50,0.3,1);
    % (rej: vector with 0 & ones... bad!)
    
    frej= find(rej==1);
    Out=sort([Irej,frej]);
    Trials2=ones(EEG.trials,1);
    tr_nr=[1:EEG.trials];
    for n=1: size(Out,2)
        if Out(n)>0
            Trials2(Out(n))=0; %for rejection
            tr_nr(Out(n))=0;  % for trigger assignment
        else
        end
    end
    %tr_nr=sort(tr_nr);
    contmat(s,3)=(EEG.trials-size(Out,2));
    
    %% baseline correction
    EEG = pop_rmbase( EEG, [-200 0]);
     %% 
    for n=1:size(tr_nr,2)
        if tr_nr(n) >0 && bt(n,1)<=bt(n,2)% skips artifact trials
        t=bt(tr_nr(n),1);
        try
                        Time_ERN(:,:,pcnt+t)=EEG.data(:,:,n);
        catch
        end
        end
    end    
    
    %% FRN
    EEG=EEGC;
    
     %% segment to feedback onset

    [EEG, indices]=pop_epoch( EEG, {'S 88'}, [-0.2 .8]);
    
    contmat(s,1)=EEG.trials;
 
    %% Auswahl der Epochen
    %  Trialnumber is always exactly 2 positions after the feedback in the urevent structure (except breaks)
    
    trigs=[{EEG.event(:).type}]';
    trig_db=zeros(size(trigs));
    bt=zeros(EEG.trials,3);
    for n=1:length(trigs)
        f=trigs{n};
        trig_db(n) = str2double(f(2:end));
    end
    err=0;
     % 
    for n=1:EEG.trials   
        if EEG.epoch(n).eventtype{1} == 'S 88' % find Feedback 
            k = EEG.epoch(n).eventurevent{1}; 
            if n < EEG.trials && (strcmp(EEG.urevent(k+1).type,'S 41') || strcmp(EEG.urevent(k+1).type,'S 51')|| strcmp(EEG.urevent(k+1).type,'S 61')) 
                if (tr==50 || tr==100 ||tr==150 || tr==200)
                f=EEG.urevent(k+3).type;
                else
                f=EEG.urevent(k+2).type; 
                end
                tr=str2double(f(2:end)); 
                bt(n,1)=tr-1;
                bt(n,2)=n; 
            elseif n==EEG.trials
                bt(n,1)=250;
                bt(n,2)=n;
            end 
            
        end  
        
     end
    
    %%
    % Artifact rejection: Tirals are not deleted but identified and not read out later. 
    % (Irej: index of rejected trials)
    
    [EEGn, Irej] = pop_eegthresh(EEG,1,[1:65] ,-150,150, -0.2,.798,1,0);
   
    %% %   >> [rej rejE] = rejtrend( signal, winsize, maxslope, minR, step);

    [rej, rejE] = rejtrend(EEG.data,500,50,0.3,1);
    % (rej: vector with 0 & ones... bad!)
    
    frej= find(rej==1);
    Out=sort([Irej,frej]);
    Trials2=ones(EEG.trials,1);
    tr_nr=[1:EEG.trials];
    for n=1: size(Out,2)
        if Out(n)>0
            Trials2(Out(n))=0; %for rejection
            tr_nr(Out(n))=0;  % for trigger assignment
        else
        end
    end
    %tr_nr=sort(tr_nr);
    contmat(s,3)=(EEG.trials-size(Out,2));
    
    %% baseline correction
    EEG = pop_rmbase( EEG, [-200 0]);
     %% 
    for n=1:size(tr_nr,2)
        if tr_nr(n) >0 && bt(n,1)<=bt(n,2)
        t=bt(tr_nr(n),1);
        try
            
            Time_FRN(:,:,pcnt+t)=EEG.data(:,:,n);
        catch
        end
        end
    end    
    
    
     %% SPN
    EEG=EEGC;
    
     %% 
    [EEG, indices]=pop_epoch( EEG, {'S 88'}, [-0.8 .8]);
    
    contmat(s,1)=EEG.trials;
 
%%
    
    trigs=[{EEG.event(:).type}]';
    trig_db=zeros(size(trigs));
    bt=zeros(EEG.trials,3);
    for n=1:length(trigs)
        f=trigs{n};
        trig_db(n) = str2double(f(2:end));
    end
    err=0;
     % 
    for n=1:EEG.trials  
        index = find(ismember(EEG.epoch(n).eventtype, 'S 88'));
            k = EEG.epoch(n).eventurevent{index}; 
            if n < EEG.trials && (strcmp(EEG.urevent(k+1).type,'S 41') || strcmp(EEG.urevent(k+1).type,'S 51')|| strcmp(EEG.urevent(k+1).type,'S 61')) % checken, ob der Trigger davor ein Ton war
                if (tr==50 || tr==100 ||tr==150 || tr==200)
                f=EEG.urevent(k+3).type;
                else
                f=EEG.urevent(k+2).type;
                end
                tr=str2double(f(2:end)); 
                bt(n,1)=tr-1; 
                bt(n,2)=n; 
            elseif n==EEG.trials
                bt(n,1)=250;
                bt(n,2)=n;
            end 
            
        %end  
        
     end
    
    %%
    % Artifact rejection: Tirals are not deleted but identified and not read out later. 
    % (Irej: index of rejected trials)
    
    [EEGn, Irej] = pop_eegthresh(EEG,1,[1:65] ,-150,150, -0.8,.798,1,0);
   
    %%

    [rej, rejE] = rejtrend(EEG.data,800,50,0.3,1);
    % (rej: vector with 0 & ones... bad!)
    
    frej= find(rej==1);
    Out=sort([Irej,frej]);
    Trials2=ones(EEG.trials,1);
    tr_nr=[1:EEG.trials];
    for n=1: size(Out,2)
        if Out(n)>0
            Trials2(Out(n))=0; %for rejection
            tr_nr(Out(n))=0;  % for trigger assignment
        else
        end
    end
    %tr_nr=sort(tr_nr);
    contmat(s,3)=(EEG.trials-size(Out,2));
    
    %% baseline correction
    EEG = pop_rmbase( EEG, [-500 -300]); % might have to be changed if it doesn't work well
     %% 
    for n=1:size(tr_nr,2)
        if tr_nr(n) >0 && bt(n,1)<=bt(n,2)% excludes artifacts
        t=bt(tr_nr(n),1);
        try
            Time_SPN(:,:,pcnt+t)=EEG.data(:,:,n);
        catch
        end
        end
    end    
    
    
    
end

%% Save matrices

save(sprintf('%sConfidence/Export/Time_FRN.mat', PATH), 'Time_FRN', '-v7.3');
save(sprintf('%sConfidence/Export/Time_ERN.mat', PATH), 'Time_ERN', '-v7.3');
save(sprintf('%sConfidence/Export/Time_SPN.mat', PATH), 'Time_SPN', '-v7.3');
