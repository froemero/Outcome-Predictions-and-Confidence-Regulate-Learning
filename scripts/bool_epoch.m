function [EEG, indices] = bool_epoch(EEG, target, cond_trig, pos, lim)
%ok. Aim of the present function is selecting epochs based on marker
%combinations, as by ABE's in BVA.
%therefore I first have to identify the position of targets in the EEG
%struct.

%%Example:
%target=224
%cond_trig=101

OLDEEG = EEG; %% maybe good having this


        %% search for triggers in EEG.event.type
        %read in markers
        trigs=[{EEG.event(:).type}]';

        
        
        %transfer Markers into doubles for comparison
        trig_db=[];
        
        n=1;
        for n=1:length(trigs)
            f=trigs{n};
            trig_db(n,1) = str2double(f(2:end));
        end
        try
        targall=size([find(trig_db==target(1));find(trig_db==target(2))],1); %returns number of targets 
        % doesn't work for only 1 target so:
        catch targall=size(find(trig_db==target(1)),1); % get only target 1
        end
%% select critical cases and read out linenumbers
trgl=[];
ctrgl=[];


ct1=size(find(trig_db==cond_trig(1)),1);
try
    ct2=size(find(trig_db==cond_trig(2)),1); % return 2nd condition triggers if any
catch ct2=0;
end

 
        try
        ftargp=sort([find(trig_db==target(1));find(trig_db==target(2))]); % searches for targets
        catch ftargp=find(trig_db==target(1)); % returns target lines
        end        
       
        for n=1:length(ftargp) % goes through all target lines
                if trig_db(ftargp(n)+pos)==cond_trig(1);
                    %checks whether there is a condition trigger at the
                    %specified position
                    %if so, it saves targetline
                    trgl=[trgl ftargp(n)];
                end
        end
  

fprintf('bool_epoch():%d triggers found\n',targall);
fprintf('bool_epoch():%d condition 1 and %d condition 2 found\n',ct1, ct2);
fprintf('bool_epoch():%d matching events selected\n', length(trgl));


%% Struktur vorbereiten

g = [];
try, g.epochfield; 	 	  catch, g.epochfield = 'type'; end; % obsolete
try, g.timeunit; 	 	  catch, g.timeunit = 'points'; end;
try, g.verbose; 	      catch, g.verbose = 'on'; end;
try, g.newname; 	      catch, g.newname = fastif(isempty(EEG.setname), '', [EEG.setname ' epochs' ]); end;
try, g.eventindices;      catch, g.eventindices = 1:length(EEG.event); end;
try, g.epochinfo;         catch, g.epochinfo = 'yes'; end;
try, if isempty(g.valuelim), g.valuelim = [-Inf Inf]; end; catch, g.valuelim = [-Inf Inf]; end;

tmpevent = EEG.event;
tmpeventlatency = [ tmpevent(:).latency ];
[tmpeventlatency Itmp] = sort(tmpeventlatency);
EEG.event = EEG.event(Itmp);  % sort by ascending time
Ievent = g.eventindices;

%% select epochs for those targets matching the criteria 
alllatencies = tmpeventlatency(trgl);


% taken from pop_epoch:

switch lower( g.timeunit )
    case 'points',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, [lim(1) lim(2)]*EEG.srate, ...
            'valuelim', g.valuelim, 'allevents', tmpeventlatency);
        tmptime = tmptime/EEG.srate;
    case 'seconds',	[EEG.data tmptime indices epochevent]= epoch(EEG.data, alllatencies, lim, 'valuelim', g.valuelim, ...
            'srate', EEG.srate, 'allevents', tmpeventlatency);
    otherwise, disp('pop_epoch(): invalid event time format'); beep; return;
end;
alllatencies = alllatencies(indices);
fprintf('pop_epoch():%d epochs generated\n', length(indices));

% update other fields
% -------------------
if lim(1) ~= tmptime(1) & lim(2)-1/EEG.srate ~= tmptime(2)
    fprintf('pop_epoch(): time limits have been adjusted to [%3.3f %3.3f] to fit data points limits\n', ...
        tmptime(1), tmptime(2)+1/EEG.srate);
end;
EEG.xmin = tmptime(1);
EEG.xmax = tmptime(2);
EEG.pnts = size(EEG.data,2);
EEG.trials = size(EEG.data,3);
EEG.icaact = [];
if ~isempty(EEG.setname)
    if ~isempty(EEG.comments)
        EEG.comments = strvcat(['Parent dataset "' EEG.setname '": ----------'], EEG.comments);
    end;
    EEG.comments = strvcat(['Parent dataset: ' EEG.setname ], ' ', EEG.comments);
end;
%try, g.newname; 	      catch, g.newname = fastif(isempty(EEG.setname), '', [EEG.setname ' epochs' ]); end;
EEG.setname = g.newname;

% count the number of events to duplicate and duplicate them
% ----------------------------------------------------------
totlen = 0;
for index=1:EEG.trials, totlen = totlen + length(epochevent{index}); end;
EEG.event(1).epoch = 0;          % create the epoch field (for assignment consistency afterwards)
if totlen ~= 0
    newevent(totlen) = EEG.event(1); % reserve array
else
    newevent = [];
end;

% modify the event structure accordingly (latencies and add epoch field)
% ----------------------------------------------------------------------
allevents = [];
count = 1;
for index=1:EEG.trials
    for indexevent = epochevent{index}
        newevent(count)         = EEG.event(indexevent);
        newevent(count).epoch   = index;
        newevent(count).latency = newevent(count).latency - alllatencies(index) - tmptime(1)*EEG.srate + 1 + EEG.pnts*(index-1);
        count = count + 1;
    end;
end;
EEG.event = newevent;
EEG.epoch = [];
EEG = eeg_checkset(EEG, 'eventconsistency');

% check for boundary events
% -------------------------
disp('pop_epoch(): checking epochs for data discontinuity');
if ~isempty(EEG.event) & isstr(EEG.event(1).type)
    tmpevent = EEG.event;
    boundaryindex = strmatch('boundary', { tmpevent.type });
    if ~isempty(boundaryindex)
        indexepoch = [];
        for tmpindex = boundaryindex
            indexepoch = [indexepoch tmpevent(tmpindex).epoch ];
        end;
        EEG = pop_select(EEG, 'notrial', indexepoch);
        % update the "indices of accepted events", too
        indices = indices(setdiff(1:length(indices),indexepoch));
    end;
end;



end



