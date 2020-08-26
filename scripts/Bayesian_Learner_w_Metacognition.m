%% Bayesian learner with Metacognition:
% Leverages insights about its own performance to optimize learning from
% feedback.


% Goal:
% implement 3 versions:
% 1) just feedback based learning
% 2) just error monitoring and feedback based learning
% 3) error monitoring with Confidence and feedback based learning

%%

clear; clc; close all

%% Task variables


tVal    =20; % target
sVal    =45; % scale formerly
% note that these will just be used as indices later and are not the actual
% values
models =[1:3];

for model = 1:length(models)
    %% Simulation parameters
    
    numIts  =500;        % number of iterations
    TT      =250;         % number of trials
    
    if model ==1 % learns from feedback only
        BayesOptimal =0; % decides whether prediction is biased towards zero during updating or just weighted by reliability
        NoConfidenceUse=1;
        useMonitoring =0;
        
    elseif model ==2 % learns from feedback using PM
        BayesOptimal =0; % decides whether prediction is biased towards zero during updating or just weighted by reliability
        NoConfidenceUse=1;
        useMonitoring =1;
        
    else % uses everything available
        BayesOptimal =1; % decides whether prediction is biased towards zero during updating or just weighted by reliability
        NoConfidenceUse=0;
        useMonitoring =1;
    end
    
    
    %% Learning variables
    
    p_t=linspace(0, 100); % range of means
    p_s=linspace(.1, 100, 50); % range of widths
    
    allEcNoise=[5, 10, 20]; % how noisy the efference copy can be
    possConf = [3, 2, 1]; % high, medium, low confidence  - to be derived from the estimate of what the reliability of the EC is
    %% Learner characteristics
    
    sigma_r=10; % how noisy execution is
    allMCLs =[0, .75, 1]; % levels of metacognitive ability (zero to 1; zero: no clue; one: knows exactly what's going on)
    
    
    %% Simulation
    
    % GOAL: infer s & t
    for subType=1:length(allMCLs) % simulating different types of people
        metaCogLevel=allMCLs(subType);  % 1 = we totally know how good our efference copy is,
        
        for it=1:numIts
            % generate an actual task
            t=p_t(tVal);
            s=p_s(sVal);
            
            pFgivAllElse=nan(length(p_t), length(p_s), length(allEcNoise)); % prepopulating probability matrix w NaNs
            
            for trial=1:TT % loop through trials
                
                % Set prior
                if trial>1
                    priorTandS=postTandS; % belief what the scale is
                else
                    % specify prior distribution over S,T:
                    priorTandS=ones(length(p_t), length(p_s)); % initially everything is equally likely
                    priorTandS=priorTandS./(sum(priorTandS(:)));
                end
                
                
                %% PRODUCE Timing based on expectation over T --> try to produce
                % most likely response given what you've learned
                margPostT=sum(priorTandS, 2);
                if  ~isfinite(margPostT)
                    keyboard
                end
                
                exp_t=sum(margPostT.*p_t'); % best guess of what the target is
                i=exp_t; % intended response (i.e. what we want to do)
                % Produce a timing (with some noise):
                r= i+normrnd(0,sigma_r);  % Actual timing production
                
                %% Get efference copy of the executed response
                % only used in Model 2 and 3
                % choose a level of efference copy accuracy:
                Q=randperm(length(allEcNoise), 1); % did I pay attention/have a good sense this time or no?
                ecNoise=allEcNoise(Q); % choose EC noise from our list of possible levels;
                if useMonitoring
                    % efference copy
                    c   = r+ normrnd(0,ecNoise);
                else
                    c = i;
                end
                
                
                %% Equation (1) Compute feedback based on the timing you produced:
                f = (r-t).*s; % actual feedback
                
                %% read-out of efference copy quality
                % only used in model 3
                
                trueECQual = 1:length(allEcNoise)==Q;
                falseECQual = randsample(find(~trueECQual), 1);
                
                
                % determine whether their assessment of Q is veridical or no
                if BayesOptimal
                    if rand <metaCogLevel
                        subQ = Q;
                    else
                        subQ = falseECQual;
                    end
                else
                    subQ = randsample([1:3], 1);% some random confidence rating
                end
                
                if ~NoConfidenceUse
                    sigma_c = allEcNoise(subQ); % ec noise, conditional on subQ, also known as confidence!
                else
                    sigma_c = mean(allEcNoise);   % prediction model has some sense how accurate the EC is overall, but not by trial
                end
                %% Updating
                % model 1 and 2
                
                if  useMonitoring
                    v = (1./(1./sigma_c.^2 + 1./sigma_r.^2)); % model 3 equation (6)
                else
                    v = (1./(1./sigma_r.^2)); % model 1 expected variance in feedback entirely driven by response noise equation (11)
                end
                
                % setting expectation depending on available information
                if BayesOptimal
                    % equation (8)
                    m = i.*v*(1/sigma_r.^2)+c.*v*(1/sigma_c.^2); % model 3 - all info is available equation (7)
                    % we weight what we wanted to do with what we think we did
                    % v*(1/productionNoise.^2) and v*(1/EcNoise.^2) always sum up
                    % to 1, so this is really, how much to rely on the efference
                    % copy and how much to believe that what one did was hight
                else
                    if useMonitoring
                        % alternative with efference copy as is, rather than weighted inference over efference copy
                        m = i.*v*(1/sigma_r.^2)+c.*v*(1/sigma_c.^2); % takes average reliability of efference copy as sigma_c
                    
                    else
                        m=i; % we think we did what we wanted to do.
                    end
                end
                
                % Building the joint distribution over N(c;r,?_c^2)N(r;i,?_r^2 )p(s,t)
                meanOfExpectedF=(m-p_t)'*p_s;
                stdOfEpectedF  = sqrt(v).*repmat(p_s, length(p_t), 1); % Now likelihood is gaussian with mean zero and scale p_S
                % equation (8)
                pFgivAllElse=normpdf(f, meanOfExpectedF, stdOfEpectedF); % hey, we actually got some feedback now!
                
                % Get posterior by multiplying prior * likelihood (and normalizing)
                postTandS=priorTandS.*pFgivAllElse;
                postTandS=postTandS./(sum(postTandS(:)));
                
                
                %% Store some things
                
                triali(trial)= i;
                
                % objective error
                trialsignedError(trial)=(r-t).*s; % how wrong was what they did?
                trialError(trial)=abs(trialsignedError(trial)); % error magnitude
                
                % latent variables:
                trialGoalError(trial)=(abs(i-t)).*s; % how wrong was what they thought they should do
                trialECnoise(trial)=Q; % Q = true category of ecnoise
                trialmetacogLevel(trial)=metaCogLevel; % is this agent well calibrated or no?
                
                % get the best guess of scale to give prediction rating
                margPostS=sum(priorTandS, 1);
                expS=sum(margPostS.*p_s); % best guess of what the scale is
                curS(trial)=expS;
                margpostS=sum(postTandS, 1);
                exppostS=sum(margpostS.*p_s);
                deltas(trial) = expS - exppostS;
                margPostT=sum(postTandS, 2);
                exp_t=sum(margPostT.*p_t'); % best guess of what the target is now
                nexti=exp_t;
                deltat(trial) = i-nexti;
                % prediction
                trialPredictedError(trial)=abs(c-i).*expS; % expected error given efference copy, intended response and best guess about scale
                % equation (9)
                trialsignedPredictedError(trial)=(c-i).*expS; % same but signed
                % confidence aka co
                trialConfidence(trial)= possConf(subQ); % confidence given belief about ec noise
                
            end
            % store trial variables across iterations
            % objective error
            allsignedErrors(it,:,subType)=trialsignedError;
            allErrors(it,:,subType)=trialError;
            % latent variables
            allGoalErrors(it,:,subType)=trialGoalError;
            alltrialECnoise(it,:,subType)=trialECnoise;
            alltrialMetacoglevel(it,:,subType)=trialmetacogLevel;
            allcurS(it,:,subType)= curS;
            % Prediction
            allsignedPredictedError(it,:,subType)=trialsignedPredictedError;
            allPredictedError(it,:,subType)=trialPredictedError;
            % Confidence
            alltrialConfidence(it,:,subType)=trialConfidence;
            alltrialtupdate(it,:,subType)= deltat;
            alltrialsupdate(it,:,subType)= deltas;
            
        end
    end
    
    
    
    
    %%
    
    AllPerf = squeeze(nanmean(allErrors));
    AllEst = squeeze(nanmean(allGoalErrors));
    AllBlock = repelem([1:5], 50);
    
    AlBlockIndex = repmat(AllBlock, [numIts,1,3]);
    
    %% Figure 1 E
    figure()
    set(gcf,'Color', [1 1 1])
    hold on
    linea=plot(allsignedPredictedError(alltrialConfidence==1  ),allsignedErrors(alltrialConfidence==1  ),'O','MarkerSize',2)
    
    lineb=plot(allsignedPredictedError(alltrialConfidence==2  ),allsignedErrors(alltrialConfidence==2  ),'O','MarkerSize',2)
    
    linec=plot(allsignedPredictedError(alltrialConfidence==3  ),allsignedErrors(alltrialConfidence==3  ),'O','MarkerSize',2)
    xx=lsline
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('Predicted Outcome','FontSize',12,'FontWeight','bold')
    ylabel('Actual Outcome','FontSize',12,'FontWeight','bold')
    ylim([-6000, 6000])
    xlim([-6000, 6000])
    hline = refline([0 0]);
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'southeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    saveas(gcf, sprintf('Prediction_Outcome_Confidence_model_%d.pdf', model))
    
    
    %% Updating
    
    figure()
    set(gcf,'Color', [1 1 1])
    subplot(3,2,1)
    hold on
    linea=plot(abs(allsignedErrors(alltrialConfidence==1  )-allsignedPredictedError(alltrialConfidence==1  )),abs(alltrialsupdate(alltrialConfidence==1  )),'O','MarkerSize',2)
    
    lineb=plot(abs(allsignedErrors(alltrialConfidence==2  )-allsignedPredictedError(alltrialConfidence==2  )),abs(alltrialsupdate(alltrialConfidence==2  )),'O','MarkerSize',2)
    
    linec=plot(abs(allsignedErrors(alltrialConfidence==3  )-allsignedPredictedError(alltrialConfidence==3  )),abs(alltrialsupdate(alltrialConfidence==3  )),'O','MarkerSize',2)
    xx=lsline
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('SPE','FontSize',12,'FontWeight','bold')
    ylabel('Scale Update','FontSize',12,'FontWeight','bold')
    ylim([0, 70])
    xlim([0, 10000])
    hline = refline([0 0]);
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'northeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    
    subplot(3,2,2)
    hold on
    linea=plot(abs(allsignedErrors(alltrialConfidence==1  )-allsignedPredictedError(alltrialConfidence==1  )),abs(alltrialtupdate(alltrialConfidence==1  )),'O','MarkerSize',2)
    
    lineb=plot(abs(allsignedErrors(alltrialConfidence==2  )-allsignedPredictedError(alltrialConfidence==2  )),abs(alltrialtupdate(alltrialConfidence==2  )),'O','MarkerSize',2)
    
    linec=plot(abs(allsignedErrors(alltrialConfidence==3  )-allsignedPredictedError(alltrialConfidence==3  )),abs(alltrialtupdate(alltrialConfidence==3  )),'O','MarkerSize',2)
    xx=lsline
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('SPE','FontSize',12,'FontWeight','bold')
    ylabel('Target Update','FontSize',12,'FontWeight','bold')
    ylim([0, 50])
    xlim([0, 10000])
    hline = refline([0 0]);
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'northeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    subplot(3,2,3)
    hold on
    linea=plot([1, 2, 3], [mean(abs(alltrialsupdate(alltrialConfidence==1  ))), mean(abs(alltrialsupdate(alltrialConfidence==2  ))),mean(abs(alltrialsupdate(alltrialConfidence==3  )))],'O','MarkerSize',5)
   
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('Confidence Level','FontSize',12,'FontWeight','bold')
    ylabel('Scale Update','FontSize',12,'FontWeight','bold')
    ylim([0.2, 0.8])
    xlim([0, 4])
    xticklabels({'','low','medium','high'})
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'northeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    subplot(3,2,4)
    hold on
    linea=plot([1, 2, 3], [mean(abs(alltrialtupdate(alltrialConfidence==1  ))), mean(abs(alltrialtupdate(alltrialConfidence==2  ))),mean(abs(alltrialtupdate(alltrialConfidence==3  )))],'O','MarkerSize',5)
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('Confidence Level','FontSize',12,'FontWeight','bold')
    ylabel('Target Update','FontSize',12,'FontWeight','bold')
    ylim([0.2, 0.8])
    xlim([0, 4])
    xticklabels({'','low','medium','high'})
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'northeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    
    subplot(3,2,5)
    hold on
    linea=plot(abs(allsignedErrors(alltrialConfidence==1  )),abs(alltrialsupdate(alltrialConfidence==1  )),'O','MarkerSize',2)
    
    lineb=plot(abs(allsignedErrors(alltrialConfidence==2  )),abs(alltrialsupdate(alltrialConfidence==2  )),'O','MarkerSize',2)
    
    linec=plot(abs(allsignedErrors(alltrialConfidence==3  )),abs(alltrialsupdate(alltrialConfidence==3  )),'O','MarkerSize',2)
    xx=lsline
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('Error Magnitude','FontSize',12,'FontWeight','bold')
    ylabel('Scale Update','FontSize',12,'FontWeight','bold')
    ylim([0, 70])
    xlim([0, 10000])
    hline = refline([0 0]);
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'northeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    subplot(3,2,6)
    hold on
    linea=plot(abs(allsignedErrors(alltrialConfidence==1  )),abs(alltrialtupdate(alltrialConfidence==1  )),'O','MarkerSize',2)
    
    lineb=plot(abs(allsignedErrors(alltrialConfidence==2  )),abs(alltrialtupdate(alltrialConfidence==2  )),'O','MarkerSize',2)
    
    linec=plot(abs(allsignedErrors(alltrialConfidence==3  )),abs(alltrialtupdate(alltrialConfidence==3  )),'O','MarkerSize',2)
    xx=lsline
    xx(1).LineWidth=3;
    xx(2).LineWidth=3;
    xx(3).LineWidth=3;
    xlabel('Error Magnitude','FontSize',12,'FontWeight','bold')
    ylabel('Target Update','FontSize',12,'FontWeight','bold')
    ylim([0, 50])
    xlim([0, 10000])
    hline = refline([0 0]);
    lgd=legend([xx], {'high', 'medium', 'low'}, 'Box', 'off');
    lgd.Location = 'northeast';
    lgd.FontSize =12;
    title(lgd,'Confidence','FontSize',14)
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    saveas(gcf, sprintf('Updating_model_%d.pdf', model))
    %% Figure 1 F
    
    figure()
    set(gcf,'Color', [1 1 1])
    hold on
    linea=plot(log(squeeze(nanmean(allGoalErrors))),'LineWidth',2)
    xlabel('Trial','FontSize',12,'FontWeight','bold')
    ylabel('Target Error Magnitude','FontSize',12,'FontWeight','bold')
    lgd=legend([linea], {'low', 'medium', 'high'}, 'Box', 'off')
    lgd.FontSize =12;
    title(lgd,sprintf('Confidence\n Calibration'),'FontSize',14)
    xlim([0,251])%
    ylim([3, 9])
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    saveas(gcf, sprintf('Estimate_model_%d.pdf', model))
    
    
    

        %% collapsed across calibration levels that don't matter
        figure()
        set(gcf,'Color', [1 1 1])
        hold on
        linea=plot(log(nanmean(squeeze(nanmean(allGoalErrors)),2)),'LineWidth',2)
        xlabel('Trial','FontSize',12,'FontWeight','bold')
        ylabel('Target Error Magnitude','FontSize',12,'FontWeight','bold')
        xlim([0,251])%
        ylim([3, 9])
        set(gca,'TickDir','out', 'Box', 'off')
        set(gcf, 'Renderer', 'Painters');
        saveas(gcf, sprintf('Estimate_model_%d_collapsed.pdf', model))
        

    
    %% scale representation
    figure()
    set(gcf,'Color', [1 1 1])
    hold on
    linea=plot(log(squeeze(nanmean(abs(allcurS - s)))),'LineWidth',2)
    xlabel('Trial','FontSize',12,'FontWeight','bold')
    ylabel('Scale Error Magnitude','FontSize',12,'FontWeight','bold')
    lgd=legend([linea], {'low', 'medium', 'high'}, 'Box', 'off')
    lgd.FontSize =12;
    title(lgd,sprintf('Confidence\n Calibration'),'FontSize',14)
    xlim([0,251])%
    ylim([0.5, 4])
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    saveas(gcf, sprintf('Estimate_scale_model_%d.pdf', model))
    
    
      %% scale representation
    figure()
    set(gcf,'Color', [1 1 1])
    hold on
    linea=plot(log(nanmean(squeeze(nanmean(abs(allcurS - s))),2)),'LineWidth',2)
    xlabel('Trial','FontSize',12,'FontWeight','bold')
    ylabel('Scale Error Magnitude','FontSize',12,'FontWeight','bold')
    xlim([0,251])%
    ylim([0.5, 4])
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    saveas(gcf, sprintf('Estimate_scale_model_%d_collapsed.pdf', model))
    
    %% scale and target updates

    figure()
    set(gcf,'Color', [1 1 1])
    subplot(1,2,1)
    hold on
    linea=plot(log((squeeze(nanmean(abs(alltrialtupdate))))),'LineWidth',2)
    xlabel('Trial','FontSize',12,'FontWeight','bold')
    ylabel('log Target Update','FontSize',12,'FontWeight','bold')
    lgd=legend([linea], {'low', 'medium', 'high'}, 'Box', 'off')
    lgd.FontSize =12;
    title(lgd,sprintf('Confidence\n Calibration'),'FontSize',14)
    xlim([0,251])%
    ylim([-5, 4])
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    subplot(1,2,2)
    hold on
    linea=plot(log(squeeze(nanmean(abs(alltrialsupdate)))),'LineWidth',2)
    xlabel('Trial','FontSize',12,'FontWeight','bold')
    ylabel('log Scale Update','FontSize',12,'FontWeight','bold')
    lgd=legend([linea], {'low', 'medium', 'high'}, 'Box', 'off')
    lgd.FontSize =12;
    title(lgd,sprintf('Confidence\n Calibration'),'FontSize',14)
    xlim([0,251])%
    ylim([-5, 4])
    set(gca,'TickDir','out', 'Box', 'off')
    set(gcf, 'Renderer', 'Painters');
    saveas(gcf, sprintf('Updates_model_%d.pdf', model))
end




