%-----------------------------------------------------------------------------------------------------------------------
%-- Logistic regression of psychophysical data
%--	ZZ 20210808
% %% ZZ 20210808 add logistic regression
% refer to Tsunada, 2019, eLife & Akrami, 2018, Nature
%-----------------------------------------------------------------------------------------------------------------------
function Psycho_logistic_ZZ(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, batch_flag)

TEMPO_Defs;
Path_Defs;


% get task information
% temp_microstim = data.misc_params(MICROSTIM, :);
if length(munique(data.misc_params(MICROSTIM,:)'))>1
    error('This is a block with Microstimulation !')
end

temp_stim_type = data.moog_params(STIM_TYPE, :, MOOG);
temp_heading = data.moog_params(HEADING, :, MOOG);
% temp_motion_coherence = data.moog_params(COHERENCE, :, MOOG);   % only 1 coherence in my task
temp_outcome = data.misc_params(OUTCOME, :);
% temp_mask_status = data.moog_params(MASK_STATUS, :, MOOG);

% If length(BegTrial) > 1 and all elements are positive, they are trials to be included.
% Else, if all elements are negative, they are trials to be excluded.
% This enable us to exclude certain trials ** DURING ** the recording more easily. HH20150410
select_trials = false(size(temp_heading));
if length(BegTrial) == 1 && BegTrial > 0 % Backward compatibility
    select_trials(BegTrial:EndTrial) = true;
elseif all(BegTrial > 0) % To be included
    select_trials(BegTrial) = true;
elseif all(BegTrial < 0) % To be excluded
    select_trials(-BegTrial) = true;
    select_trials = ~ select_trials;
else
    disp('Trial selection error...');
    keyboard;
end

% trials = 1:length(select_trials);

stim_type = temp_stim_type(select_trials);
heading = temp_heading(select_trials);
% motion_coherence = temp_motion_coherence(select_trials);
outcome = temp_outcome(select_trials);
% mask_status = temp_mask_status(select_trials);

unique_stim_type = munique(stim_type');
unique_heading = munique(heading');
% unique_motion_coherence = munique(motion_coherence');
% unique_mask_status = munique(mask_status');

one_repetition = length(unique_heading) * length(unique_stim_type);
repetitionN = floor(length(heading)/one_repetition);

LEFT = 1;
RIGHT = 2;

event_in_bin = squeeze(data.event_data(:,:,select_trials))';

choice_per_trial = LEFT * squeeze(sum(event_in_bin == IN_T2_WIN_CD,2))' + RIGHT*squeeze(sum(event_in_bin == IN_T1_WIN_CD, 2))';
if length(unique(choice_per_trial)) > 2
    disp('Neither T1 or T2 chosen / More than one target chosen.  This should not happen! File must be bogus.');
    fprintf('%g cases...\n', sum(choice_per_trial==3));
    beep;
end

for condition = 1:3
    if sum(unique_stim_type == condition) == 0
%         correct_rate{condition} = nan;
        beta{condition} = nan(2,1);
        p_value{condition} = nan(2,1);
        beta_with_hist{condition} = nan(4,1);
        p_value_with_hist{condition} = nan(4,1);
        fit_data_psycho_cum{condition} = nan(4,1);
%         lapse = nan;
        continue;
    end
        
    trials_select_cc = logical(stim_type==condition);  % trials number for the current condition
    trial_ind = find(trials_select_cc ==1);  % which trial is in this condition
    
%%%% ========= Psychometric Function Fitting =========== %%%%%%
    for h = 1:length(unique_heading)
        trials_select_ckh = trials_select_cc & (heading == unique_heading(h));
        
        rightward_trials  = (trials_select_ckh & choice_per_trial == RIGHT);
        rightward_num = sum(rightward_trials);
        
        fit_psycho_cum_lapse{condition}(h,1) = unique_heading(h);
        fit_psycho_cum_lapse{condition}(h,2) = rightward_num;
        fit_psycho_cum_lapse{condition}(h,3) = sum(trials_select_ckh);
    end

%%%%%%% =======  Lapse ======= %%%%%%%
    %     % lapse calculated by its definition, but left and right lapses are default same
%     lapse= (sum(y_choice(heading_cc==min(unique_heading))==RIGHT) + sum(y_choice(heading_cc==max(unique_heading))==LEFT)) / 2*repetitionN;

% === Lapse is determined by fitting now   ZZ 20210913
% inspired by the Pisupati et al., 2021, elife 
% and the toolkit used for fitting is from Prins & Kingdom, 2018
[params_psycho] = cum_gauss_lapse_comparison(fit_psycho_cum_lapse{condition}, batch_flag); 

    
    %%%%%% ======  Constructing regressor matrix =====%%%%%%%%%%
    heading_cc = heading(trial_ind);    % heading of the current conditions
    outcome_pt = nan(1, length(trial_ind));  % previous trial outcome
    rewarded_pt = outcome_pt;   unrewarded_pt = outcome_pt;
    choice_pt = outcome_pt;  % previous trial choice
    heading_pt = rewarded_pt;
    
    % Previous trial information, choice bias in current trial conditioned
    % on the outcome of the previous trial
    if sum(trial_ind-1 == 0) >0
        rewarded_pt(1) = 0;        % 1, reward right choice; -1, reward left choice; 0, no reward;
        unrewarded_pt(1) = 0;    % 1, not reward right choice; -1, not reward left choice; 0, rewarded;
        outcome_pt(2:end) = outcome(trial_ind(2:end)-1);
        choice_pt(2:end) = choice_per_trial(trial_ind(2:end)-1);
        heading_pt(2:end) = heading(trial_ind(2:end)-1);
    else
        outcome_pt = outcome(trial_ind-1);
        choice_pt = choice_per_trial(trial_ind-1);
        heading_pt = heading(trial_ind-1);

    end
    
    % previous-trial rewarded regressor
    rewarded_pt(outcome_pt==5) = 0;
    rewarded_pt(outcome_pt==0 & choice_pt==RIGHT) = 1;    % 1, reward right choice; -1, reward left choice; 0, no reward;
    rewarded_pt((outcome_pt==0 & choice_pt==LEFT)) = -1;
    % previous-trial unrewarded regressor
    unrewarded_pt(outcome_pt==0) = 0;
    unrewarded_pt(outcome_pt==5 & choice_pt==RIGHT) = 1;   % 1, not reward right choice; -1, not reward left choice; 0, rewarded;
    unrewarded_pt(outcome_pt==5 & choice_pt==LEFT) = -1;
    
    
    % ==========  Logistic glmfit input ======== %
    x_regres = [heading_cc' rewarded_pt' unrewarded_pt'];
    
    y_choice = choice_per_trial(trial_ind);   
    y_right = y_choice==RIGHT;  y_left= y_choice==LEFT; 
    y_regres = double(y_right'); 
    
    % Do not consider reward history
    [bb, stats] = logistic_with_lapse(x_regres(:,1),y_regres,params_psycho(3:4)); 
    beta{condition} = bb;
    p_value = stats.p; 
    
    [bb,stats] = logistic_with_lapse(x_regres,y_regres,params_psycho(3:4)); 
    beta_with_hist{condition} = bb; 
    p_value_with_hist{condition} = stats.p; 
    
    fit_data_psycho_cum{condition} = params_psycho;
end


%% Data Saving

if ~isempty(batch_flag)
    config.batch_flag = batch_flag;
    
    % Time
    time = data.htb_header{1}.date;
    space_pos = strfind(time, ' ');
    date = [time(space_pos(1)+1:space_pos(2)-1) '-' time(space_pos(2)+1:space_pos(3)-1) '-' time(space_pos(4)+1:end)];
%     outpath = ['Z:\Data\Tempo\Batch\' batch_flag(1:end-2) '\'];
%     
%     % Check directory
%     if ~exist(outpath,'dir')
%         mkdir(outpath);
%     end
%     
%     dot_pos = strfind(FILE,'.');
%     if ~isempty(dot_pos)
%         fileName = FILE(1:dot_pos-1);
%     end
%     
%     savefilename = [outpath [fileName '_' num2str(SpikeChan)] '_' 'logist'];
    
%     result = PackResult(fileName, SpikeChan, repetitionN, unique_stim_type,...
%         beta,p_value,beta_with_hist,p_value_with_hist,fit_data_psycho_cum);
    
    result = PackResult(date,FILE, SpikeChan, repetitionN, unique_stim_type,...
        beta,p_value,beta_with_hist,p_value_with_hist,fit_data_psycho_cum);

    
    config.suffix = 'logistc';
    
    config.save_figures = [];
    
    config.sprint_once_marker = [];
    config.sprint_once_contents = []; 
    
    config.sprint_loop_marker = [];
    config.sprint_loop_contents = [];
%     save(savefilename, 'result');
    SaveResult(config,result);
end

end