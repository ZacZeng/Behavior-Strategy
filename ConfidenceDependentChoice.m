% This part is aim to separate the whole psychometric curve according to
% privious correct choice and visualize Confidence-dependent choice updating
% Referring to Lak et al., 2020, eLife
% When running this code, the current folder need contain all behavior data

function ConfidenceDependentChoice(~)

pathname = uigetdir(cd, 'Choose a folder');
if pathname ==0
    msgbox('You did not choose a correct folder');
    return;
else
    cd(pathname);
end

% confidence_dependent_choice();

%     function confidence_dependent_choice(~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

method = 0;   % 0: Maximum likelihood;  1: Square error
tolerance = 10;
LEFT = 1; RIGHT =2;
% I only choose behavior data after monkey well-trianed and pool them as
% ground mean
grand_heading = [-10 -5 -2.5 -1.25 0 1.25 2.5 5 10];  % Maybe having different heading sets, but I coarsely align them
hhi = -10: 0.05:10;  % for psychometric curve plotting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pool all well-trained sessions
fileName_list = dir ( '*_Psycho.mat*');

% Initiation
psycho_fit = arrayfun(@(x) zeros(9,3), [1:3]', 'UniformOutput', 0);  % For grand psychometric curve fitting, 3 conditions default
psycho_fit_1st_half = psycho_fit; psycho_fit_2nd_half = psycho_fit;  % Divide sessioins as two halves
psycho_fit_hist = repmat(psycho_fit,1,3);  %  History-dependent grand psychometric curves
cur_depend_on_pre_unit = arrayfun(@(x) zeros(9,3), [1:10]', 'UniformOutput', 0);  % current trial depdent on condition and correct choice of previous trial
cur_depend_on_pre = arrayfun(@(x) cur_depend_on_pre_unit, reshape(1:9, 3,3), 'UniformOutput', 0);  % 3 pre_conds * 3 cur_conds
cur_depend_on_pre_1st = cur_depend_on_pre;  cur_depend_on_pre_2nd = cur_depend_on_pre;
fit_hist_psycho = cell(length(fileName_list),1);  temp_hist_psy_fit = fit_hist_psycho;

% Pool all sessions
for sess = 1:length(fileName_list)
    
    load(fileName_list(sess).name);
    raw = result.raw;
    
    unique_conditions = unique(raw(:,1));
    unique_headings = unique(raw(:, 3));
    pre_headings = [unique_headings(1:5); 0 ;unique_headings(6:end)];   % Divide pre 0 heading into left choice and right choice
    
    cur_cond_pre_cond_ind = cell(length(unique_conditions), length(unique_conditions));
    
    if length(unique_conditions) ~= 3  % I only consider sessions containing 3 stimulus conditions, of course this is unnecessary if we pool all sessions
        continue;
    end
    
    if length(unique_headings) ~= 9  % I only consider sessions containing 9 headings, of course this is alse unnecessary
        continue;
    end
    
    if result.repetitionN < 15 % This constraint may be necessary if we divide sessions into two parts
        continue;
    end
    
    cur_trials = raw(2:end, :);     % current trials considered need minus 1 if we look back the previous trial
    pre_trials = raw(1:end-1, :);
    
    % Divide current trials into two halves
    sess_1st_half = cur_trials(1:floor(length(cur_trials)/2), :);
    sess_2nd_half = cur_trials(length(sess_1st_half)+1:end, :);
    
    % ===================== Grand mean as ground thruth================================
    % Pool all sessions regardless of previous-trial history
    for cc = 1:length(unique_conditions)    % cc: current condition
        this_condition = raw(:,1)==unique_conditions(cc);
        this_condition_1st_half = sess_1st_half(:,1)==unique_conditions(cc);  % compare the grand mean curve with the grand first half and the rest half.
        this_condition_2nd_half = sess_2nd_half(:,1)==unique_conditions(cc);
        
        for ch = 1:length(unique_headings)  % ch: current headings
            this_heading = raw(:,3)==unique_headings(ch);
            this_heading_1st_half = sess_1st_half(:,3)==unique_headings(ch);
            this_heading_2nd_half = sess_2nd_half(:,3)==unique_headings(ch);
            
            temp_hist_psy_fit{sess}{cc}(ch,2) = sum(raw(this_condition & this_heading, 4) == RIGHT);
            temp_hist_psy_fit{sess}{cc}(ch,3) =sum(this_condition & this_heading);
            
            psycho_fit{cc}(ch,2) = psycho_fit{cc}(ch,2) +  temp_hist_psy_fit{sess}{cc}(ch,2);
            psycho_fit{cc}(ch,3) = psycho_fit{cc}(ch,3) + temp_hist_psy_fit{sess}{cc}(ch,3);
            
            psycho_fit_1st_half{cc}(ch,2) = psycho_fit_1st_half{cc}(ch,2) + sum(sess_1st_half(this_condition_1st_half & this_heading_1st_half, 4)==RIGHT);
            psycho_fit_1st_half{cc}(ch,3) = psycho_fit_1st_half{cc}(ch,3) + sum(this_condition_1st_half & this_heading_1st_half);
            
            psycho_fit_2nd_half{cc}(ch,2) = psycho_fit_2nd_half{cc}(ch,2) + sum(sess_2nd_half(this_condition_2nd_half & this_heading_2nd_half, 4)==RIGHT);
            psycho_fit_2nd_half{cc}(ch,3) = psycho_fit_2nd_half{cc}(ch,3) + sum(this_condition_2nd_half & this_heading_2nd_half);
            
        end
        temp_hist_psy_fit{sess}{cc}(:,1)= grand_heading;
        psycho_fit{cc}(:,1) = grand_heading;
        psycho_fit_1st_half{cc}(:,1) = grand_heading;
        psycho_fit_2nd_half{cc}(:,1) = grand_heading;
    end
    % ================================================================================
    
    % ==== Confidence-dependent Psychometric curve of current trials ======
    % Divide trials according to pre condition & cur condtion &
    % pre heading & cur heading & correct pre choice
    fit_hist_psycho{sess}(:, 1) = temp_hist_psy_fit{sess};
    
    for cc = 1:length(unique_conditions)  % current condition
        cur_cond_ind = cur_trials(:, 1)==unique_conditions(cc);
        
        for choice = LEFT:RIGHT  % previous trial choice
            pre_correct_choice_ind = pre_trials(:,4)==choice & pre_trials(:,5)==0;  % only considering correct previous trials
            
            for ch = 1:length(unique_headings) % ch, current heading
                cur_heading_ind = cur_trials(:, 3)==unique_headings(ch);
                hist_ind = pre_correct_choice_ind & cur_heading_ind & cur_cond_ind;
                
                psycho_fit_hist{cc,choice+1}(ch,2) = psycho_fit_hist{cc,choice+1}(ch,2) + sum(cur_trials(hist_ind, 4)==RIGHT);
                psycho_fit_hist{cc,choice+1}(ch,3) = psycho_fit_hist{cc,choice+1}(ch,3) + sum(hist_ind);
                
                
                for pc = 1:length(unique_conditions)  % previous condition
                    pre_cond_ind = pre_trials(:, 1)==unique_conditions(pc);
                    
                    for ph = 1:length(pre_headings) % ph, previous heading
                        if pre_headings(ph) == 0
                            pre_head_correct_ind = (pre_trials(:, 3) ==pre_headings(ph) & pre_trials(:,5) == 0 & pre_trials(:,4)==(ph-4));
                        else
                            pre_head_correct_ind = (pre_trials(:, 3)==pre_headings(ph) & pre_trials(:, 5)==0);  % only consider those whose previous trials are correct
                        end
                        %                         pre_head_correct_ind = (pre_trials(:, 3)==pre_headings(ph) & pre_trials(:, 5)==0);  % only consider those whose previous trials are correct
                        %                         if ph == ceil(length(unique_headings)/2)
                        %                             pre_head_correct_ind =
                        cur_cond_pre_cond_head_correct = pre_cond_ind & pre_head_correct_ind & cur_cond_ind;
                        
                        pc_cur_trial_selected = pre_correct_choice_ind & cur_cond_ind&cur_heading_ind;  % for history dependent curve, only need previous choice
                        fit_hist_psycho{sess}{cc,choice+1}(ch,2) = sum(cur_trials(pc_cur_trial_selected,4)==RIGHT);
                        fit_hist_psycho{sess}{cc,choice+1}(ch,3) = sum(pc_cur_trial_selected);
                        
                        confidence_trial_ind = cur_heading_ind & cur_cond_pre_cond_head_correct;   % for confidence dependent curve, need previous heading & choice
                        confidence_trial_ind_1st = confidence_trial_ind(1:floor(length(confidence_trial_ind)/2));
                        confidence_trial_ind_2nd = confidence_trial_ind(length(confidence_trial_ind_1st)+1:end);
                        
                        cur_depend_on_pre{cc,pc}{ph}(ch,2) = cur_depend_on_pre{cc,pc}{ph}(ch,2) + sum(cur_trials(confidence_trial_ind,4)==RIGHT);
                        cur_depend_on_pre{cc,pc}{ph}(ch,3) = cur_depend_on_pre{cc,pc}{ph}(ch,3) + sum(confidence_trial_ind);
                        
                        cur_depend_on_pre_1st{cc,pc}{ph}(ch,2) = cur_depend_on_pre_1st{cc,pc}{ph}(ch,2) + sum(sess_1st_half(confidence_trial_ind_1st,4)==RIGHT);
                        cur_depend_on_pre_1st{cc,pc}{ph}(ch,3) = cur_depend_on_pre_1st{cc,pc}{ph}(ch,3) + sum(confidence_trial_ind_1st);
                        cur_depend_on_pre_2nd{cc,pc}{ph}(ch,2) = cur_depend_on_pre_2nd{cc,pc}{ph}(ch,2) + sum(sess_2nd_half(confidence_trial_ind_2nd,4)==RIGHT);
                        cur_depend_on_pre_2nd{cc,pc}{ph}(ch,3) = cur_depend_on_pre_2nd{cc,pc}{ph}(ch,3) + sum(confidence_trial_ind_2nd);
                        
                        cur_depend_on_pre{cc,pc}{ph}(:,1) = grand_heading;
                        cur_depend_on_pre_1st{cc,pc}{ph}(:,1) = grand_heading;
                        cur_depend_on_pre_2nd{cc,pc}{ph}(:,1) = grand_heading;
                        
                    end
                end
            end
            fit_hist_psycho{sess}{cc,choice+1}(:,1) = grand_heading;    % 3 rows, different cur conditions; 3 columns, different history: 1st grand mean, 2nd prevous left choice correct, 3rd pre right choice correct
            psycho_fit_hist{cc,choice+1}(:,1) = grand_heading;
            
        end
    end
end
psycho_fit_hist(:,1) = psycho_fit;


%% Dataset for fitting
psy_fit = PackResult(psycho_fit, psycho_fit_1st_half, psycho_fit_2nd_half);  % grand mean comparing
psycho_fit_hist;  % for grand mean, history dependent curves
fit_hist_psycho;  % for each session, history dependent curves for different current conditions (left or right choice history)
fit_hist_psycho_logistic = fit_hist_psycho;  % for logistic test for each session
confidence_fit = PackResult(psycho_fit, cur_depend_on_pre, psycho_fit_1st_half, cur_depend_on_pre_1st, psycho_fit_2nd_half, cur_depend_on_pre_2nd);  % confidence dependent curve plotting (headings and conditions history)
% confidence_1st_half = PackResult( psycho_fit_1st_half, cur_depend_on_pre_1st);  % confidence dependent curve plotting & dividing sessions into two halves
% confidence_2nd_half = PackResult( psycho_fit_2nd_half, cur_depend_on_pre_2nd);

% ------------------------------------------------------- Data fitting ----------------------------------------------
%== Grand mean comparing
psy_field = fieldnames(psy_fit);

for cc = 1: length(unique_conditions) % cc, cur condition
    for ds = 1:length(psy_field) % ds, datasets
        psy_fit.(psy_field{ds}){cc}(:,2) = psy_fit.(psy_field{ds}){cc}(:,2) ./ psy_fit.(psy_field{ds}){cc}(:,3);  %  rightward choice%
        
        [bb, tt] = cum_gaussfit_max1(psy_fit.(psy_field{ds}){cc}, method, 0);  % MLE fitting
        
        Psy_perf{cc, ds} = [bb tt];
    end
end


%== History-dependent grand mean comparing
for cc = 1:size(psycho_fit_hist,1)
    for choice = 1:size(psycho_fit_hist,2)
        psycho_fit_hist{cc,choice}(:,2) = psycho_fit_hist{cc,choice}(:,2) ./ psycho_fit_hist{cc,choice}(:,3);
        
        [bb, tt] = cum_gaussfit_max1(psycho_fit_hist{cc,choice}, method, 0);
        
        Psy_hist_grand_perf{cc,choice} = [bb tt];
    end
end



%== Previous choice dependent curve for each session
comp_pair = {[1 2], [1 3]};  % 1,grand; 2, pre left; 3, pre right
for sess = 1: length(fit_hist_psycho)
    if ~isempty(fit_hist_psycho{sess})
        for cc = 1: length(unique_conditions)
            
            %-- For logistic test
            for choice = 1: 3
                fit_hist_psycho_logistic{sess}{cc,choice}(:, 2) = fit_hist_psycho_logistic{sess}{cc,choice}(:, 2) ./ fit_hist_psycho_logistic{sess}{cc,choice}(:, 3);
                % interpolation when monkeys strangely bias to one target
                if choice ==3
                    if sum(fit_hist_psycho_logistic{sess}{cc,2}(:,2)<=0) > length(grand_heading)-3 || sum(fit_hist_psycho_logistic{sess}{cc,2}(:,2)>=1) > length(grand_heading)-3,...
                            sum(fit_hist_psycho_logistic{sess}{cc,2}(:,3)<=0) > length(grand_heading)-3 || sum(fit_hist_psycho_logistic{sess}{cc,3}(:,2)>=1) > length(grand_heading)-3
                        xxi = min(grand_heading) : 1: max(grand_heading);
                        for choice = 1: 3
                            temp{sess}{cc,choice}(:,1) = xxi';
                            temp_logistic{sess}{cc,choice}(:,2) = interp1(fit_hist_psycho_logistic{sess}{cc,choice}(:,1), fit_hist_psycho_logistic{sess}{cc,choice}(:,2), xxi', 'pchip');
                            
                            % correction of interpolation data
                            temp_logistic{sess}{cc,choice}(temp_logistic{sess}{cc,choice}(:,2)<0,2) = 0;
                            temp_logistic{sess}{cc,choice}(temp_logistic{sess}{cc,choice}(:,2)>1,2) = 1;
                            
                            temp_logistic{sess}{cc,choice}(:,3) = max(fit_hist_psycho_logistic{sess}{cc,1}(:, 3));
                            
                        end
                        
                        fit_hist_psycho_logistic{sess}(cc,:) = temp_logistic{sess}(cc,:);
                    end
                end
            end
            
            %-- Curves fitting
            for choice = 1: 3   % pc, previous choice, 2-left, 3-right
                fit_hist_psycho{sess}{cc,choice}(:, 2) =  fit_hist_psycho{sess}{cc,choice}(:, 2) ./  fit_hist_psycho{sess}{cc,choice}(:, 3);
                
                [bb, tt] = cum_gaussfit_max1( fit_hist_psycho{sess}{cc,choice});
                
                Psy_hist_perf{sess}{cc,choice} = [bb tt];
                
                % for logistic test
                % manually constructing dataset
                fit_hist_psycho_logistic{sess}{cc,choice}(:, 3) =  max(fit_hist_psycho_logistic{sess}{cc,1}(:, 3));    % manully set repetition the same as grand mean, only for logistic test
                fit_hist_psycho_logistic{sess}{cc,choice}(:, 2) = fit_hist_psycho_logistic{sess}{cc,choice}(:, 2) .* fit_hist_psycho_logistic{sess}{cc,choice}(:, 3);
            end
            
            for cp = 1:length(comp_pair)  % cp, comparing pairs
                % Constructing dataset for each session for logistic test
                yy(:,1) = [fit_hist_psycho_logistic{sess}{cc,comp_pair{cp}(1)}(:,1); fit_hist_psycho_logistic{sess}{cc,comp_pair{cp}(2)}(:,1)];
                yy(:,2) = [ones(length(yy(:,1))/2, 1); zeros(length(yy(:,1))/2, 1)];
                yy(:,3) = [fit_hist_psycho_logistic{sess}{cc,comp_pair{cp}(1)}(:,2); fit_hist_psycho_logistic{sess}{cc,comp_pair{cp}(2)}(:,2)];
                yy(:,4) = [fit_hist_psycho_logistic{sess}{cc,comp_pair{cp}(1)}(:,3); fit_hist_psycho_logistic{sess}{cc,comp_pair{cp}(2)}(:,3)];
                
                % use logit fit, probit fit might be better, need to compare sometimes
                [b, dev, stats] = glmfit([yy(:,1) yy(:,2) yy(:,1).*yy(:,2)], [yy(:,3) yy(:,4)], 'binomial', 'link', 'probit');
                
                p_value.bias{sess}(cc,cp) = stats.p(3);
                p_value.slope{sess}(cc,cp) = stats.p(4);
            end
            
        end
    end
end


% == Confidence-dependent curves and updating (Refer to Lak et al., 2020, eLife)
Difficult = 1;  Easy = 2;
% Psychometric curves fitting
fN = fieldnames(confidence_fit); %fN, fieldname
for pair = 1:3
    for cc = 1: length(psycho_fit_1st_half)   % cc, current condition
        confidence_fit.(fN{pair*2-1}){cc}(:,2) = confidence_fit.(fN{pair*2-1}){cc}(:,2) ./ confidence_fit.(fN{pair*2-1}){cc}(:,3);  % change confidence_fit but reserve raw psycho_fit
        for pc = 1:size(cur_depend_on_pre,2)  % pc, previous condition
            for ph = 1:length(pre_headings)
                % Fitting
                confidence_fit.(fN{pair*2}){cc,pc}{ph}(:,2) = confidence_fit.(fN{pair*2}){cc,pc}{ph}(:,2) ./ confidence_fit.(fN{pair*2}){cc,pc}{ph}(:,3);
                [bb, tt] = cum_gaussfit_max1(confidence_fit.(fN{pair*2}){cc,pc}{ph}, method, 0);
                confidn_psy_perf{pair}{cc,pc}{ph} = [bb tt];
                
                % Confidence-dependent choice updating
                choice_updating{pair}{cc,pc}(:,ph) = confidence_fit.(fN{pair*2}){cc,pc}{ph}(:,2) - confidence_fit.(fN{pair*2-1}){cc}(:,2); % row, cur heading, column, pre heading
            end
            
            % Choice updating as a function of previous headings for current easy and
            % difficult trials
            ph_depen_cu{pair}{cc,pc}(Difficult, :) = nanmean(choice_updating{pair}{cc,pc}(3:7,:), 1);
            ph_depen_cu{pair}{cc,pc}(Easy, :) = nanmean(choice_updating{pair}{cc,pc}(setdiff(1:9, 3:7),:), 1);
            
        end
    end
end


%% ================================ Plotting =========================================================
%==== Dataset for plotting =====
Psy_perf; % Grand mean comparing, 3 stim conditions * 3 grand levels (all trials, first half of a session, second half of a session )
Psy_hist_perf; % Curves fitting for each session
p_value; % p_values for bias and slope of Psy_hist_perf for each session
Psy_hist_grand_perf; % History-dependent grand mean comparing, 3 stim conditions * 3 pre_trial conditions (grand, pre left choice, pre right choice)
confidn_psy_perf; % psychometric curves depend on previous conditions and headings  % can be compared with Psy_perf
choice_updating; % confidence-dependent choice updating
ph_depen_cu; % previous-trial difficulty dependent choice updating

outer_color = {'b'; 'r'; 'g'; 'k'};  % outer loop, ves-b, vis-r, comb-g, grand mean-k
titleX = {'Vest', 'Vis', 'Comb'};
titleY = titleX;
xxi = min(grand_heading) : 0.05 : max(grand_heading);
pre_headings = [grand_heading(1:5) 0 grand_heading(6:end)];

%  -------------  Grand mean comparing --------------------------
% As negtive control
curve_color = [0 0 0.5; 0.5 0 0; 0 0.5 0];  % vest; vis; comb;
symbo =  ['o'; '^';'v'];
lgd = {'ALL'; '1st half'; '2nd half'};

set(figure(1921), 'name', 'Grand Mean Comparing', 'pos', [27 63 1800 500]); clf(1921);
h = tight_subplot(1, length(unique_conditions), [.1 .03], [.1 .2]);

for cc = 1: length(unique_conditions)
    textN = 0; % text index
    set(gcf, 'CurrentAxes', h(cc));
    for ds = 1:3 % different datasets
        plot(xxi, cum_gaussfit(Psy_perf{cc, ds}, xxi),'Color',curve_color(cc,:)*(ds-1),'Linewidth',3,'MarkerFaceColor',curve_color(cc,:)*(ds-1));
        
        textN = textN +1;
        text(max(grand_heading)*0.15, 0.25-textN*0.07, sprintf('%5.2f   %5.2f', Psy_perf{cc, ds}(1), Psy_perf{cc, ds}(2)),'color', curve_color(cc,:)*(ds-1), 'fontsize', 12);
        
        hold on;
        xlabel('Heading Angles');
        ylim([0, 1]);
        
    end
    box off;
    legend(lgd, 'Location','northwest','AutoUpdate','off');
    if cc == 1
        ylabel('Proportion of Rightward Choice')
    end
    
    for ds = 1:3
        plot(grand_heading, psy_fit.(psy_field{ds}){cc}(:,2), symbo(ds), 'Color', curve_color(cc,:)*(ds-1), 'MarkerSize', 9, 'MarkerFaceColor', curve_color(cc,:)*(ds-1));
    end
    
    plot([min(grand_heading) max(grand_heading)], [0.5 0.5], 'k--');
    plot([0 0], [0 1], 'k--');
    
    text(max(grand_heading)*0.15, 0.25, 'bias   threshold', 'color','k', 'fontsize', 12);
    title(titleX{cc}, 'Color', outer_color{cc});
    
end
SetFigure();



%------------------------- Bias and threshold develop across sessions -------------
empty_ind = cellfun(@(x) isempty(x), Psy_hist_perf);
temp_Psy_hist_perf = Psy_hist_perf(~empty_ind);
cell_elem = cellfun(@(x) cell2mat(x), temp_Psy_hist_perf, 'UniformOutput', false)';
Psy_hist_perf_2plot = cell2mat(cell_elem);
sessN2plot = length(Psy_hist_perf_2plot)/length(unique_conditions);

% As negtive control to indicate stable behavior of monkey
Bias = 1; Threshold = 2;

set(figure(1931), 'name', ' Bias and threshold develop across sessions', 'pos', [27 63 1800 500]); clf(1931);
h = tight_subplot(1,2, [.1 .1], [.1 .2]);

for option = Bias : Threshold
    set(gcf, 'CurrentAxes', h(option));
    
    for cc = 1:length(unique_conditions)
        plot(1:sessN2plot, Psy_hist_perf_2plot(cc:3:end, option), 'Color',outer_color{cc},'Linewidth', 2); hold on;
    end
    legend(titleX, 'Location', 'northeast', 'AutoUpdate','off');
    
    if option == Bias
        axis tight;
        ylims = ylim; ylims = max(abs(ylims));
        ylim(1.05*[-ylims ylims]);
        xlabel('Sessions'); ylabel('Bias');
    end
    
    if option == Threshold
        axis tight;
        ylims = ylim; ylims = max(ylims);
        ylim([0 1.05*ylims]);
        xlabel('Sessions'); ylabel('Threshold');
    end
    
end
SetFigure();


%--------------------- Grand history curves ----------------------------
curve_color = [0 0 0.5; 0.5 0 0; 0 0.5 0];  % vest; vis; comb;
symbo =  ['o'; '^';'v'];
lgd = {'ALL'; 'Pre Left Choice'; 'Pre Right Choice'};

set(figure(1941), 'name', 'Grand history-dependent curves', 'pos', [27 63 1800 500]); clf(1941);
h = tight_subplot(1, length(unique_conditions), [.1 .03], [.1 .2]);

for cc = 1: length(unique_conditions)
    textN = 0; % text index
    set(gcf, 'CurrentAxes', h(cc));
    for pc = 1:3 % different pre-trial history
        plot(xxi, cum_gaussfit(Psy_hist_grand_perf{cc, pc}, xxi),'Color',curve_color(cc,:)*(pc-1),'Linewidth',3,'MarkerFaceColor',curve_color(cc,:)*(pc-1));
        
        textN = textN +1;
        text(max(grand_heading)*0.15, 0.25-textN*0.07, sprintf('%5.2f   %5.2f', Psy_hist_grand_perf{cc, pc}(1), Psy_hist_grand_perf{cc, pc}(2)),'color', curve_color(cc,:)*(pc-1), 'fontsize', 12);
        
        hold on;
        xlabel('Heading Angles');
        ylim([0, 1]);
        
    end
    box off;
    legend(lgd, 'Location','northwest','AutoUpdate','off');
    if cc == 1
        ylabel('Proportion of Rightward Choice')
    end
    
    for pc = 1:3
        plot(grand_heading, psycho_fit_hist{cc,pc}(:,2), symbo(pc), 'Color', curve_color(cc,:)*(pc-1), 'MarkerSize', 9, 'MarkerFaceColor', curve_color(cc,:)*(pc-1));
    end
    
    plot([min(grand_heading) max(grand_heading)], [0.5 0.5], 'k--');
    plot([0 0], [0 1], 'k--');
    
    text(max(grand_heading)*0.15, 0.25, 'bias   threshold', 'color','k', 'fontsize', 12);
    title(titleX{cc}, 'Color', outer_color{cc});
    
end
SetFigure();



%------- Session-by-session bias and threshold comparision ------
% == For different pre-trial history
comp_pair;
Psy_hist_perf_2plot;
temp_p_bias = p_value.bias(~empty_ind);  temp_p_thres = p_value.slope(~empty_ind);
p_bias = cell2mat(temp_p_bias'); p_thres = cell2mat(temp_p_thres');

set(figure(1951), 'name', 'Session-by-Session history-dependent performance comparison', 'pos', [27 63 1200 900]); clf(1951);
h = tight_subplot(2, length(comp_pair), [.1 .1], [.1 .1]);

for option = Bias : Threshold
    maxi_option(option) = max(max(abs(Psy_hist_perf_2plot(:,[1 2 3]*option))));
    for cp = 1: length(comp_pair)
        set(gcf, 'CurrentAxes', h(option+2*(cp-1)))
        %         subplot(2,2,option+2*(cp-1));
        for cc = 1: length(unique_conditions)
            temp_compare = Psy_hist_perf_2plot(cc:3:end, option+2*(comp_pair{cp}-1));
            
            [hh, pp]= ttest(temp_compare(:,1), temp_compare(:,2));   % paired t-test
%             [~, ~, pp]= crosstab(temp_compare(:,1), temp_compare(:,2));   % chi-square test 

            hypoth{option}(cc, cp) = hh;
            population_p{option}(cc, cp) = pp;
            
            % Mark sessions with significant difference
            if option == Bias
                sig_ind = p_bias(cc:3:end,cp) < 0.05;
                temp_compare_sig = temp_compare(sig_ind, :);
            else
                sig_ind = p_thres(cc:3:end,cp) < 0.05;
                temp_compare_sig = temp_compare(sig_ind, :);
            end
            
            plot(temp_compare(:, 1), temp_compare(:,2), 'o', 'color', outer_color{cc},'MarkerSize',8, 'LineWidth', 0.75); hold on;
            plot(temp_compare_sig(:,1), temp_compare_sig(:,2), 'o', 'color', outer_color{cc},'MarkerSize',8,'MarkerFaceColor',outer_color{cc}, 'LineWidth', 0.75);
        end
        
        if option == Bias
            axis([-maxi_option(option) maxi_option(option) -maxi_option(option) maxi_option(option)]*1.05);
            xlims = xlim; ylims =ylim;
            plot([xlims(1) xlims(2)], [ylims(1) ylims(2)], 'k-');
            xticks(ceil(xlims(1)):2:floor(xlims(2)));  yticks (ceil(ylims(1)):2:floor(ylims(2)));
            if cp == 1
                xlabel('Bias of All'); ylabel('Bias of Previous Left Choice');
            else
                xlabel('Bias of All'); ylabel('Bias of Previous Right Choice');
            end
        end
        
        if option == Threshold
            axis([0 maxi_option(option) 0 maxi_option(option)]*1.05);
            xlims = xlim; ylims =ylim;
            plot([xlims(1) xlims(2)], [ylims(1) ylims(2)], 'k-');
            xticks(0:1:floor(xlims(2)));  yticks (0:1:floor(ylims(2)));
            if cp == 1
                xlabel('Threshold of All'); ylabel('Threshold of Previous Left Choice');
            else
                xlabel('Threshold of All'); ylabel('Threshold of Previous Right Choice');
            end
        end
        
        axis square;
        
    end
end
SetFigure();
set(findall(gcf,'type','axes','-not','tag','legend'),'box','on'); % Box on



%------------ Population bias and threshold comparison --------------
% For different pre-trial history
comp_pair;
Psy_hist_perf_2plot;
population_p;
pc_color = {'k', 'c', 'm'};  % all, pre-left, pre-right
lgd = {'ALL'; 'Pre Left Choice'; 'Pre Right Choice'};

set(figure(1961), 'name', 'Population history-dependent performance comparison', 'pos', [27 63 1800 800]); clf(1961);
h = tight_subplot(2,3,[.1 .03], [.1 .1]);

for option = Bias : Threshold
    if option == Bias
        edges = linspace(-maxi_option(option),maxi_option(option),10);
    else
        edges = linspace(0, maxi_option(option), 10);
    end
    
    for cc = 1: length(unique_conditions)
        set(gcf, 'CurrentAxes', h(option+2*(cc-1)));
        
        for pc = 1:3 % previous choice
            popu2plot{option,cc}(:, pc) = Psy_hist_perf_2plot(cc:3:end, option+(pc-1)*2);
            
            histogram(popu2plot{option,cc}(:, pc), edges, 'EdgeColor', pc_color{pc}, 'DisplayStyle', 'stairs', 'LineWidth',2, 'Normalization', 'probability');
            hold on;
        end
        legend(lgd,'Location','northeast', 'AutoUpdate', 'off');
        ylim([0 0.8]);
        
        popu_mean{option, cc} = mean(popu2plot{option, cc});
%         popu_mean{option, cc} = median(popu2plot{option, cc});

        for pc = 1:3
            plot([popu_mean{option, cc}(1,pc) popu_mean{option, cc}(1,pc)], [0 0.8], 'color', pc_color{pc}, 'LineStyle', '--', 'Linewidth', 2.5);
        end
        
        for cp = 1:length(comp_pair)
            xlims = xlim;
            %             text(xlims(2)-1, 0.5, 'p value', 'color', 'k');
            text(xlims(2), 0.5-0.07*cp, [num2str(comp_pair{cp}(1)) ' vs. ' num2str(comp_pair{cp}(2)) ' :  ' sprintf('%0.3e', population_p{option}(cc,cp))],...
                'color', outer_color{cc}, 'HorizontalAlignment', 'right');
        end
        
        if option == Bias
            title(titleX{cc}, 'color', outer_color{cc});
        end
        
        if cc == 2
            if option == Bias
                xlabel('Bias');
            else
                xlabel('Threshold');
            end
        end
        
    end
end
SetFigure();



%-------------------- Confidence-dependent psycho curves ----------------
confidn_psy_perf;
psy_fit; psy_field;

lgd = {'ALL'; '1st half'; '2nd half'};
heading_color = colormap(winter);
inner_color = heading_color(end-round(linspace(1,64,length(pre_headings)))+1, :);    % inner loop, different heading history
pre_heading = {'-10', '-5', '-2.5', '-1.25', '-0', '+0', '1.25','2.5','5','10'};
pre_heading_str = cellfun(@(x) num2str(x), pre_heading, 'UniformOutput', false);

for ds = 1:3 % all, 1st half , 2nd half
    set(figure(1961+ds*10),'name',['Confidence-dependent Psychometric Curves --' lgd{ds}], 'pos',[27 63 1500 900]); clf(1961+ds*10);
    h = tight_subplot(3,3,[.03 .03], [.05 .1], [.1 .1]);
    
    for cc = 1:length(unique_conditions)
        
        for pc = 1:length(unique_conditions)
            set(gcf, 'CurrentAxes', h(cc+3*(pc-1)));
            
            for ph = 1:length(pre_headings)
                plot(xxi, cum_gaussfit(confidn_psy_perf{ds}{cc,pc}{ph}, xxi), 'color', inner_color(ph,:), 'Linewidth', 2); hold on;
                
                ylim([0, 1]);
                xlim([-11, 11]); xticks(min(grand_heading):2:max(grand_heading));
                
            end
            
            plot(xxi, cum_gaussfit(Psy_perf{cc,1}, xxi), 'color', 'k', 'Linewidth', 2, 'DisplayName', 'Grand Mean');
            plot([0 0], [0 1], 'k--');
            plot(xlim, [0.5 0.5], 'k--');
            
            if cc+3*(pc-1) == 7
                legend(pre_heading_str, 'Location', 'southeast', 'AutoUpdate', 'off');
            end
            
            plot(psy_fit.(psy_field{ds}){cc}(:,1), psy_fit.(psy_field{ds}){cc}(:,2), 'k.', 'MarkerSize', 10);
            
            for ph = 1: length(grand_heading)
                plot(cur_depend_on_pre{cc,pc}{ph}(:,1), cur_depend_on_pre{cc,pc}{ph}(:,2)./cur_depend_on_pre{cc,pc}{ph}(:,3),...
                    '.', 'color', inner_color(ph,:), 'MarkerSize', 10);
            end
            
            if pc == 1
                YaxisLabel = sprintf('Current Condition: %s', titleY{cc});
                ylabel({YaxisLabel; 'Rightward Rate'}, 'color', outer_color{4}, 'FontWeight', 'bold', 'FontSize', 12);
            end
            
            if cc == 1
                title(['Previous Condition: ', titleX{pc}],'color',outer_color{pc},'FontWeight','bold','FontSize',15)
            end
            
            if cc == 3
                xlabel('Current Heading','FontWeight','bold','FontSize',12);
            else
                xticklabels({});
            end
        end
    end
    SetFigure();
end



% ----------------- Confidence-dependent choice updating --------------------
choice_updating;
% using similar color as Lak paper
white = [255 255 255]/255; % for 0 updating
purple = [160 32 240]/255; % for positive updating
orange = [255 128 0]/255; % for negtive updating
color_data_neg = flipud((0:49)'/50 *(orange-white)+white);  color_data_posi = (0:49)'/50*(purple-white)+white;
color_map = [color_data_neg;color_data_posi];
cur_heading = grand_heading;

for ds = 1:3
    set(figure(1991+ds*10),'name',['Confidence-dependent Choice Updating --' lgd{ds}], 'pos',[27 63 1500 900]); clf(1991+ds*10);
    h = tight_subplot(3,3, [.03 .03]);
    
    temp_maxi_updating = cellfun(@(x) max(max(abs(x))), choice_updating{ds});
    maxi_updating{ds} = max(max(temp_maxi_updating));
    
    for cc = 1:length(unique_conditions)
        for pc = 1:length(unique_conditions)
            set(gcf, 'CurrentAxes', h(cc+(pc-1)*3));
            
            scale= choice_updating{ds}{cc, pc}/maxi_updating{ds};
            sign_ind = sign(choice_updating{ds}{cc, pc});
            %           color_data_0 = scale(sign_ind==0)* white;
            
            hm = heatmap(gcf,pre_heading,cur_heading,choice_updating{ds}{cc, pc}, 'colormap', color_map,...
                'ColorLimits', [-maxi_updating{ds} maxi_updating{ds}], 'ColorbarVisible', 'off', 'GridVisible', 'off', 'CellLabelColor','none');
            hm.FontSize = 13;
            
            if pc == 1
                hm.YLabel = {['Current Condition: ', titleY{cc}] ;'Current Heading'};
            end
            
            if cc ==1
                hm.Title = ['Previous Condition: ', titleX{pc}];
                if pc == 3
                    hm.ColorbarVisible = 'on';
                end
            end
            
            if cc == 3
                hm.XLabel = 'Previous Heading' ;
            end
            
        end
    end
end



% ---------------- Confidence-dependent choice updating curves --------------
ph_depen_cu;
dot_color = {'k', [0.75 0.75 0.75]}; % dots color
line_color = {'k', 'r'};  % k, pre-left heading, gray, pre-right heading
markers = {'o','s'};
lgd = {'Hard', 'Easy','Left', 'Right'};
nam = {'ALL'; '1st half'; '2nd half'};

for ds = 1:3
    set(figure(2021+ds*10),'name',['Confidence-dependent choice updating curves --', nam{ds}],'pos',[27 63 1500 900]); clf(2021+ds*10);
    h= tight_subplot(3,3,[.05 .05],[.05 .1], [.1 .05]);
    
    temp_ph_maxi = cellfun(@(x) max(max(abs(x))), ph_depen_cu{ds});
    ph_maxi(ds) = max(max(temp_ph_maxi));
    
    for cc = 1: length(unique_conditions)
        for pc = 1: length(unique_conditions)
            set(gcf, 'CurrentAxes', h(cc+(pc-1)*3));
            
            for i = Difficult : Easy
                plot(pre_headings, ph_depen_cu{ds}{cc,pc}(i,:),...
                    markers{i}, 'color', dot_color{i}, 'MarkerSize', 10, 'MarkerFaceColor', dot_color{i});
                hold on;
            end
            
            for i = Difficult : Easy
                plot(pre_headings(1:length(pre_headings)/2),ph_depen_cu{ds}{cc,pc}(i,1:length(pre_headings)/2), line_color{1}, 'LineWidth', 2);
                plot(pre_headings(length(pre_headings)/2+1:end),ph_depen_cu{ds}{cc,pc}(i,length(pre_headings)/2+1:end), line_color{2}, 'LineWidth', 2);
            end
            ylim([-ph_maxi(ds) ph_maxi(ds)]);
            xlim([-11, 11]); xticks(min(grand_heading):2:max(grand_heading));
            
            if cc ==1
                title(['Previous Condition: ', titleX{pc}]);
                if pc == 3
                    legend(lgd, 'Location', 'northeast', 'AutoUpdate', 'off');
                end
            end
            
            if pc == 1
                ylabel({['Current Condition: ', titleY{cc}] ;'Updating'});
            end
            
            if cc == 3
                xlabel('Previous Heading');
            end
        end
    end
    SetFigure();
end



%% Saving
saveas(1921, 'Grand Mean Comparing Psycho', 'svg'); saveas(1921, 'Grand Mean Comparing Psycho', 'jpg'); 
saveas(1931, 'Bias and Threshold Developing', 'svg'); saveas(1931, 'Bias and Threshold Developing', 'jpg'); 
saveas(1941, 'Grand History-dependent Psycho', 'svg'); saveas(1941, 'Grand History-dependent Psycho', 'jpg'); 
saveas(1951, 'Sessions History-dependent Perf Comparing', 'svg'); saveas(1951, 'Sessions History-dependent Perf Comparing', 'jpg'); 
saveas(1961, 'Population History-dependent Perf Comparing', 'svg'); saveas(1961, 'Population History-dependent Perf Comparing', 'jpg'); 
saveas(1971, 'Confidence-dependent Psycho_ALL', 'svg'); saveas(1971, 'Confidence-dependent Psycho_ALL', 'jpg'); 
saveas(1981, 'Confidence-dependent Psycho_1st', 'svg'); saveas(1981, 'Confidence-dependent Psycho_1st', 'jpg'); 
saveas(1991, 'Confidence-dependent Psycho_2nd', 'svg'); saveas(1991, 'Confidence-dependent Psycho_2nd', 'jpg'); 
saveas(2001, 'Confidence-dependent Choice Updating_ALL', 'svg'); saveas(2001, 'Confidence-dependent Choice Updating_ALL', 'jpg'); 
saveas(2011, 'Confidence-dependent Choice Updating_1st', 'svg'); saveas(2011, 'Confidence-dependent Choice Updating_1st', 'jpg'); 
saveas(2021, 'Confidence-dependent Choice Updating_2nd', 'svg'); saveas(2021, 'Confidence-dependent Choice Updating_2nd', 'jpg'); 
saveas(2031, 'Confidence-dependent Choice Updating Curves_ALL', 'svg'); saveas(2031, 'Confidence-dependent Choice Updating Curves_ALL', 'jpg'); 
saveas(2041, 'Confidence-dependent Choice Updating Curves_1st', 'svg'); saveas(2041, 'Confidence-dependent Choice Updating Curves_1st', 'jpg'); 
saveas(2051, 'Confidence-dependent Choice Updating Curves_2nd', 'svg'); saveas(2051, 'Confidence-dependent Choice Updating Curves_2nd', 'jpg'); 

end