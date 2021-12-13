function [params, params_no_lapse, Dev, pDev,DevSim,convergedGoF] = cum_gauss_lapse_comparison(dataset, batch_flag) 

% This script aims to fit psychometric function using cumulative gasussian
% function and compares those models that with or without lapse assumption
% and output four fitted paramters 
% The toolkit I used was developed by Prins & Kingdom, 2018,
% naming Palamedes
% ZZ 20210914
% M1, model 1, full model where lapses are variable
% M2, model 2, restricted lapses, lambda + gamma < 0.1
% M3, model 3, fixed lapses, fixed lambda and gamma
% M4, model 4, no lapses, lambda = gamma = 0

% == Input
% dataset: n*3 matrix, the first column is independent term, different
% headings in my task; the second column is correct choice number in different
% heading conditions, rightward choice number in my task; the last column is
% total trial number of each heading condition

% == Output
% params: 4 fitted parameters (mu, sigma, lambda, gamma)
% mu: bias of psychometric function
% sigma: slope 
% gamma: left asymptote, zero if unnecessary
% lambda: right asymptote, zero if unnecessary

% Dev: Deviance of the fitted model (compared with the saturated model)
% pDev: Population of simulations whose deviance are larger than ture Dev,
           % larger pDev indicate better fitting of our models
% DevSim: Deviance of each simulation
% convergedGoF: Whether each simulation converge or not

%% Data rearrangement
if iscell(dataset)
    dataset = cell2mat(dataset);
end

if size(dataset, 2) > size(dataset, 1)
    dataset = dataset' ;
end

headings = dataset(:, 1);
trial_num = dataset(:, 3);

rightward_num = dataset(:,2);
if sum(rightward_num > 1.1) <1
    if all(rightward_num==0)
        warning('Check whether monkey always choose one target...');
    end
    rightward_num = rightward_num .* trial_num; 
end
 

%% Predefinition
% Cumulative gaussian function
PF = @PAL_CumulativeNormal;

% The ways lapse rates are handled
% nAPLE: jointly fits thresholds, slopes and lapse rates
% jAPLE is identical to nAPLE except that it assumes that the higest
% stimulus intensity is at an asymptotic level and thus that an error
% observed at this intensity can only be due to lapses
lapseFit = 'nAPLE';   
% lapseFit = 'jAPLE';

% Nelder-Mead search options
% just like fminsearch
options = PAL_minimize('options');
options.TolX = 1e-09;
options.TolFun = 1e-09;
options.MaxIter = 10000;
options.MaxFunEvals = 10000;

% Parameter search range to search for initial values
paramsFree{1} = [1 1 1 1];  % 1: free, 0: fixed
numParams{1} = sum(paramsFree{1});
% searchGrid(1).alpha = [min(headings):.05:max(headings)];
searchGrid(1).alpha = [-5:.1:5];  % bias of nearly all of sessions fall into this range
searchGrid(1).beta = 10.^ [-1:.05:1];     % inverse of threshold
searchGrid(1).gamma = [0:.01:1];
searchGrid(1).lambda =[0:.01:1]; 

% paramsFreeM2 = [1 1 1 1];  % 1: free, 0: fixed
% numParamsM2 = sum(paramsFreeM2);
% searchGridM2.alpha = [min(headings):.05:max(headings)];
% searchGridM2.beta = [-1:.05:1] .^ 10;     % inverse of threshold
% searchGridM2.gamma = [0:.005:.06];
% searchGridM2.lambda = [0:.005:.06];

paramsFree{2} = [1 1 0 1];  % 1: free, 0: fixed
numParams{2} = sum(paramsFree{2});
% searchGrid(2).alpha = [min(headings):.05:max(headings)];
searchGrid(2).alpha = [-5:.1:5];  % bias of nearly all of sessions fall into this range
searchGrid(2).beta = 10.^ [-1:.05:1] ;     % inverse of threshold
searchGrid(2).gamma = 0;
searchGrid(2).lambda = [0:.01:1];

paramsFree{3} = [1 1 1 0];  % 1: free, 0: fixed
numParams{3} = sum(paramsFree{3});
searchGrid(3).alpha = [-5:.1:5];  % bias of nearly all of sessions fall into this range
% searchGrid(3).alpha = [min(headings):.05:max(headings)];
searchGrid(3).beta = 10.^ [-1:.05:1];     % inverse of threshold
searchGrid(3).gamma = [0:.01:1];
searchGrid(3).lambda = 0;

paramsFree{4} = [1 1 0 0];
numParams{4} = sum(paramsFree{4});
% searchGrid(4).alpha = [min(headings):.05:max(headings)];
searchGrid(4).alpha = [-5:.1:5];  % bias of nearly all of sessions fall into this range
searchGrid(4).beta = 10.^ [-1:.05:1] ;     % inverse of threshold
searchGrid(4).gamma = 0;
searchGrid(4).lambda = 0;


[LLSaturated, numParamsSaturated] = PAL_PFML_LLsaturated(rightward_num, trial_num);

[paramsFitted{1},LL{1},exitflag(1)] = PAL_PFML_Fit(headings, rightward_num, trial_num, searchGrid(1), paramsFree{1}, PF,...
    'guessLimits',[0 1],'lapseLimits', [0 1], 'searchOptions', options,'lapseFit',lapseFit); 

% [paramsFittedM2,LLM2,exitflagM2] = PAL_PFML_Fit(headings, rightward_num, trial_num, searchGridM2, paramsFreeM2, PF,...
%      'guessLimits',[0 .06],'lapseLimits', [0 .06],'searchOptions', options,'lapseFit',lapseFit,'gammaEQlambda',1); 

 [paramsFitted{2},LL{2},exitflag(2)] = PAL_PFML_Fit(headings, rightward_num, trial_num, searchGrid(2), paramsFree{2}, PF,...
    'lapseLimits', [0 1], 'searchOptions', options,'lapseFit',lapseFit); 

[paramsFitted{3},LL{3},exitflag(3)] = PAL_PFML_Fit(headings, rightward_num, trial_num, searchGrid(3), paramsFree{3}, PF,...
    'guessLimits',[0 1], 'searchOptions', options,'lapseFit',lapseFit); 

[paramsFitted{4},LL{4},exitflag(4)] = PAL_PFML_Fit(headings, rightward_num, trial_num, searchGrid(4), paramsFree{4}, PF,...
    'searchOptions', options,'lapseFit',lapseFit); 

params_no_lapse = paramsFitted{4}; 


AIC = nan(1,4);
exit_ind = find(exitflag);
if ~isempty(exit_ind)
    for i= exit_ind
        AIC(i) = -2 * LL{i} + 2*numParams{i};
    end
    
    if numel(exit_ind) > 1
        [~,ind] = sort(AIC(exit_ind));        %  choose model with min(AIC)
        params = paramsFitted{exit_ind(ind(1))};
        paramsFree_GOF = paramsFree{exit_ind(ind(1))};
        searchGrid_GOF = searchGrid(exit_ind(ind(1)));

    else
        params = paramsFitted(exit_ind);
        paramsFree_GOF = paramsFree{exit_ind};
        searchGrid_GOF = searchGrid(exit_ind);

    end  
else
    params = parmasFitted{4};
    paramsFree_GOF = paramsFree{4};
    searchGrid_GOF = searchGrid(4);
end

if nargout > 2 
    %Determine Goodness-of-Fit of Lesser model
    numMCsimuls = 1000;   % number of Monte Carlo simulations
    [Dev, pMCGOF, DevSim, convergedGoF] = PAL_PFML_GoodnessOfFit(headings',rightward_num',trial_num',...
        params,paramsFree_GOF,numMCsimuls,PF,'searchGrid', searchGrid_GOF);
    
    pDev = pMCGOF; % the greater this value is, the better the model-fitting is
    
end


% %% Bootstrap simulation
% B = 1000;  % 1000 simulations
% 
% if batch_flag
%     % Parametric bootstrap to obtain SDs for parameters
%     [SD, paramsSim, LLSim, converged] = PAL_PFML_BootstrapParametric(headings,trial_num,...
%         paramsFitted, paramsFree, B, PF, 'lapseLimits',[0 1], 'searchGrid', searchGridM1, 'lapseFit', lapseFit);
%     
%     % Goodness of fit by Monte Carlo
%     [Dev, pDev, DevSim, converged] = PAL_PFML_GoodnessOfFit(headings, rightward_num, trial_num,...
%         paramsFitted, paramsFree, B, PF, 'searchGrid', searchGridM1, 'lapseLimits', [0 1], 'lapseFit', lapseFit);
% end
% 

    






