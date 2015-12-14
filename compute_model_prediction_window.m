function [prob1_model prob2_model prob1_actl prob2_actl] = compute_model_prediction_window(padded, plot_mode, window_size)
% COMPUTE_MODEL_PREDICTION_WINDOW computes a moving window model prediction
% [note] currently only works for RLbeta_alt4_v6
cd \projects\SocialInflu\sit_stan\

load '_data\data_129.mat'
load '_outputs\param_beta4_alt6'

%% some checking
if nargin < 1
  error('No data');
elseif nargin < 2
    plot_mode = 0;
    window_size = 10;
elseif nargin < 3
    % number of trials per bin, possibly adjust this
    % (it might be quite noisy per subject anyway!)
    window_size = 10;
end

prob1_model = nan(size(data,1), 100);
prob2_model = nan(size(data,1), 100);
prob1_actl  = nan(size(data,1), 100);
prob2_actl  = nan(size(data,1), 100);

%% compute the moving-window model prediction

for s = 1:size(data,1) % subjects
    
    subData = data{s};
    [~,~,~,model] = RevLearn_RLbeta4_alt6(parMat(s,:), subData);
    
    % extract choice
    choice1 = model.c1; % 1st choice, 1 or 2
    choice1(choice1==2) = 0;
    choice2 = model.sw; % choice switch, 0 or 1
    
    % extract action probs of each trial
    prob1 = model.prob1(:,1); 
    prob2 = model.prob2;
    
    % sort probability
    [prob1_sorted, ind1] = sort(prob1); 
    [prob2_sorted, ind2] = sort(prob2); 
    
    % sort choice
    choice1_sorted = choice1(ind1);
    choice2_sorted = choice2(ind2);
    
    if size(prob1,1) == 1
        my_window = ones(1,window_size);  % row_vector
    else
        my_window = ones(window_size,1);  % vector
    end
    
    if padded == 0        
        p_choice1_sorted = conv(choice1_sorted,my_window,'same');
        p_choice2_sorted = conv(choice2_sorted,my_window,'same');
        
    elseif padded == 1  % try this if the edges are noisy (the result is better)
        
        if size(choice1_sorted,1) == 1
            choice1_padded = [zeros(1,floor(window_size/2)) choice1_sorted
                ones(1,floor(window_size/2))];
            choice2_padded = [zeros(1,floor(window_size/2)) choice2_sorted
                ones(1,floor(window_size/2))];
        else
            choice1_padded = [zeros(floor(window_size/2),1); choice1_sorted;
                ones(floor(window_size/2),1)];
            choice2_padded = [zeros(floor(window_size/2),1); choice2_sorted;
                ones(floor(window_size/2),1)];
        end
        
        p_choice1_sorted = conv(choice1_padded,my_window,'same');
        p_choice1_sorted(1:floor(window_size/2)) = [];
        p_choice1_sorted(end-floor(window_size/2)+1:end) = [];
        
        p_choice2_sorted = conv(choice2_padded,my_window,'same');
        p_choice2_sorted(1:floor(window_size/2)) = [];
        p_choice2_sorted(end-floor(window_size/2)+1:end) = [];
        
    end
    
    prob1_model(s,:) = prob1_sorted;
    prob2_model(s,:) = prob2_sorted;
    prob1_actl(s,:)  = p_choice1_sorted / window_size;
    prob2_actl(s,:)  = p_choice2_sorted / window_size;
    
end

if plot_mode
   figure
   
   %%% --- model prediction plot of choice1
   subplot(1,2,1)
   
   
   %%% --- model prediction plot of choice2, switch (0) or stay (1)
   subplot(1,2,2)
   
    
    
    
    
    
    
    
    
    
    
    
    
    
end













