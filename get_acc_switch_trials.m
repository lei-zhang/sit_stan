function [acc_sw, acc_sw_nh] = get_acc_switch_trials
% obtain predictive accuracy only on switch trials
cd \projects\SocialInflu\sit_stan\

load '_data\data_129.mat'
load '_data\data3_129.mat'
load '_outputs\param_beta4_alt6'
load '_outputs\param_beta4_alt6_nh.mat'
load '_outputs\loo_RLbeta4_a6.mat'
load '_outputs\ppc_RLbeta4_a6.mat'
load '_outputs\ppc_RLbeta4_a6_pointwise.mat'
load '_outputs\ppc_RLbeta4_a6_pointwise_nh.mat'

ns = size(data3,3);
acc_sw = ones(ns,1);
acc_sw_nh = ones(ns,1);

for s = 1:ns
   subdata = squeeze(data3(:,:,s));
   sw_ind = find(subdata(:,5)==1);
   acc_sw(s,1) = mean(ppc_RLbeta4_a6_pointwise(s, sw_ind));
   acc_sw_nh(s,1) = mean(ppc_RLbeta4_a6_pointwise_nh(s, sw_ind));
end

acc_sw( isnan(acc_sw) ) = 0;
acc_sw_nh( isnan(acc_sw_nh) ) = 0;

% 111's nh fit is better than hierarchical fit
% mean diff = [ 0.2883    0.3200]






