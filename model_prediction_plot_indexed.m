% plot the binned model prediction, for those whose LOOICs are low, medium, and high

cd \projects\SocialInflu\sit_stan\

load '_data\data_129.mat'
load '_data\data3_129.mat'
load '_outputs\param_beta4_alt6'
load '_outputs\loo_RLbeta4_a6.mat'
load '_outputs\ppc_RLbeta4_a6.mat'
load '_outputs\ppc_RLbeta4_a6_switchTrial.mat'
load '_outputs\missingCnt.mat'
load '_outputs\switchCnt.mat'


% obtain the 1/3 and 2/3 quantiles
qtl =  quantile(loo_RLbeta4_a6,[1/3,2/3]); 

ind_l = find( loo_RLbeta4_a6 <= qtl(1) );
ind_m = find( loo_RLbeta4_a6 > qtl(1) & loo_RLbeta4_a6 <= qtl(2) );
ind_h = find( loo_RLbeta4_a6 > qtl(2)) ;

%% plot for all the subjects
[fit1, fit2] = compute_choice_prediction_RLbeta4v6(data, parMat);

f1 = figure; set(f1,'color',[1 1 1]);
subplot(1,2,1); plot_model_prediction(fit1, 'choice1');
subplot(1,2,2); plot_model_prediction(fit2, 'choice2');

%% plot for subjects with low LOOIC (good model prediction)
[fit1, fit2] = compute_choice_prediction_RLbeta4v6(data, parMat, ind_l);
loo_l = mean(loo_RLbeta4_a6(ind_l));
ppc_l = mean(ppc_RLbeta4_a6(ind_l));

f2 = figure; set(f2,'color',[1 1 1]);
subplot(1,2,1); plot_model_prediction(fit1, 'choice1');
subplot(1,2,2); plot_model_prediction(fit2, 'choice2');

%% plot for subjects with medium LOOIC (fare model prediction)
[fit1, fit2] = compute_choice_prediction_RLbeta4v6(data, parMat, ind_m);
loo_m = mean(loo_RLbeta4_a6(ind_m));
ppc_m = mean(ppc_RLbeta4_a6(ind_m));

f3 = figure; set(f3,'color',[1 1 1]);
subplot(1,2,1); plot_model_prediction(fit1, 'choice1');
subplot(1,2,2); plot_model_prediction(fit2, 'choice2');

%% plot for subjects with high LOOIC (poor model prediction)
[fit1, fit2] = compute_choice_prediction_RLbeta4v6(data, parMat, ind_h);
loo_h = mean(loo_RLbeta4_a6(ind_h));
ppc_h = mean(ppc_RLbeta4_a6(ind_h));

f4 = figure; set(f4,'color',[1 1 1]);
subplot(1,2,1); plot_model_prediction(fit1, 'choice1');
subplot(1,2,2); plot_model_prediction(fit2, 'choice2');

%% quantile split, plot
looVec1 = [loo_l, loo_m, loo_h];
ppcVec1 = [ppc_l, ppc_m, ppc_h];

figure; bar(1:3, looVec1); title('LOOIC'); set(gca,'XTickLabel',{'looic_l','looic_m','looic_h'})
figure; bar(1:3, ppcVec1); title('ppc'); set(gca,'XTickLabel',{'looic_l','looic_m','looic_h'})

%% plot trial series for those whose LOOIC are smaller than 30
% ind_ extremely low
ind_el = find(loo_RLbeta4_a6 <= 30);
loo_el = mean(loo_RLbeta4_a6(ind_el));
ppc_el = mean(ppc_RLbeta4_a6(ind_el));

[fit1, fit2] = compute_choice_prediction_RLbeta4v6(data, parMat, ind_el);

f2 = figure; set(f2,'color',[1 1 1]);
subplot(1,2,1); plot_model_prediction(fit1, 'choice1');
subplot(1,2,2); plot_model_prediction(fit2, 'choice2');

plot_trialSeries_RVSL_subj(ind_el)

%% plot trial series for those whose LOOIC are larger than 100
% ind_ extremely high
ind_eh = find(loo_RLbeta4_a6 >= 100);
loo_eh = mean(loo_RLbeta4_a6(ind_eh));
ppc_eh = mean(ppc_RLbeta4_a6(ind_eh));

[fit1, fit2] = compute_choice_prediction_RLbeta4v6(data, parMat, ind_eh);

f2 = figure; set(f2,'color',[1 1 1]);
subplot(1,2,1); plot_model_prediction(fit1, 'choice1');
subplot(1,2,2); plot_model_prediction(fit2, 'choice2');

plot_trialSeries_RVSL_subj(ind_eh)

%% absolute value split, plot
looVec2 = [loo_el, loo_eh];
ppcVec2 = [ppc_el, ppc_eh];

figure; bar(1:2, looVec2); title('LOOIC'); set(gca,'XTickLabel',{'looic_el', 'looic_eh'})
figure; bar(1:2, ppcVec2); title('ppc'); set(gca,'XTickLabel',{'looic_el', 'looic_eh'})



