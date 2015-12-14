function plot_trialSeries_RVSL_subj_RLbeta4(subj)
% This founction plots 1st choice, switch not not, outcomes and missing data point

load '\projects\SocialInflu\sit_stan\_data\data3_129.mat'
load '\projects\SocialInflu\sit_stan\_outputs\param_beta4_alt6'
load '\projects\SocialInflu\sit_stan\_outputs\ppc_RLbeta4_a6_pointwise.mat'

delete('trialSeries.ps');

for j = 1:length(subj)
    
    
    % --- load data ------------------------------------------
    data = squeeze(data3(:,:, subj(j) ));
    nt   = size(data,1);
    choice1  = data(:,3);
    missInd  = data(:,41);
    swch     = data(:,5);
    otcm2    = data(:,14);
    reversal = data(:,2);
    
    % --- get internal model variables -----------------------
    [~,~,~,model] = RevLearn_RLbeta4_alt6(parMat(subj(j),:), data);
    myValue = model.myValue;
    otherValue = model.otherValue;
    with    = model.with;
    agst    = model.agst;
    wgtWith = model.wgtWith;
    wgtAgst = model.wgtAgst;
    myValue_beta = model.myValue_beta;
    otherValue_beta = model.otherValue_beta;
    valfun1 = model.valfun1;
    valfun2 = model.valfun2;
    
    %%% --- plot -----
    f(j) = figure;
    set(f(j),'color',[1 1 1], 'position', [20 20 1450 1100]);
    
    %%
    % --- plot choice1's history
    plot(1:nt, data(:,3), 'k:', 'linewidth', 1)
    hold on
    
    nt_1 = find(data(:,3) == 1);
    nt_2 = find(data(:,3) == 2);
    plot(nt_1, ones(length(nt_1),1), 'ko', 'MarkerSize',5, 'MarkerFaceColor', 'm')
    plot(nt_2, 2*ones(length(nt_2),1), 'ko', 'MarkerSize',5, 'MarkerFaceColor', 'g')
    
    
    % --- plot missing choice
    missChoice = choice1(logical(missInd));
    missTrial  = find(missInd==1);
    plot(missTrial,missChoice,'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
    
    % --- plot 2nd outcome
    plot(find(otcm2==1),  otcm2(otcm2==1)*2.5, 'g.')
    plot(find(otcm2==-1), otcm2(otcm2==-1)*(-2.5), 'r.')
    
    % --- plot switch or not
    plot(find(swch==1), swch(swch==1)*3, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b')
    
    % --- plot my Value of option1 and option2
    myV_onset = -3;
    plot(1:nt, myValue(:,1) + myV_onset,'m.:', 'linewidth', 1.5)
    plot(1:nt, myValue(:,2) + myV_onset,'g.:', 'linewidth', 1.5)
    H1 = text(-7,myV_onset, 'my value');
    % set(H1, 'rotation',90, 'position', [-3 -1 0])
    
    % --- plot other Value of option1 and option2
    othV_onset = -7;
    plot(1:nt, otherValue(:,1) + othV_onset, 'm.:', 'linewidth', 1.5)
    plot(1:nt, otherValue(:,2) + othV_onset, 'g.:', 'linewidth', 1.5)
    H2 = text(-8.2, othV_onset, 'other value');
    % set(H2, 'rotation',90, 'position', [-3 -5 0])
    
    % text(101, -1.8, 'option 1: magenta')
    % text(101, -2.3, 'option 2: green')
    
    % --- plot beta-weighted value of option1 and option2
    plot(1:nt, myValue_beta(:,1) + myV_onset,'mo-', 'linewidth', 2, 'MarkerSize',4, 'MarkerEdgeColor',[.5 .5 .5], 'MarkerFaceColor', 'm')
    plot(1:nt, myValue_beta(:,2) + myV_onset,'go-', 'linewidth', 2, 'MarkerSize',4, 'MarkerEdgeColor',[.5 .5 .5], 'MarkerFaceColor', 'g')
    
    plot(1:nt, otherValue_beta(:,1) + othV_onset,'mo-', 'linewidth', 1.5, 'MarkerSize',4, 'MarkerEdgeColor',[.5 .5 .5], 'MarkerFaceColor', 'm')
    plot(1:nt, otherValue_beta(:,2) + othV_onset,'go-', 'linewidth', 1.5, 'MarkerSize',4, 'MarkerEdgeColor',[.5 .5 .5], 'MarkerFaceColor', 'g')
    
    % --- plot valfun1 (valfun1 = beta1 * myValue + beta2 * otherValue )
    valfun1_onset = -0.5; % -8.5;
    plot(1:nt, valfun1(:,1) + valfun1_onset,'mo-', 'linewidth', 2, 'MarkerSize',5, 'MarkerFaceColor', 'm')
    plot(1:nt, valfun1(:,2) + valfun1_onset,'go-', 'linewidth', 2, 'MarkerSize',5, 'MarkerFaceColor', 'g')
    H3 = text(-16,valfun1_onset, 'beta_1*myV + beta_2*othV');
    %set(H3, 'rotation',90, 'position', [-3 -9.5 0])
    
    % --- plot valfun2 (valfun2 = beta3 + beta4*valdiff + beta5*wgtWith + beta6*wgtAgainst)
    valfun2_onset = 11;
    plot(1:nt, valfun2 + valfun2_onset,'bo-', 'linewidth', 2, 'MarkerSize',5, 'MarkerFaceColor', 'b')
    line(xlim, [valfun2_onset, valfun2_onset], 'color','k','lineStyle',':')
    H4 = text(-6,valfun2_onset, 'valfun_2');
    %set(H4, 'rotation',90, 'position', [-3 -11 0])
    
    % --- plot with /against
    w_a_onset = -12.5;
    plot(1:nt, with + w_a_onset, 'g.-')
    plot(1:nt, agst + w_a_onset, 'r.-')
    H5 = text(-16.5, w_a_onset + 0.5, 'unweighted with /against');
    
    % --- plot weighted with / agsint
    w_w_a_onset = -14.5;
    plot(1:nt, wgtWith + w_w_a_onset, 'g.-')
    plot(1:nt, wgtAgst + w_w_a_onset, 'r.-')
    H6 = text(-15.5, w_w_a_onset + 0.5, 'weighted with / against');
    
    line(xlim,[w_a_onset + 1.5, w_a_onset + 1.5], 'color','k','lineStyle','-')
    % line(xlim,[-7 -7], 'color','k','lineStyle','-')
    text(101, w_a_onset, 'with: green')
    text(101, w_a_onset -1, 'against: red')
    
    %% --- plot acc (predictive accuracy)           
    acc = ppc_RLbeta4_a6_pointwise(subj(j),:);
    y_pos = 4;
    plot(1:nt, acc+y_pos, 'b-.', 'linewidth', 2.5)
    
    % line(xlim,[y_pos+1.5, y_pos+1.5], 'color','k','lineStyle','-')
    line(xlim,[y_pos+0.5, y_pos+0.5 ], 'color','k','lineStyle',':')
    
    H7 = text(-14,y_pos+0.5, 'predictive accuracy');
    % set(H7, 'rotation',90, 'position', [-5 y_pos-1.5 0])
    set(gca, 'YTick', [w_w_a_onset, w_w_a_onset+1, w_a_onset, w_a_onset+1,  ...
        othV_onset, myV_onset, valfun1_onset, 1 2 2.5 3, y_pos y_pos+1, valfun2_onset], 'YTickLabel', ...
        {'0','1','0','1','0','0', '0','1st choice: option1','1st choice: option2','outcome','switch', '0.00','1.00', '0'} )
    
    ylim([-16, valfun2_onset + 3])
    
    hold off
    
    %%
    
    
    % --- plot rewardProb reversal
    rev_ind = find(reversal == 1);
    for r = 1:length(rev_ind)
        line([rev_ind(r), rev_ind(r)], ylim, 'color','r','lineStyle','--')
    end
    
    title(sprintf('subject No. %d; PARAM: lr = %3.2f, disc = %3.2f, beta1-6 = [%3.2f %3.2f %3.2f %3.2f %3.2f %3.2f]', ...
        subj(j), parMat(subj(j),:)))
    xlabel('trials')
    
    hold off
    
    % --- save plots into file
    print('-f', '-dpsc2','-append','-loose','-r150', 'trialSeries.ps')
    
end
