% compare RLbeta_alt4_v6 fitting result: hierarchical vs. individual

load '_outputs\param_beta4_alt6'
load '_outputs\param_beta4_alt6_nh.mat'
load '_outputs\acc_sw.mat'

%% scatter plot of the model fitting parameters
f1 = figure;
set(f1,'color',[1 1 1]);
set(f1,'position',[50 50 1400 800]);

np = size(parMat,2);

par_name = {'lr','disc','beta1','beta2','beta3','beta4','beta5','beta6'};

for j = 1:np
   subplot(2,4,j)
   scatter(parMat(:,j), parMat_nh(:,j),6)
   ylabel('individual fitting')
   xlabel('hierarchical fitting')
   title(par_name{j})
   if j ==1
       xlim([0 .8])
       ylim([0 .8])
   elseif j==2
       xlim([0 .8])
       ylim([0 .8])
   elseif j == 3
       xlim([-10 10])
       ylim([-10 10])
   elseif j == 4
       xlim([-1 3])
       ylim([-1 3])
   elseif j == 5
       xlim([-5.5 -1])
       ylim([-5.5 -1])
   elseif j==6
       xlim([-7 -0])
       ylim([-7 -0])
   elseif j == 7
       xlim([-8 2])
       ylim([-8 2])
   elseif j == 8
       xlim([-6 6])
       ylim([-6 6])
   end
   axis square
end

%% scatter plot of the predictive accuracy only on switch trials
f1 = figure;
set(f1,'color',[1 1 1]);
set(f1,'position',[50 50 400 400]);
scatter(acc_sw, acc_sw_nh,5)
ylabel('individual fitting')
xlabel('hierarchical fitting')
title('predictive accuracy on switch trials')
xlim(ylim)
line([0 .8],[0 .8])
axis square

%% figure out the 'big improvement'
thres = -0.10;
ind =  find((acc_sw - acc_sw_nh) < thres);

np = size(parMat,2);
par_name = {'lr','disc','beta1','beta2','beta3','beta4','beta5','beta6'};

f1 = figure;
set(f1,'color',[1 1 1]);
set(f1,'position',[50 50 1400 800]);

for j = 1:np
   subplot(2,4,j)
   scatter(parMat(:,j), parMat_nh(:,j),6)
   hold on
   scatter(parMat(ind,j), parMat_nh(ind,j),6, 'r', 'filled')
   hold off
   ylabel('individual fitting')
   xlabel('hierarchical fitting')
   title(par_name{j})
   if j ==1
       xlim([0 .8])
       ylim([0 .8])
   elseif j==2
       xlim([0 .8])
       ylim([0 .8])
   elseif j == 3
       xlim([-10 10])
       ylim([-10 10])
   elseif j == 4
       xlim([-1 3])
       ylim([-1 3])
   elseif j == 5
       xlim([-5.5 -1])
       ylim([-5.5 -1])
   elseif j==6
       xlim([-7 -0])
       ylim([-7 -0])
   elseif j == 7
       xlim([-8 2])
       ylim([-8 2])
   elseif j == 8
       xlim([-6 6])
       ylim([-6 6])
   end
   axis square
end












