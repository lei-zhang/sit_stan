load('F:\projects\SocialInflu\sit_stan\_outputs\looVec.mat')
looVec(4) = [];

modelName = {'RLbeta\_alt4','RLbeta\_alt3','RLbeta\_alt2','RLbeta\_alt1',...
    'RLcoh','RLnc','RL'};
np = {'8','9','7','7','4','2','2'};

f1 = figure;
set(f1,'color',[1 1 1],'position', [10 10 900 500])

b = barh(7:-1:1,log(looVec),'FaceColor', [0 74 147]/255);
xlabel('LOOIC', 'FontSize', 20)
ylabel('models', 'FontSize', 20)
set(gca,'YTickLabel', modelName, 'XLim', [8000,12000],'XTick', 8000:2000:12000,...
   'YLim', [0 7.5],'YTick', 1:7, 'fontsize',15, 'box','off')

for j = 1:8
    text(-1100,j,np{j},'fontsize',18);
end
text(-1800,9,'#. of params', 'fontsize',15);
