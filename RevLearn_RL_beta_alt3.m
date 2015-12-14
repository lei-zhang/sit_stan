%% prepare data and initilize parameters

load('../_data/data3_129.mat')

ns = size(data3,3);
nt = size(data3,1);

c1 = squeeze(data3(:,3,:))';
c2 = squeeze(data3(:,10,:))';
r  = squeeze(data3(:,14,:))';

otherC1   = zeros(ns,nt,4);
otherC2   = zeros(ns,nt,4);
pOth_C2_v1   = zeros(ns,nt,4);
pOth_C2_v2   = zeros(ns,nt,4);

for s = 1:ns
    otherC1(s,:,:) = data3(:,6:9,s);
    otherC2(s,:,:) = data3(:,55:58,s);    
end

%% separate reinforcer model
load('../_outputs/actprob_alt3_v1.mat')
load('../_outputs/actprob_alt3_v2.mat')
for s = 1:ns  
    for t = 1:nt
        if c2(s,t) == 1
            pOth_C2_v1(s,t,:) = prob_sC2(s,t,:);
            pOth_C2_v2(s,t,:) = prob_sC2_v2(s,t,:);
        elseif c2(s,t) == 2
            pOth_C2_v1(s,t,:) = prob_oC2(s,t,:);
            pOth_C2_v2(s,t,:) = prob_oC2_v2(s,t,:);
        end
    end
end


%% calculate accuracy

pOth_C2_v1 = pOth_C2_v1(:,1:99,:);
pOth_C2_v2 = pOth_C2_v2(:,1:99,:);
otherC1 = otherC1(:,2:100,:);
otherC1(otherC1==2) = 0;

cm = zeros(ns,4); % correlation matrix

for s = 1:ns
    
    a  = squeeze(pOth_C2_v1(s,:,:))>0.5;
    b  = squeeze(pOth_C2_v2(s,:,:))>0.5;    
    c  = squeeze(otherC1(s,:,:));
    
    for o= 1:4
        % acc(s,o) = sum(sum([a(:,o) b(:,o)],2)==2) / length(a);
        acc1(s,o) = (sum(a(:,o)==c(:,o))) / length(a);
        acc2(s,o) = (sum(b(:,o)==c(:,o))) / length(b);
    end
    
%     An=bsxfun(@minus,a,mean(a,1));
%     Bn=bsxfun(@minus,b,mean(b,1));
%     An=bsxfun(@times,An,1./sqrt(sum(An.^2,1)));
%     Bn=bsxfun(@times,Bn,1./sqrt(sum(Bn.^2,1)));
%     cm(s,:)=sum(An.*Bn,1);
    
end



%% plot
figure
hist(mean(acc1,2));
hold on
hist(mean(acc2,2));
h = findobj(gca,'Type','patch');
set(h(1),'Facecolor',[1 0 0],'EdgeColor','k');
set(h(2),'Facecolor',[0 0 1],'EdgeColor','k');

return;
