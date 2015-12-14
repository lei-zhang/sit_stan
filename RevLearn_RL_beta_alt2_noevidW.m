%% prepare data and initilize parameters

ns = size(data3,3);
nt = size(data3,1);

c2 = squeeze(data3(:,10,:))';

otherC1   = zeros(ns,nt,4);
wOthers   = zeros(ns,nt,4);
cfsC2     = zeros(ns,nt,4);
cfoC2     = zeros(ns,nt,4);

prob_sC2  = zeros(ns,nt,4);
prob_oC2  = zeros(ns,nt,4);
pOth_C2   = zeros(ns,nt,4);

prob_sC2_ori  = zeros(ns,nt,4);
prob_oC2_ori  = zeros(ns,nt,4);
prob_sC2_ori_med  = zeros(ns,nt,4);
prob_oC2_ori_med  = zeros(ns,nt,4);
pOth_C2_ori   = zeros(ns,nt,4);
pOth_C2_ori_med   = zeros(ns,nt,4);

wProb_sC2 = zeros(ns,nt,4);
wProb_oC2 = zeros(ns,nt,4);
with    = zeros(ns,nt);
against = zeros(ns,nt);

for s = 1:ns
    otherC1(s,:,:) = data3(:,6:9,s);
    wOthers(s,:,:) = data3(:,51:54,s);
    cfsC2(s,:,:)   = data3(:,61:64,s);
    cfoC2(s,:,:)   = data3(:,65:68,s);
    prob_sC2_ori(s,:,:) = data3(:,73:76,s);
    prob_oC2_ori(s,:,:) = data3(:,77:80,s);
    prob_sC2_ori_med(s,:,:) = data3(:,101:104,s);
    prob_oC2_ori_med(s,:,:) = data3(:,105:108,s);
end

%% model

for s = 1:ns
    
    for t = 1:nt
        if c2(s,t) == 1
            pOth_C2(s,t,:) = prob_sC2(s,t,:);
            pOth_C2_ori(s,t,:) = prob_sC2_ori(s,t,:);
            pOth_C2_ori_med(s,t,:) = prob_sC2_ori_med(s,t,:);
        elseif c2(s,t) == 2
            pOth_C2(s,t,:) = prob_oC2(s,t,:);
            pOth_C2_ori(s,t,:) = prob_oC2_ori(s,t,:);
            pOth_C2_ori_med(s,t,:) = prob_oC2_ori_med(s,t,:);
        end
    end
end

%% calculate accuracy

pOth_C2_ori = pOth_C2_ori(:,1:99,:);
pOth_C2_ori_med = pOth_C2_ori_med(:,1:99,:);
otherC1 = otherC1(:,2:100,:);
otherC1(otherC1==2) = 0;

cm = zeros(ns,4); % correlation matrix
acc2 = zeros(ns,4);

for s = 1:ns
    m  = squeeze(pOth_C2_ori_med(s,:,:))>0.5; 
    b  = squeeze(pOth_C2_ori(s,:,:))>0.5;
    c  = squeeze(otherC1(s,:,:));
    
    for o= 1:4
        acc2(s,o) = (sum(b(:,o)==c(:,o))) / length(b);
        acc3(s,o) = (sum(m(:,o)==c(:,o))) / length(b);

        
    end
    
    
end

disp(mean(mean(acc2)))