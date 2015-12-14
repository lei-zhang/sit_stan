# sit_stan
project 'sit': Bayesian hierarchical fitting with Stan

This folder contains all the existing sit(social influence task)-related models. Below is the model index, in a chronological order based on development history.

RevLearn_

 - RL: simple RL model with Rescorla Wagner update rule.
 
 - RLnc: ficticious RL model, also model the non-chosen choice (-reward - value[non-chosen]).
 - RLnc_2lr: add an additional learning rate to the non-chosen update
 - RLnc_cfa: add an additional conterfactual attention parameter, (-cfa*reward)
 - RLnc_2lr_cfa: putting 2lr and cfa together
 
 - RLcoh: reweight the values according to decision type (with vs. against)
 - RLcoh_2lr
 - RLcoh_cfa
 - RLcoh_2lr_cfa
 
 - RLcumrew: reweight the values according to coplayers' reward history (cumulative reward)
 - RLcumrew_2lr
 - RLcumrew_cfa
 - RLcumrew_2lr_cfa
 
 - RLcoh_modvalue: the coherence reweight goes into the value directly, the same as RLcoh
 - RLcoh_modprob_tempin: the coherence reweight only goes into the softmax transformation, temperature*value + reweight
 - RLcoh_modprob_tempout: the coherence reweight only goes into the softmax transformation, temperature*(value + reweight)
 
 - RLbeta_alt1: weighted sum of otherRewards according to subjects' preference towards co-players
 
 - RLbeta_alt2: use beta_cdf to quantify others' choice tendency, and then used to model otherValue
 - RLbeta_alt2_c_v1: with evidence weight
 - RLbeta_alt2_c_v2: without evidence weight
 
 - RLbeta_alt3: treat each co-player as a separate reinforcer, use a simple RL update to quantify otherValue
 - RLbeta_alt3_p1_v1: step1, get lr and tau for the others, single lr and tau
 - RLbeta_alt3_p2_v1: 1 lr for chosen and non-chosen choices, with c_rep 
 - RLbeta_alt3_p1_v2: step1, get lr and tau for the others, 4 lr and tau
 - RLbeta_alt3_p2_v2: 1 lr for chosen and non-chosen choices, without c_rep; essentially, this should be the same as _p2_v1
 
 - RLbeta_alt4: (preference) weighted sum of cumulative otherRewards to represent otherValue
 - RevLearn_RLbeta_alt4_c_w_v1_1lr: only uses current preference weights and current with/against
 - RevLearn_RLbeta_alt4_c_w_v2_1lr: preference weight sum up to one --> [3 2 1 1]/7
 - RevLearn_RLbeta_alt4_c_w_v3_1lr: v3_without weight on otherValue when update
 - RevLearn_RLbeta_alt4_c_w_v4_1lr: use actual weight and with/against
 - RevLearn_RLbeta_alt4_c_w_v5_1lr: beta1-4 separately model myValue and otherValue (chosen vs. non-chosen)
 - RevLearn_RLbeta_alt4_c_w_v6_1lr: actual weight but only current with/against
 - RevLearn_RLbeta_alt4_c_v6_1lr  : use actual with/against, rather than weighted with/against 
 - RevLearn_RLbeta_alt4_c_w_v6_2lr: actual weight but only current with/against, 2 learning rate
 - RevLearn_RLbeta_alt4_c_w_v7_1lr: based on v1, windows size = 4
 - RevLearn_RLbeta_alt4_c_w_v8_1lr: same beta for myV, separate beta for otherV
 - RevLearn_RLbeta_alt4_c_w_v9_1lr: based on v1, windows size = 5
 - RevLearn_RLbeta_alt4_c_w_v10_1lr: based on v1, adding a cfa parameter

 - RevLearn_RLbeta_alt4_c_w_v11_1lr: based on v6, normalize the discounted cumulative reward (cr), before using for update otherValue (divided by sum)
 - RevLearn_RLbeta_alt4_c_w_v12_1lr: based on v6, use [-1 1] as the discounted reward, rather than [0 1], normalized cr with SOFTMAX (sumup to one)
 - RevLearn_RLbeta_alt4_c_w_v13_1lr: based on v6, valdiff = valfun1(c1) - valfun(~c1), instead of valdiff = myValue(c1) - myValue(~c1)
 - RevLearn_RLbeta_alt4_c_w_v14_1lr: based on v6, normalize the discounted cumulative reward, 2 * SOFTMAX - 1, --> [-1 1]
 - RevLearn_RLbeta_alt4_c_w_v15_1lr_ind: based on v6, fit the model non-hierarchically, (individual fitting)
 - RevLearn_RLbeta_alt4_c_w_v16_1lr: based on v6, valdiff = bet * ( V(c1) - V(~c1) )
 - RevLearn_RLbeta_alt4_c_w_v17_1lr: based on v6, use 2 beta for Valdiff (+ and -) --> beta4 * (vdiff >=0) * vdiff + beta5 * (vdiff<0) * vdiff
 - RevLearn_RLbeta_alt4_c_w_v18_1lr: based on v6, use inv_logit for normalizing cr
 - RevLearn_RLbeta_alt4_c_w_v19_1lr: based on v6, use inv_logit for normalizing otherValue, not otherCR
 - RevLearn_RLbeta_alt4_c_w_v20_1lr: based on v6, separate beta for general bias: when nWigh >= nAgst, and when aWith < nAgst
 - RevLearn_RLbeta_alt4_c_w_v21_1lr: based on v12 & v20, separate beta for general bias, using [-1 1] for cr, inv_logit for normalizing
 - RevLearn_RLbeta_alt4_c_w_v22_1lr: based on v6, using other's accuracy as the weight for updating otherValue
 - RevLearn_RLbeta_alt4_c_w_v23_1lr: based on v6, adding two temperatures for categorical_logit() and for bernoulli_logit()

 - RevLearn_RLbeta_alt4_c_w_v24_1lr: based on v19, use (inv_logit *2 -1) for normalizing otherValue, not otherCR
 - RevLearn_RLbeta_alt4_c_w_v25_1lr: --> v24, + 1)different valdiff betas, 2)accuracy as weight, 3)valdiff=diff(valfun), 4) diff(w/a) 2 betas, 5) [-1 1] for cr
 - RevLearn_RLbeta_alt4_c_w_v26_1lr: --> v24, + 1)valdiff=diff(valfun), 2) diff(w/a) 2 betas, 
 - RevLearn_RLbeta_alt4_c_w_v27_1lr: --> v24, + 1)valdiff=diff(valfun), 2) diff(w/a) 2 betas, 3) [-1 1] for cr
 - RevLearn_RLbeta_alt4_c_w_v28_1lr: based on v24, + 1)valdiff=diff(valfun), 2)only uses nAgainst
 - RevLearn_RLbeta_alt4_c_w_v29_1lr: based on v26, without general bias
 - RevLearn_RLbeta_alt4_c_w_v30_1lr: based on v28, without general bias

 - RevLearn_RLbeta_alt4_c_w_v31_1lr: based on v6, add (1) 'inv_logit *2 -1' for normalizing otherValue, (2) valdiff=diff(valfun)
 - RevLearn_RLbeta_alt4_c_w_v32_1lr: based on v31, use something like the RLcoh, 1/2 for second choice, not 0/1 switch 
 - RevLearn_RLbeta_alt4_c_w_v33_1lr: based on v32, without beta[3,s] * valfun2[:]
 - RevLearn_RLbeta_alt4_c_w_v34_1lr: based on v31, leave out the general bias
 - RevLearn_RLbeta_alt4_c_w_v35_1lr: based on v31, valdiff = bet1 * ( valfun1(c1) - valfun1(~c1) )
 - RevLearn_RLbeta_alt4_c_w_v36_1lr: based on v31, otherReward2, use [-1 1] for the current trial, use [0 1] for the past trial
 - RevLearn_RLbeta_alt4_c_w_v37_1lr: based on v36, 4) diff(w/a) 2 betas
 - RevLearn_RLbeta_alt4_c_w_v38_1lr: based on v37, valdiff = bet1 * valdiff

 
 - RLbeta_alt5: (preference) weighted sum of cumulative otherRewards to represent otherValue, bind beta5 and beta6 together
 - RevLearn_RLbeta_alt5_c_w_v1_1lr: no 'general bias', run: vfun2 <- beta3*vDiff + beta4*(wgtWigh - wgtAgst)
 - RevLearn_RLbeta_alt5_c_w_v3_1lr: with 'general bias', run: vfun2 <- beta3 + beta4*vDiff + beta5*(wgtWigh - wgtAgst)


 
 
 
 
 - _w: weighted coherence, i.e. weight .* with, weight .* against
 - _n: normalised other Value, i.e. othV(c2) / (othV(c2), othV(~c2))
 
 
 
