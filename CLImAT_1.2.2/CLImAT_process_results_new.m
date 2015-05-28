function [p_states,num_SNP,aCN,p_s,states_pre,state_seq,het_seq] = CLImAT_process_results_new(depend_table)
%------over-all information of the cancer sample------
%p_states: proportions of all hidden states
%num_SNP: total number of SNPs investigated
%aCN: averaged copy number
%p_s: proportions of chromosomal regions

global data_pos_ds_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep

state_seq=[];
het_seq=[];

Num_US = 20;
tv_S = depend_table(:,2)~=0;
US_indx = depend_table(tv_S,1);%uniqe state mapped to states in CLImAT
CN_mapping = depend_table(tv_S,3)'; %make sure it's 1xN

%initialize output parameters
p_states = [];
num_SNP = [];
aCN = [];
p_s = [];
if nargout > 4%
    states_pre = [];
end

%initialize intermediate variables
current_num_SNP = 0;
exp_num_states = [];
if nargout > 4
    current_states_pre = [];
end

weight_len = 0;
total_len = 0;
for i = 1:length(gamma_sep) %for the ith chromosome
    post_probs = gamma_sep{i};
    data_pos = data_pos_ds_sep{i};

    %---handle p_states and num_SNP---
    if isempty(exp_num_states) %initialization
        exp_num_states = zeros(size(post_probs,1),1);
    end
    exp_num_states = exp_num_states+sum(post_probs,2);
    current_num_SNP = current_num_SNP+size(post_probs,2);

    [temp,MAP_state] = max(post_probs,[],1);
    state_pre_seg = CLImAT_segment_results(MAP_state,0);
    temp = data_pos(state_pre_seg(:,2))-data_pos(state_pre_seg(:,1))+1;
    weight_len  = weight_len + temp*CN_mapping(state_pre_seg(:,3))';
    total_len = total_len + sum(temp);
    %---handle MAP states---
    if nargout > 4 %output predicted MAP states for both clones
        state_seq = [state_seq {MAP_state}];
        condi_probs = zeros(size(MAP_state));
        condi_probs_fluct = zeros(size(MAP_state));
        for k = 1:length(MAP_state)
            condi_probs(k) = condi_probs_sep{i}(MAP_state(k),k);
            condi_probs_fluct(k) = condi_probs_fluct_sep{i}(MAP_state(k),k);
        end
        het_seq = [het_seq {condi_probs > (1-condi_probs-condi_probs_fluct)}];
        for k = 1:size(state_pre_seg,1)
            current_states_pre = [current_states_pre;[i,state_pre_seg(k,1:2),...
                US_indx(state_pre_seg(k,3)),0,1]];
        end
    end %if nargout > 4%output predicted MAP states

end % for i=1:length(gamma) %for the ith chromosome

%---handle p_states and num_SNP---
current_p_states = exp_num_states/current_num_SNP;
p_states = [p_states current_p_states];
num_SNP = [num_SNP current_num_SNP];
ploidy = weight_len/total_len;

%---handle aCN---
aCN = [aCN ploidy];

%---handle p_s---
normal_tv = (depend_table(tv_S,1) == 3) & (depend_table(tv_S,2 )== 1); %normal state
p_s = [p_s [sum(current_p_states((depend_table(tv_S,1) <= Num_US) & (~normal_tv)));0;0;...
        sum(current_p_states(normal_tv));sum(current_p_states(depend_table(tv_S,1) > Num_US))]];

%---handle MAP states---
if nargout > 4%output predicted MAP states
    states_pre = [states_pre {current_states_pre}];
end
