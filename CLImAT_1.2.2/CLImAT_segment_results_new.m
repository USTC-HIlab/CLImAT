function states_seg = CLImAT_segment_results_new(cn_seq, mcn_seq, score_seq, tv_select_seq)
% segmentation based on copy number, B allelic frequency
% 03/23/2014 by Zhenhua

states_seg = [];
if isempty(cn_seq)
    return;
end
pre_cn = cn_seq(1);
pre_mcn = mcn_seq(1);
pre_status = tv_select_seq(1);
start_indx = 1;

for i = 2:length(cn_seq)
    if pre_status == 1 && tv_select_seq(i) == 1
        if ~(pre_cn == cn_seq(i) && pre_mcn == mcn_seq(i))
            indx = start_indx:i-1;
            tv = tv_select_seq(indx) == 1;
            states_seg = [states_seg; start_indx i-1 pre_cn pre_mcn ai_status(pre_cn, pre_mcn) mean(score_seq(indx(tv)))];
            pre_cn = cn_seq(i);
            pre_mcn = mcn_seq(i);
            start_indx = i;
            clear indx tv;
        end
    end
    if pre_status == 1 && tv_select_seq(i) == 0
        if pre_cn ~= cn_seq(i)
            indx = start_indx:i-1;
            tv = tv_select_seq(indx) == 1;
            states_seg = [states_seg; start_indx i-1 pre_cn pre_mcn ai_status(pre_cn, pre_mcn) mean(score_seq(indx(tv)))];
            pre_cn = cn_seq(i);
            pre_mcn = mcn_seq(i);
            pre_status = 0;
            start_indx = i;
            clear indx tv;
        end
    end
    if pre_status == 0 && tv_select_seq(i) == 1
        if pre_cn ~= cn_seq(i)
            states_seg = [states_seg; start_indx i-1 pre_cn pre_mcn ai_status(pre_cn, pre_mcn) mean(score_seq(start_indx:i-1))];
            pre_cn = cn_seq(i);
            start_indx = i;         
        end
        pre_mcn = mcn_seq(i);
        pre_status = 1;
    end
    if pre_status == 0 && tv_select_seq(i) == 0
        if pre_cn ~= cn_seq(i)
            states_seg = [states_seg; start_indx i-1 pre_cn pre_mcn ai_status(pre_cn, pre_mcn) mean(score_seq(start_indx:i-1))];
            pre_cn = cn_seq(i);
            pre_mcn = mcn_seq(i);
            start_indx = i;
        end
    end
end

if start_indx <= length(cn_seq)
    indx = start_indx:length(cn_seq);
    tv = tv_select_seq(indx) == 1;
    if sum(tv) > 0
        score = mean(score_seq(indx(tv)));
    else
        score = mean(score_seq(indx));
    end
    states_seg = [states_seg; start_indx length(cn_seq) pre_cn pre_mcn ai_status(pre_cn, pre_mcn) score];
end

end

function AI = ai_status(cn ,mcn)

if cn < 2
    AI = 1;
elseif cn == mcn
    AI = 2;
else
    AI = 3;
end

end
