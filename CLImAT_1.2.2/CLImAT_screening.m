function [LL_all,CLImAT_paras,p_states,num_SNP,aCN,p_s] = CLImAT_screening(stepsize1,init_CLImAT_paras, depend_table,thres1,max_iter1,verbose1)

%---------------------run the algorithm------------------------------
%1xN cell vectors
prior_all = init_CLImAT_paras{1};
transmat_all = init_CLImAT_paras{2};
o_all = init_CLImAT_paras{3};
lambda_all = init_CLImAT_paras{4};
w_all = init_CLImAT_paras{5};
p_all = init_CLImAT_paras{6};
indivec_all = init_CLImAT_paras{7};
otherparas_all = init_CLImAT_paras{8};

LL_all = [];
CLImAT_paras = cell(1,8);
t_all = 0; 
if nargout >2
    p_states = [];
    num_SNP = [];
    aCN = [];
    p_s = [];
end

for i=1:length(o_all)
%     if verbose1
%         tic
%     end
    
%     disp(['parameter ' num2str(i)]);
    
    %1x1 cell
    init_CLImAT_paras(1) = prior_all(i);
    init_CLImAT_paras(2) = transmat_all(i);
    init_CLImAT_paras(3) = o_all(i);
    init_CLImAT_paras(4) = lambda_all(i);
    init_CLImAT_paras(5) = w_all(i);
    init_CLImAT_paras(6) = p_all(i);
    init_CLImAT_paras(7) = indivec_all(i);
    init_CLImAT_paras(8) = otherparas_all(i);
    
    [LL, prior, transmat, o, lambda, w, p, nrIterations] = CLImAT_EM_Newton_single_clone(stepsize1,init_CLImAT_paras,depend_table,thres1,max_iter1,verbose1);
        
    LL_all = [LL_all LL(end)];
    CLImAT_paras{1} = [CLImAT_paras{1} {prior}];
    CLImAT_paras{2} = [CLImAT_paras{2} {transmat}];
    CLImAT_paras{3} = [CLImAT_paras{3} {o}];
    CLImAT_paras{4} = [CLImAT_paras{4} {lambda}];
    CLImAT_paras{5} = [CLImAT_paras{5} {w}];
    CLImAT_paras{6} = [CLImAT_paras{6} {p}];
    CLImAT_paras{7} = [CLImAT_paras{7} init_CLImAT_paras(7)];
    CLImAT_paras{8} = [CLImAT_paras{8} init_CLImAT_paras(8)];
    
    if nargout >2
        [temp1,temp2,temp3,temp4] = CLImAT_process_results_new(depend_table);
        p_states = [p_states temp1];
        num_SNP = [num_SNP temp2];
        aCN = [aCN temp3];
        p_s = [p_s temp4];
    end

    if verbose1
%         t = toc;
        fprintf(1, '\n');
        disp('--------------- screening report -----------------')
        disp(['run ' num2str(i) ' done, w:' num2str(reshape(w,1,[])) ', lambda:' num2str(lambda) ', o:' num2str(o)]);
        disp(['p:' num2str(reshape(p,1,[]))]);
        disp(['LL:' num2str(LL(end),'%5.1f')]);
        disp('--------------- screening report -----------------')
%         t_all = t_all+t;
    end
    
end