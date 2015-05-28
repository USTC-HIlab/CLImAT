function [LL_all,CLImAT_paras,best_indx] = ...
    CLImAT_main(screening_method,init_CLImAT_paras,depend_table,stepsize_ds,thres_EM,max_iter,verbose,fn_nosuffix)
%--------------------------------------------------------------------%
%------------------>       version 1.2       <---------------------
%--------------------------------------------------------------------%
%2014/10/16 by Zhenhua
%CLImAT main function,basically everything is done here
%--------------------------- screening -------------------------
global clamp_thres
global mc_w
global NoSolutionFlag
clamp_thres = 1-1e-10;
mc_w = 0.8;

if screening_method == 1
elseif screening_method == 2 %d-sampling screening->w0, screen
    %-------------------------------------------------------------------
    %               ---> d-sampling screening <---
    %-------------------------------------------------------------------
    %initialize parameters
    thres_del = 0.07;%0.0006

    init_CLImAT_paras = CLImAT_Init_paras(init_CLImAT_paras,depend_table,0);
    [LL,CLImAT_paras,p_states,num_SNP,aCN] = CLImAT_screening...
        (stepsize_ds,init_CLImAT_paras,depend_table,thres_EM,max_iter,verbose);
    %-------------------------------------------------------------------
    %               ---> summarize ds results <---
    %-------------------------------------------------------------------
    %get detailed info for further ananlysis
    %CN           0 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 0
    N_genotype = [0 1 1 1 2 2 2 2 2 3 3 3 3 3 3 3 4 4 4 4 0];
    tv_del = depend_table(depend_table(:,2)~=0,3)<1;
    US_indx = depend_table(depend_table(:,2)~=0,1);
    
    p_total_del = sum(p_states(tv_del,:),1);
    DF_bias = (1-1./(N_genotype(US_indx)+1))*p_states;

    tv = (p_total_del<thres_del) & (aCN < 4.3);
    NoSolutionFlag=false;
    if ~any(tv)
        warning('Can not find a feasible solution with pre-defined criteria!');
        tv=~tv;
        NoSolutionFlag=true;
    end
    indx = find(tv);
    
    [temp,mindx] = max(LL(tv));
    best_indx = indx(mindx);
    indx = best_indx;
    thres1_ratio = 0.0007; %0.00015
    %correction for possible signal noise
    %     ratio = (log(num_SNP)/2).*(DF_bias(best_indx)-DF_bias)./(LL(best_indx)-LL+eps);
    ratio = zeros(1,length(LL));
    for i=1:length(LL)
        temp1 = LL(best_indx)-LL(i);
        temp2 =(log(num_SNP(i))/2).*(DF_bias(best_indx)-DF_bias(i));
        if abs(temp1)<=1 && abs(temp2)<=0.01
            ratio(i) = 0;
        else
          ratio(i) = temp2/temp1;
        end
    end
    
    candi_best = best_indx;
    ds_candi1 = (ratio>=thres1_ratio) & tv;
    if any(ds_candi1)
        [temp,mindx] = max(ratio(ds_candi1));
        indx = find(ds_candi1);
        best_indx = indx(mindx);
    end
    
%     % for test, save intermediate results
%     %-----------------------------------------------------------------------------------
%     fid = fopen('Model.details','a+');
%     fprintf(fid,'----------------------------------------------------------\n');
%     for i = 1:length(ratio)
%         if i == best_indx
%             fprintf(fid,'2\t%f\t%f\t%f\t%f\n',p_total_del(i),aCN(i),ratio(i),LL(i));
%         elseif i == candi_best
%             fprintf(fid,'1\t%f\t%f\t%f\t%f\n',p_total_del(i),aCN(i),ratio(i),LL(i));
%         else
%             fprintf(fid,'0\t%f\t%f\t%f\t%f\n',p_total_del(i),aCN(i),ratio(i),LL(i));
%         end
%     end
% %     fprintf(fid,'\n');
%     fprintf(fid,'----------------------------------------------------------\n');
%     fclose(fid);
%     %-----------------------------------------------------------------------------------
    
    %-------------------------------------------------------------------
    init_CLImAT_paras = CLImAT_Init_paras(CLImAT_paras,depend_table,best_indx);
    [LL_all,CLImAT_paras] = CLImAT_screening(1,init_CLImAT_paras,depend_table,5*thres_EM,5,verbose);
    best_indx = 1;
    %----------------------------------------------------------------------
elseif screening_method == 3 
else
    error('screening method is not recognized!');
end


function CLImAT_paras = CLImAT_Init_paras(init_CLImAT_paras,depend_table,best_indx)
%this function is used to initialize/process parameters for CLImAT training,
%best_indx is used to indicate which parameter configuration (ususally
%multiple generated in previous screening procedure) are selected. If
%best_indx equals 0, parameters will be initialized
global rd_median
global clamp_thres
global tumor_range
% % S = sum(depend_table(:,2)~=0); %number of states in the table

CLImAT_paras = cell(1,8);
%parameter initialization
if best_indx == 0
    %---w---, ###variable in grid searching###
    if isempty(init_CLImAT_paras{5}) %B w
        w_temp = [0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
        lambda_temp = [rd_median rd_median*2/3 rd_median/2];
        o_temp = 0;


        tv = w_temp >= (1-tumor_range(2)) & w_temp <= (1-tumor_range(1));
        if sum(tv) == 0
            %if no solution candidate, find the closest one, maybe two
            tv = (abs(w_temp-(1-tumor_range(2))) == min(abs(w_temp-(1-tumor_range(2)))));
            tv = tv | (abs(w_temp-(1-tumor_range(1))) == min(abs(w_temp-(1-tumor_range(1)))));
        end
        w_temp = w_temp(tv);

        w_all = repmat(w_temp,1,length(o_temp)*length(lambda_temp));   
        lambda_all = repmat(lambda_temp,length(w_temp)*length(o_temp),1);
        lambda_all = lambda_all(:)';
        tv = lambda_all == rd_median*2/3;
        lambda_all(tv) = rd_median*2./(3-w_all(tv));         
        tv = lambda_all == rd_median/2;
        lambda_all(tv) = rd_median./(2-w_all(tv));         
        o_all = repmat(o_temp,length(w_temp),length(lambda_temp));
        o_all = o_all(:)';

        w_0 = w_all;
    else %
        w_0 = init_CLImAT_paras{5};
    end
    N = size(w_0,2);
    CLImAT_paras{5} = mat2cell(w_0,size(w_0,1),ones(1,size(w_0,2)));

    %---o---
    if isempty(init_CLImAT_paras{3})
        o_0 = o_all;
    elseif length(init_CLImAT_paras{3}) == 1
        o_0 = repmat(init_CLImAT_paras{3},1,N);
    else
        o_0 = init_CLImAT_paras{3};
    end
    CLImAT_paras{3} = mat2cell(o_0,size(o_0,1),ones(1,size(o_0,2)));

    %---lambda---
    if isempty(init_CLImAT_paras{4})
        lambda_0 = lambda_all;
    else
        lambda_0 = init_CLImAT_paras{4};
    end
    CLImAT_paras{4} = mat2cell(lambda_0,size(lambda_0,1),ones(1,size(lambda_0,2)));

    %---pi---
    if isempty(init_CLImAT_paras{1})
        S = sum(depend_table(:,2)~=0);
        prior_0 = 1/(S)*ones(S,1);
        CLImAT_paras{1} = repmat({prior_0},1,N);
    else
        prior_0 = init_CLImAT_paras{1};
        CLImAT_paras{1} = repmat({prior_0},1,N);
    end

    %---A---
    if isempty(init_CLImAT_paras{2})
        S = sum(depend_table(:,2)~=0);
        transmat_0 = norm_trans(ones(S,S),clamp_thres);
%             transmat_0 = ones(S,S)/S;
        CLImAT_paras{2} = repmat({transmat_0},1,N);
    else
        transmat_0 = init_CLImAT_paras{2};
        CLImAT_paras{2} = repmat({transmat_0},1,N);
    end
    
    %---p--- 
    if isempty(init_CLImAT_paras{6})
        tv = depend_table(:,2) ~= 0;
        cn_u = unique(depend_table(tv,3));
        p_0 = ones(length(cn_u),1)*0.5;
    else
        p_0 = init_CLImAT_paras{6};
    end
    CLImAT_paras{6} = repmat({p_0},1,N);
    
    %---indicator vector---
    if isempty(init_CLImAT_paras{7}) %indicator vector: '1' for update '0' fixed
        adj_all = ones(1,6);
    else
        adj_all = init_CLImAT_paras{7};
    end
    CLImAT_paras{7} = repmat({adj_all},1,N);

    if isempty(init_CLImAT_paras{8}) %parameters for observation function and EM algorithm
%         other_paras = [{[1 1.55 1.55 1.55]} 0 0];
        other_paras = [{[1 1.7 1.7 1.7]} 0 0]; %normal prior(0,2,4,6 copies)
    else
        other_paras = init_CLImAT_paras{8};
    end
    CLImAT_paras{8} = repmat({other_paras},1,N);

else %parse the results from previous screening
    for i = 1:length(best_indx)
        %--pi--
        CLImAT_paras{1} = [CLImAT_paras{1} init_CLImAT_paras{1}(best_indx(i))];
        % % %         S = length(init_CLImAT_paras{1}{best_indx(i)});
        % % %         CLImAT_paras{1} = [CLImAT_paras{1} {1/(S)*ones(S,1)}];
        %--A--
        CLImAT_paras{2} = [CLImAT_paras{2} init_CLImAT_paras{2}(best_indx(i))];
        % % %         S = length(init_CLImAT_paras{2}{best_indx(i)});
        % % %         CLImAT_paras{2} = [CLImAT_paras{2} {norm_trans(ones(S,S),clamp_thres)}];
        %--o--
        CLImAT_paras{3} = [CLImAT_paras{3} init_CLImAT_paras{3}(best_indx(i))];
        %--lambda--
        CLImAT_paras{4} = [CLImAT_paras{4} init_CLImAT_paras{4}(best_indx(i))];
        % % %         CLImAT_paras{4} = [CLImAT_paras{4} {0}];
        %--w--
        CLImAT_paras{5} = [CLImAT_paras{5} init_CLImAT_paras{5}(best_indx(i))];
        % % %         CLImAT_paras{7} = [CLImAT_paras{7} {[(0.03^2);(0.03^2)]}];
        %--p--
        CLImAT_paras{6} = [CLImAT_paras{6} init_CLImAT_paras{6}(best_indx(i))];
        %--indicator vector--
        CLImAT_paras{7} = [CLImAT_paras{7} init_CLImAT_paras{7}(best_indx(i))];
        %--parameters for observation function and EM algorithm--
%         CLImAT_paras{8} = [CLImAT_paras{8} {[{[1 1.55 1.55 1.55]} 0 0]}];
        CLImAT_paras{8} = [CLImAT_paras{8} {[{[1 1.7 1.7 1.7]} 0 0]}];
    end
end
