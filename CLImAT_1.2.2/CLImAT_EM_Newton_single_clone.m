function [LL, prior, transmat, o, lambda, w, p, nrIterations] = ...
    CLImAT_EM_Newton_single_clone(stepsize1,init_CLImAT_paras,depend_table, thresh, max_iter,verbose)

%2014/06/10, this new EM algorithm uses univariate method to update 
%parameters sperately, Newton-Raphson method is adopted
global data_rd_sep
global data_baf_sep
global clamp_thres

previous_loglik = -inf;
converged = 0;
num_iter = 1;
LL = [];

if ~iscell(data_baf_sep)
    error('baf data should be stored in cells!');
end

if ~iscell(data_rd_sep)
    error('rd data should be stored in cells!');
end

prior = init_CLImAT_paras{1};
transmat = init_CLImAT_paras{2};
o = init_CLImAT_paras{3};
lambda = init_CLImAT_paras{4};
w = init_CLImAT_paras{5};
p = init_CLImAT_paras{6};
normal_prior = cell2mat(init_CLImAT_paras{8}(1));

while (num_iter <= max_iter) && ~converged
    % perform EM algorithm
    [loglik, exp_num_trans, exp_num_visits1,o_u,lambda_u,w_u,p_u] = ...
        CLImAT_compute_ess_single_clone(stepsize1,prior,transmat,o,lambda,w,p,depend_table,normal_prior);
    
    converged = em_converged_m(loglik, previous_loglik, verbose,thresh);
    
    % update parameters
    if init_CLImAT_paras{7}(1)
        prior = norm_trans(exp_num_visits1',0)';
    end
    if init_CLImAT_paras{7}(2) && ~isempty(exp_num_trans)
        % clamp_thres = 1-1e-4;
        transmat = norm_trans(exp_num_trans,clamp_thres);
    end
    if init_CLImAT_paras{7}(3) %update o here
        o = o_u;
    end
    if init_CLImAT_paras{7}(4) %update lambda here
        lambda = lambda_u;
    end
    if init_CLImAT_paras{7}(5) %update w here
        w = w_u;
    end
    if init_CLImAT_paras{7}(6) %update p here
        p = p_u;
    end
    
%     if verbose
%         disp(['w:' num2str(reshape(w,1,[]))]);
%         disp(['p:' num2str(reshape(p,1,[]))]);
%         disp(['lambda:' num2str(lambda) ', o:' num2str(o)]);
%         fprintf(1, 'iteration %d, loglik = %f\n', num_iter, loglik);
%     end
    
    num_iter =  num_iter + 1;
    previous_loglik = loglik;
    LL = [LL loglik];
end
nrIterations = num_iter - 1;

%--------------------------------------------------------------------------
function [loglik, exp_num_trans, exp_num_visits1,o_u,lambda_u,w_u,p_u] = ...
    CLImAT_compute_ess_single_clone(stepsize1,prior,transmat,o,lambda,w,p,depend_table,normal_prior)
global data_rd_sep
global data_baf_sep
global data_pos_sep
global data_pos_ds_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rd_sep); % each row is a sample
S_all = length(transmat); % number of all states 
exp_num_trans = zeros(S_all,S_all);
exp_num_visits1 = zeros(S_all,1);

%-----------------------E step-----------------------------
gamma_sep = [];
condi_probs_sep = [];
condi_probs_fluct_sep = [];
data_pos_ds_sep = [];
loglik = 0;
N = 0; % the size of the whole data set

for ex=1:numex %
    if stepsize1 >1 %down_screening
        indx_ds = 1:stepsize1:length(data_rd_sep{ex});
        obs_baf = data_baf_sep{ex}(:,indx_ds);
        obs_rd = data_rd_sep{ex}(indx_ds);
        obs_pos = data_pos_sep{ex}(indx_ds);
        N = N+length(indx_ds);
    else %no ds
        obs_baf = data_baf_sep{ex};
        obs_rd = data_rd_sep{ex};
        obs_pos = data_pos_sep{ex};
        N = N+length(obs_rd);
    end

    %condi_probs: Pi(G=k|S=j,O)
    [obslik,condi_probs,condi_probs_fluct] = CLImAT_get_obslik_single_clone(Chromosomes(ex),obs_baf,obs_rd,o,lambda,w,p,depend_table,normal_prior);
    %     [alpha, beta, gamma, current_ll, xi_summed] = fwdback(prior, transmat, obslik);
    %[alpha, gamma, current_ll, beta, xi_summed] = Forward_Backward_Algorithm(prior, transmat, obslik);
    [alpha, gamma, current_ll, beta, xi_summed] = Forward_Backward_Algorithm(prior, transmat, obslik);
    %         [path,current_ll] = viterbi_path(prior, transmat, obslik);
    
    %nonstationary HMM
%     transmat = CLImAT_init_transmat(obs_pos, S_all);
%     [~, gamma, current_ll] = nonStationary_Fwd_Back_Algorithm(prior, transmat, obslik, 0);
    
    clear obslik alpha beta;
    loglik = loglik +  current_ll;
    exp_num_trans = exp_num_trans + xi_summed;
    exp_num_visits1 = exp_num_visits1 + gamma(:,1);
    
    gamma_sep = [gamma_sep {gamma}];
    clear gamma;
    condi_probs_sep = [condi_probs_sep {condi_probs}];
    clear condi_probs;
    condi_probs_fluct_sep = [condi_probs_fluct_sep {condi_probs_fluct}];
    clear condi_probs_fluct;
    data_pos_ds_sep = [data_pos_ds_sep {obs_pos}];
    clear obs_pos;
end

%-----------------------M step-----------------------------
%update global parameters
max_iter = 10;

%update w
w_tol = 0.01;
w_u = CLImAT_update_w_single_clone(stepsize1,o,lambda,w,p,depend_table,w_tol,max_iter);

%update p
p_tol = 0.005;
p_u = CLImAT_update_p_single_clone(stepsize1,o,lambda,w_u,p,depend_table,p_tol,max_iter);

%update lambda
lambda_tol = 0.5;
lambda_u = CLImAT_update_lambda_single_clone(stepsize1,o,lambda,w_u,p_u,depend_table,lambda_tol,max_iter);
lambda_u = round(lambda_u);

% update o
o_tol = 0.5;
o_u = CLImAT_update_o_single_clone(stepsize1,o,lambda_u,w_u,p_u,depend_table,o_tol,max_iter);
o_u = round(o_u);
% o_u = o;

%--------------------------------------------------------------------------
function w_u = CLImAT_update_w_single_clone(stepsize1,o,lambda,w,p,depend_table,w_tol,max_iter)
global data_rd_sep
global data_baf_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global tumor_range
global Chromosomes

numex = length(data_rd_sep); % each row is a sample
ns = 2; %copy number of stromal cells
mus = 0.5;
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3)'; %row vector of copy numbers of different entries
Muc = depend_table(tv,4)';
cn_u = unique(Nc);
p_temp = zeros(1,length(Nc));
for i = 1:length(cn_u)
	tv = cn_u(i) == Nc;
	p_temp(tv) = p(i);
end

iter = 0;
while 1	
	Y = w*ns+(1-w)*Nc;
	Z = w*ns*mus+(1-w)*Nc.*Muc;
	lambda_c = lambda*Y/2+o;
	Y_d = ns-Nc;
	Z_d = ns*mus-Nc.*Muc;
	YZ_d = Y-Z;
	YZ_d(YZ_d == 0) = eps;
	Muc_d = mus-Muc;
	
	part1 = lambda*Y_d.*(1-p_temp)./(2*p_temp);
	part2 = lambda_c.*(1-p_temp)./p_temp;
	part3 = ns*Muc_d.*Nc./Y;
	tv = part2 <= 0;
	part2(tv) = eps;
	tv = part3 <= 0;
	part3(tv) = eps;
   
    %first order differential
    ELL_D_D_1 = 0;
    ELL_B_D_1 = 0;
    %second order differential
    ELL_D_D_2 = 0;
    ELL_B_D_2 = 0;

    for ex=1:numex %
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        if stepsize1 >1 %down_screening
            indx_ds = 1:stepsize1:length(data_rd_sep{ex});
            obs_baf = data_baf_sep{ex}(:,indx_ds);
            obs_rd = data_rd_sep{ex}(indx_ds);
        else %no ds
            obs_baf = data_baf_sep{ex};
            obs_rd = data_rd_sep{ex};
        end
        post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));
        post_probs = gamma_sep{ex}(1:length(Y),:).*condi_probs_sep{ex}(1:length(Y),:); 
%         post_probs = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));
        for i=1:length(Y) % 05/28/2010, now only calcualte heter entries in the dependent table
            %RD
            %first order differential
            ELL_D_D_1 = ELL_D_D_1 + sum(post_probs_not_fluct(i,:)*part1(i).*(psi(obs_rd+part2(i))+log(1-p_temp(i))-psi(part2(i))));
            %second order differential
            ELL_D_D_2 = ELL_D_D_2 + sum(post_probs_not_fluct(i,:)*part1(i)^2.*(psi(1,obs_rd+part2(i))-psi(1,part2(i))));
            %rbc
            %first order fifferential
            ELL_B_D_1 = ELL_B_D_1 + sum(post_probs(i,:).*(part3(i)*(obs_baf(1,:)/Z(i)-(obs_baf(2,:)-obs_baf(1,:))/YZ_d(i))));
            %second order differential
			ELL_B_D_2 = ELL_B_D_2 + sum(post_probs(i,:).*(part3(i)*(-1*obs_baf(1,:)*Z_d(i)/Z(i)^2-(obs_baf(2,:)-obs_baf(1,:))*(Z_d(i)-Y_d(i))/YZ_d(i)^2-...
                        Y_d(i)/Y(i)*(obs_baf(1,:)/Z(i)-(obs_baf(2,:)-obs_baf(1,:))/YZ_d(i)))));
        end
    end
    ELL_ALL_D_1 = ELL_D_D_1+ELL_B_D_1;    
    ELL_ALL_D_2 = ELL_D_D_2+ELL_B_D_2;
    
%     w_u = w-ELL_ALL_D_1/ELL_ALL_D_2;
%     if isnan(w_u)
%         w_u = w;
%     end
%     if w_u < 1-tumor_range(2)
%         w_u = 1-tumor_range(2);
%     end
%     if w_u > 1-tumor_range(1)
%         w_u = 1-tumor_range(1);
%     end
    
    w_adj = -ELL_ALL_D_1/ELL_ALL_D_2;
    
    % now determine if the update violates the constrains
    % let w_u = w+c*w_adj, c is the minimal coefficient from a list of
    % coefficients c_all that activiate but do not violate the constrains
    % 1-tumor_range(2)<=w<=1-tumor_range(1)
    c_p = 1;
    
    %-w<= tumor_range(2)-1 => w+c*w_adj>=1-tumor_range(2)
    if -w_adj > 0 % not a feasible direction
        temp = (1-tumor_range(2)-w)/w_adj;
        if temp < 1
%             disp(['Constrain: w>' num2str(1-tumor_range(2)) ') is active!!']);
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    %w<=1-tumor_range(1) => w+c*w_adj<=1-tumor_range(1)
    if w_adj > 0 % not a feasible direction
        temp = (1-tumor_range(1)-w)/w_adj;
        if temp < 1
%             disp(['Constrain: w<' num2str(1-tumor_range(1)) ' is active!!']);
            if temp < c_p
                c_p = temp;
            end
        end
    end

    w_u = w+c_p*w_adj;
    
    iter = iter+1;
    
    if abs(w_u-w)<w_tol || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        w = w_u;
    end
end

%--------------------------------------------------------------------------
function o_u = CLImAT_update_o_single_clone(stepsize1,o,lambda,w,p,depend_table,o_tol,max_iter)
global data_rd_sep
global gamma_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rd_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3)'; %row vector of copy numbers of different entries
Y = w*ns+(1-w)*Nc;
cn_u = unique(Nc);
p_temp = zeros(1,length(Y));
for i = 1:length(cn_u)
    tv = cn_u(i) == Nc;
    p_temp(tv) = p(i);
end
part1 = (1-p_temp)./p_temp;

iter = 0;
while 1  
    %first order differential
    ELL_D_D_1 = 0;
    %second order differential
    ELL_D_D_2 = 0;
    
    lambda_c = lambda*Y/2+o;
    for ex=1:numex %
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        if stepsize1 >1 %down_screening
            indx_ds = 1:stepsize1:length(data_rd_sep{ex});
            obs_rd = data_rd_sep{ex}(indx_ds);
        else %no ds
            obs_rd = data_rd_sep{ex};
        end
        post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));       
        for i=1:length(Y) % 05/28/2010, now only calcualte heter entries in the dependent table
            %RD
            %first order differential
			ELL_D_D_1 = ELL_D_D_1 + sum(post_probs_not_fluct(i,:)*part1(i).*(psi(obs_rd+lambda_c(i)*part1(i))...
                +log(1-p_temp(i))-psi(lambda_c(i)*part1(i))));
            %second order differential
            ELL_D_D_2 = ELL_D_D_2 + sum(post_probs_not_fluct(i,:)*part1(i)^2.*(psi(1,obs_rd+lambda_c(i)*part1(i))...
                -psi(1,lambda_c(i)*part1(i))));
        end
    end
    
    o_adj = -ELL_D_D_1/ELL_D_D_2;
    c_p = 1;
    %-o<=-1 => o+c*o_adj>=1
    if -o_adj > 0 % not a feasible direction
        temp = (1-o)/o_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    %o<lambda/10 => o+c*o_adj<lambda/10
    if o_adj > 0 % not a feasible direction
        temp = (lambda/10-o)/o_adj;
        if temp < 1
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    o_u = o+c_p*o_adj;
    iter = iter+1;
    
    if abs(o_u-o)<o_tol || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        o = o_u;
    end
end

%--------------------------------------------------------------------------
function lambda_u = CLImAT_update_lambda_single_clone(stepsize1,o,lambda,w,p,depend_table,lambda_tol,max_iter)
global data_rd_sep
global gamma_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rd_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3)'; %row vector of copy numbers of different entries
Y = w*ns+(1-w)*Nc;
cn_u = unique(Nc);
p_temp = zeros(1,length(Y));
for i = 1:length(cn_u)
    tv = cn_u(i) == Nc;
    p_temp(tv) = p(i);
end
part1 = (1-p_temp)./p_temp;

iter = 0;
while 1  
    %first order differential
    ELL_D_D_1 = 0;
    %second order differential
    ELL_D_D_2 = 0;
    
    lambda_c = lambda*Y/2+o;

    for ex=1:numex %
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        if stepsize1 >1 %down_screening
            indx_ds = 1:stepsize1:length(data_rd_sep{ex});
            obs_rd = data_rd_sep{ex}(indx_ds);
        else %no ds
            obs_rd = data_rd_sep{ex};
        end
        post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));       
        for i=1:length(Y) % 05/28/2010, now only calcualte heter entries in the dependent table
            %RD
            %first order differential
			ELL_D_D_1 = ELL_D_D_1 + sum(post_probs_not_fluct(i,:)*Y(i)*part1(i)/2.*(psi(obs_rd+lambda_c(i)*part1(i))...
                +log(1-p_temp(i))-psi(lambda_c(i)*part1(i))));
            %second order differential
            ELL_D_D_2 = ELL_D_D_2 + sum(post_probs_not_fluct(i,:)*(Y(i)*part1(i)/2)^2.*(psi(1,obs_rd+lambda_c(i)*part1(i))...
                -psi(1,lambda_c(i)*part1(i))));
        end
    end
    lambda_adj = -ELL_D_D_1/ELL_D_D_2;
    c_p = 1;
    %-lambda<=-1 => lambda+c*lambda_adj>=1
    if -lambda_adj > 0 % not a feasible direction
        temp = (1-lambda)/lambda_adj;
        if temp < 1
%             disp('Constrain: lambda>0 is active!!');
            if temp < c_p
                c_p = temp;
            end
        end
    end
    
    lambda_u = lambda+c_p*lambda_adj;
    
    iter = iter+1;
	
    if abs(lambda_u-lambda)<lambda_tol || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        lambda = lambda_u;
    end
end

%--------------------------------------------------------------------------
function p_u = CLImAT_update_p_single_clone(stepsize1,o,lambda,w,p,depend_table,p_tol,max_iter)
global data_rd_sep
global gamma_sep
global condi_probs_fluct_sep
global Chromosomes

numex = length(data_rd_sep); % each row is a sample
ns = 2; %copy number of stromal cells
Num_US = 20; % the number of unique states regulated by global parameters
tv = (depend_table(:,1)<=Num_US);
Nc = depend_table(tv,3)'; %row vector of copy numbers of different entries
Y = w*ns+(1-w)*Nc;
cn_u = unique(Nc);
lambda_c = Y*lambda/2+o;

iter = 0;
p_u = zeros(size(p));
unConverged = 1:length(p);
while 1  
    %first order differential
    ELL_D_D_1 = zeros(length(p),1);
    %second order differential
    ELL_D_D_2 = zeros(length(p),1);
	p_temp = zeros(1,length(Y));
	for i = 1:length(cn_u)
		tv = cn_u(i) == Nc;
		p_temp(tv) = p(i);
	end
    part1 = lambda_c.*(1-p_temp)./p_temp;
	tv = part1 <= 0;
	part1(tv) = eps;
    for ex=1:numex %
        if Chromosomes(ex) == 23 || Chromosomes(ex) == 24
            continue;
        end
        if stepsize1 >1 %down_screening
            indx_ds = 1:stepsize1:length(data_rd_sep{ex});
            obs_rd = data_rd_sep{ex}(indx_ds);
        else %no ds
            obs_rd = data_rd_sep{ex};
        end
        post_probs_not_fluct = gamma_sep{ex}(1:length(Y),:).*(1-condi_probs_fluct_sep{ex}(1:length(Y),:));       
        
        for i=unConverged % 05/28/2010, now only calcualte heter entries in the dependent table
			indx = find(cn_u(i) == Nc);
			for j = indx
				%RD
				%first order differential
				ELL_D_D_1(i) = ELL_D_D_1(i) + sum(post_probs_not_fluct(j,:).*(lambda_c(j)/p_temp(j)^2*(psi(part1(j))...
					-psi(obs_rd+part1(j))-log(1-p_temp(j)))+(obs_rd-lambda_c(j))/p_temp(j)));
				%second order differential
				ELL_D_D_2(i) = ELL_D_D_2(i) + sum(post_probs_not_fluct(j,:)/p_temp(j)^2.*(2*lambda_c(j)/p_temp(j)*(psi(obs_rd+part1(j))+log(1-p_temp(j))-psi(part1(j)))...
					+lambda_c(j)*(lambda_c(j)/p_temp(j)^2*psi(1,obs_rd+part1(j))+1/(1-p_temp(j))-lambda_c(j)/p_temp(j)^2*psi(1,part1(j)))-obs_rd+lambda_c(j)));
			end
        end
    end
    p_adj = -ELL_D_D_1./ELL_D_D_2;
    iter = iter+1;
    tv = ones(1,length(unConverged)) > 0;
    j = 0;
    for i=unConverged
        j = j+1;
        c_p = 1;
        %-p(j)<=-0.001 => p(j)+c*p_adj(j)>=0.001
        if -p_adj(i) > 0 % not a feasible direction
            temp = (0.01-p(i))/p_adj(i);
            if temp < 1
%                 disp(['Constrain: p(' num2str(i) ')>0 is active!!']);
                if temp < c_p
                    c_p = temp;
                end
            end
        end
        %p(j)<=0.999 => p(j)+c*p_adj(j)<=0.999
        if p_adj(i) > 0 % not a feasible direction
            temp = (0.99-p(i))/p_adj(i);
            if temp < 1
    %             disp(['Constrain: p(' num2str(i) ')<1 is active!!']);
                if temp < c_p
                    c_p = temp;
                end
            end
        end
        p_u(i) = p(i)+c_p*p_adj(i);
        if isnan(p_u(i))
            p_u(i) = p(i);
        end
             
        if abs(p_u(i)-p(i))<p_tol
            tv(j) = 0;
        end
    end
    unConverged = unConverged(tv);
    
    if isempty(unConverged) || iter>max_iter
%         disp(['The Newton-Raphson method is done in ' num2str(iter) ' iterations'])
        break;
    else
        p = p_u;
    end
end

