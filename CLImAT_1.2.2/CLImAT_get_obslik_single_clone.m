function [obslik,condi_probs,condi_probs_fluct] = CLImAT_get_obslik_single_clone(chr,data_baf,data_rd,o,lambda,w,p,depend_table,normal_prior)
%obslik:
%condi_probs
% magic_num = 0.02;

global gender

max_rd = max(data_rd);
N = length(data_rd); %number of data points
w_all = w;%in single clone model, there is only one p

ns = 2; %copy number of stromal cells
tv_S = depend_table(:,2)~=0;
Nc = depend_table(tv_S,3)'; %vector of copy numbers of different entries
mus = 0.5;% rbc mean of stromal cells
Muc = depend_table(tv_S,4)'; %vector of rbc means of different entries
    
%----------------------------------------------------------------------%
%-----------------calculate all of the probablities--------------------%
% fluct_prob = 0.01;
% Het_prior = 0.4;
% Hom_prior = 0.6;
Het_prior = 0.5;
Hom_prior = 0.5;
Num_US = 20; % the number of unique states regulated by global parameters
tv_S = depend_table(:,2)~=0;
cn_u = unique(Nc);

if gender == 1 && (chr == 23 || chr == 24)
    Y = w_all(depend_table(tv_S,2)')+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S);
    Z = w_all(depend_table(tv_S,2)')+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S).*Muc(tv_S);
else
    Y = w_all(depend_table(tv_S,2)')*ns+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S);
    Z = w_all(depend_table(tv_S,2)')*ns*mus+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S).*Muc(tv_S);
end
% Z_homo = w_all(depend_table(tv_S,2)')*ns*mus+(1-w_all(depend_table(tv_S,2)')).*Nc(tv_S);
p_het = Z./Y;
% p_homo = Z_homo./Y;
% p_homo(1) = 0.97;
US_indx = depend_table(tv_S,1);
lambda_c = lambda*Y/2+o;

obslik = zeros(sum(tv_S),N);
condi_probs = zeros(sum(tv_S&(depend_table(:,1)<=Num_US)),N);
condi_probs_fluct = zeros(sum(tv_S&(depend_table(:,1)<=Num_US)),N);

fluct_prob = 0.5;

obslik_bc_Homo = CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),0.997);%0.997
for i=1:length(Y)
    if US_indx(i)<=Num_US %normal state 
        
        if gender == 1 && (chr == 23 || chr == 24) && Muc(i) ~= 1
            obslik(i,:) = (1-fluct_prob)*eps+fluct_prob./((data_baf(2,:)+1)*max_rd);
            condi_probs(i,:) = (1-fluct_prob)*eps./obslik(i,:);
            condi_probs_fluct(i,:) = (fluct_prob./((data_baf(2,:)+1)*max_rd))./obslik(i,:);
        else
           % RD
            j = cn_u == Nc(i);
            obslik_rd = CLImAT_eval_pdf_RD(data_rd,lambda_c(i),p(j));
            %rbc
            temp = find([1 3 8 15] == US_indx(i));
            if ~isempty(temp)
                obslik_bc_Het = normal_prior(temp)*CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),p_het(i)); 
            else
                obslik_bc_Het = CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),p_het(i));
            end  
            obslik_bc = Hom_prior*obslik_bc_Homo+Het_prior*obslik_bc_Het; 
            obslik(i,:) = (1-fluct_prob)*obslik_rd.*obslik_bc+fluct_prob./((data_baf(2,:)+1)*max_rd);
            condi_probs(i,:) = (1-fluct_prob)*obslik_rd.*(Het_prior*obslik_bc_Het)./obslik(i,:);
            condi_probs_fluct(i,:) = (fluct_prob./((data_baf(2,:)+1)*max_rd))./obslik(i,:);
        end
        
    end
end

