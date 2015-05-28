function score = CLImAT_get_reliabilityScore(chr,data_baf,data_rd,het_tv,o,lambda,w,p,cn,mcn)
% 2013/12/28 by Zhenhua Yu 

global gender

ns = 2; %copy number of stromal cells
mus = 0.5;% rbc mean of stromal cells

if gender == 1 && (chr == 23 || chr == 24)
    Y = w+(1-w)*cn;
    Z = w+(1-w)*mcn;
else
    Y = w*ns+(1-w)*cn;
    Z = w*ns*mus+(1-w)*mcn;
end

p_het = Z/Y;
lambda_c = lambda*Y/2+o;
if cn == 0
    j = 1;
elseif cn > 7
    j = 8;
else
    j = cn+1;
end
temp = lambda_c-p(j)/(1-p(j));
if temp < 0
    temp = 0;
end

if sum(het_tv) > 0
    score_bc = CLImAT_eval_pdf_BC(data_baf(1,het_tv),data_baf(2,het_tv),p_het)./CLImAT_eval_pdf_BC(round(data_baf(2,het_tv)*p_het),data_baf(2,het_tv),p_het);
    score_bc(score_bc > 1) = 1;
    score_rd = CLImAT_eval_pdf_RD(data_rd(het_tv),lambda_c,p(j))/CLImAT_eval_pdf_RD(temp,lambda_c,p(j));
    score_rd(score_rd > 1) = 1; 
    score = mean(score_bc.*score_rd);
else
    score_bc = CLImAT_eval_pdf_BC(data_baf(1,~het_tv),data_baf(2,~het_tv),0.997)./CLImAT_eval_pdf_BC(round(data_baf(2,~het_tv)*0.997),data_baf(2,~het_tv),0.997);
    score_bc(score_bc > 1) = 1;
    score_rd = CLImAT_eval_pdf_RD(data_rd(~het_tv),lambda_c,p(j))/CLImAT_eval_pdf_RD(temp,lambda_c,p(j));
    score_rd(score_rd > 1) = 1;
    score = mean(score_bc.*score_rd);
end

if cn > 0 && mcn/cn == 0.5
    score = score*1.2;
end
