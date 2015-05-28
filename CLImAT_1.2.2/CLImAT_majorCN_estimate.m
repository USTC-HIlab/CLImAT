function mcn = CLImAT_majorCN_estimate(data_baf, cn, w)

Hom_prior = 0.5;
Het_prior = 1-Hom_prior;
Y = w*2+(1-w)*cn;
i = ceil(cn/2);
LL = zeros(1, cn-i+1);
obslik_bc_Homo = CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),0.997);% homo
for mcn = i:cn
    Z = w+(1-w)*mcn;
    obslik_bc_Het = CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),Z/Y); % het
    obslik_bc = Hom_prior*obslik_bc_Homo+Het_prior*obslik_bc_Het;
    LL(mcn-i+1) = sum(log(obslik_bc));
end

[temp, indx] = max(LL);
mcn = indx+i-1;


end