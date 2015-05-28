function het_seq = CLImAT_het_estimate(chr, data_baf, cn, mcn, w)

global gender

if gender == 1 && (chr == 23 || chr == 24)
    Y = w+(1-w)*cn;
    Z = w+(1-w)*mcn;
else
    Y = w*2+(1-w)*cn;
    Z = w+(1-w)*mcn;
end
obslik_bc_Homo = CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),0.997);% homo
if cn > 0 && mcn/cn == 0.5
    obslik_bc_Het = 1.55*CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),Z/Y); % het
else
    obslik_bc_Het = CLImAT_eval_pdf_BC(data_baf(1,:),data_baf(2,:),Z/Y); % het
end
het_seq = obslik_bc_Het > obslik_bc_Homo;

end