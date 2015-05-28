function [data_bc, data_tc] = CLImAT_BAF_tQN(data_bc, data_tc)
% Quantile normalization of BAF
% 10/16/2014 by Zhenhua

flag1 = 0;
flag2 = 0;
if size(data_bc,2) > 1
    data_bc = data_bc';
    flag1 = 1;
end
if size(data_tc,2) > 1
    data_tc = data_tc';
    flag2 = 1;
end

data_ac = data_tc-data_bc;
a_fre = data_ac./data_tc;
b_fre = data_bc./data_tc;

% average distance before tQN
tv = b_fre < 0.9 & b_fre > 0.1;
dist1 = sum(abs(sort(a_fre(tv))-sort(b_fre(tv))))/sum(tv);

[temp1,indx1] = sort(a_fre);
[temp2,indx2] = sort(b_fre);

temp = [temp1 temp2];
mean_values = mean(temp,2);
temp1 = zeros(length(a_fre),1);
temp2 = zeros(length(b_fre),1);
temp1(indx1) = mean_values;
temp2(indx2) = mean_values;

% thre = 0.9;
thres = 0.5:0.05:1.5;
dists = zeros(1,length(thres));
for i = 1:length(thres)
    t1 = temp1;
    t2 = temp2;
    tv = t1./(a_fre+eps) > thres(i);
    t1(tv) = thres(i)*a_fre(tv);
    tv = t2./(b_fre+eps) > thres(i);
    t2(tv) = thres(i)*b_fre(tv);

    baf = t2./(t1+t2);
    aaf = 1-baf;

    tv = baf < 0.9 & baf > 0.1;
    x = sort(aaf(tv));
    y = sort(baf(tv));
%     figure(i);
%     clf;
%     plot(x, y, '.');
    dists(i) = sum(abs(x-y))/length(x);
end

[min_dist, I] = min(dists);
best_thre = thres(I);
if dist1 > min_dist
    tv = temp1./(a_fre+eps) > best_thre;
    temp1(tv) = best_thre*a_fre(tv);
    tv = temp2./(b_fre+eps) > best_thre;
    temp2(tv) = best_thre*b_fre(tv);
    a_fre = temp1;
    b_fre = temp2;
    baf = b_fre./(a_fre+b_fre);
    tv = baf > 0 & baf < 1;
    d = round((data_tc(tv).*baf(tv)-data_bc(tv))./(1-baf(tv)));
    data_bc(tv) = data_bc(tv)+d;
    data_tc(tv) = data_tc(tv)+d;
end

if flag1
    data_bc = data_bc';
end
if flag2
    data_tc = data_tc';
end

end