function [data_chr_all, data_pos_all, data_bc_all, data_tc_all, data_rd_all] = CLImAT_load_preprocessData(df_file)
% 2014/12/09 by Zhenhua

global window
global gender

global tc_th
global gc_th
global map_th

rd_th = [0.1 99.5];
hom_th = 0.01;

%this format works for Turin data set
Chr_indx = 1;
pos_indx = 2;
bc_indx = 3;
tc_indx = 4;
rd_indx = 5;
gc_indx = 6;
map_indx = 7;

fid = fopen(df_file, 'r');
if fid == -1
    error(['Can not open file ' df_file]);
end
fline = fgetl(fid);
result = regexp(fline,'#window\s*size\s*=\s*(\S+)','tokens','once');
window = str2double(result{1});
results = textscan(fid, '%f%f%f%f%f%f%f', 'HeaderLines', 1, 'treatAsEmpty', {'NA', 'na'});
data_chr_all = results{Chr_indx};% chromosome
data_pos_all = results{pos_indx};% position
data_bc_all = results{bc_indx};% B alellic count
data_tc_all = results{tc_indx};% alellic total count
data_rd_all = results{rd_indx};% read count
data_gc_all = results{gc_indx};% GC-content
data_map_all = results{map_indx};% mappability
clear results;
fclose(fid);

tv = data_chr_all == 24;
if sum(tv) > 0
    gender = 1;
else
    gender = 0;
end
clear tv;

fprintf(1,'Total %d windows are loaded from file "%s".\n',length(data_chr_all),df_file);
if length(data_chr_all) < 200000
    fprintf(1,'Warning: the number of windows loaded is too small, check the fromat of data file and whether data are completely loaded!\n');
end

% initial filtering
Chromosomes = intersect(unique(data_chr_all),1:24); %
low_p = prctile(data_rd_all,rd_th(1));
high_p = prctile(data_rd_all,rd_th(2));
tv = ismember(data_chr_all,Chromosomes) & (data_gc_all > gc_th & data_map_all > map_th(1) & data_map_all < map_th(2)) & (data_rd_all > low_p) & (data_rd_all < high_p) & (data_tc_all >= tc_th(1) & data_tc_all <= tc_th(2));
data_chr_all = data_chr_all(tv);
data_pos_all = data_pos_all(tv);
data_bc_all = data_bc_all(tv);
data_tc_all = data_tc_all(tv);
data_rd_all = data_rd_all(tv);
data_gc_all = data_gc_all(tv);
data_map_all = data_map_all(tv);
clear tv;

% sort by position
for i = 1:length(Chromosomes)
    tv = ismember(data_chr_all, Chromosomes(i));
    temp = sortrows([data_pos_all(tv) data_bc_all(tv) data_tc_all(tv) data_rd_all(tv) data_gc_all(tv) data_map_all(tv)], 1);
    data_pos_all(tv) = temp(:,1);
    data_bc_all(tv) = temp(:,2);
    data_tc_all(tv) = temp(:,3);
    data_rd_all(tv) = temp(:,4);
    data_gc_all(tv) = temp(:,5);
    data_map_all(tv) = temp(:,6);
end
clear tv temp;

% separate corrections of mappability and GC-content, time efficient
% Map-Correction
data_rd_all = MAP_Correction(data_rd_all, data_map_all);
% GC-Correction
data_rd_all = GC_Correction(data_rd_all, data_gc_all);
clear data_map_all data_gc_all;

high_p = prctile(data_rd_all,rd_th(2));
tv = data_rd_all < high_p;
data_chr_all = data_chr_all(tv);
data_pos_all = data_pos_all(tv);
data_bc_all = data_bc_all(tv);
data_tc_all = data_tc_all(tv);
data_rd_all = data_rd_all(tv);

% Quantile normalization of BAF
[data_bc_all, data_tc_all] = CLImAT_BAF_tQN(data_bc_all, data_tc_all);

% make sure that the number of heterozygous positions is comparative with the number of homozygous positions
tv = binopdf(data_tc_all-data_bc_all,data_tc_all,0.01) < hom_th & binopdf(data_bc_all,data_tc_all,0.01) < hom_th;
if sum(tv) < sum(~tv)
    indx = find(~tv); 
    step = max(1,floor(length(indx)/(0.5*sum(tv))));
    ds = 1:step:length(indx);
    tv(indx(ds)) = 1;
    clear indx ds
end

data_chr_all = data_chr_all(tv);
data_pos_all = data_pos_all(tv);
data_bc_all = data_bc_all(tv);
data_tc_all = data_tc_all(tv);
data_rd_all = data_rd_all(tv);
clear tv;

end

