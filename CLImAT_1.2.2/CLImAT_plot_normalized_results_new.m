function CLImAT_plot_normalized_results_new(Datafile,resultsfile,plotsdir,barcode)
%this function is used to plot 22 figures given the data as well
%as the results from NormCNV, all these figures are stored in a specified
%directory
%%%
%----------------------read results  ------------------------%
fid = fopen(resultsfile,'r');
if fid == -1
    error(['Can not open result file: ' resultsfile]);
end

%get estimated global parameters from the first row of the result file
lambda = [];
o = [];
w = [];

while 1
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'StartPos')),break,end
    %Lambda
    result1 = regexp(tline,'Copy neutral read depth:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        lambda = str2double(result1{1});
    end
    %RD shift
    result1 = regexp(tline,'Shift of read depth:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        o = str2double(result1{1});
    end
    %w1 and w2
    result1 = regexp(tline,'Proportion of abnormal cells in the sample:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        w = str2double(result1{1});
    end
end
%report errors if these values are not parsed successfully
if isempty(lambda)
    error(['Can not read estimated Read depth for normal copy from ',resultsfile]);
end
if isempty(o)
    error(['Can not read estimated RD shift from ',resultsfile]);
end
if length(w)<1
    error(['Can not read estimated tumor proportion from ',resultsfile]);
end

%then read the results
results = textscan(fid,'%f %f %f %f %f %*f %f %*f %*s %*s %*f %*f','treatAsEmpty', {'NA', 'na'});
fclose(fid);
chr_seg = results{1};
pstart_seg = results{2};
pend_seg = results{3};
cn_seg = results{4};
baf_seg = results{5};
% AI_seg = results{6};
score_seg = results{6};
clear results;

%load data
eval(['load ' Datafile]);

%--------------- plot figures ---------------------%

% HOMD, HEMD, NEUT, AMP
rd_colors = [%15 111 50;
             153 192 13;
             33 172 71;
             115 30 225;
% 			 81 71 154;
             231 31 25]./255;
% HOMD, HEMD, HET, NLOH, ALOH
baf_colors = [%15 111 50;
              153 192 13;
              33 172 71;
              38 68 180;
              20 123 57;
              84 115 150]./255;
for i = reshape(intersect(unique(data_chr_all),1:30),1,[]) %onlyl include autosome reshape(unique(Chr),1,[])
    %process SNParray data
    tv = ismember(data_chr_all,i);
%     data_rd = CLImAT_preprocess_rd(data_rd_all(tv));
    data_rd = data_rd_all(tv);
    %remove GC content bias
    data_bc = data_bc_all(tv);
    data_tc = data_tc_all(tv);
%     [data_bc, data_tc] = CLImAT_preprocess_baf(data_bc, data_tc);
    data_baf = data_bc./data_tc;
    data_pos = data_pos_all(tv);
    min_pos = min(data_pos)-100;
    max_pos = max(data_pos)+100;
    
    indx1 = find(chr_seg==i);
    
    %plot
    clf reset
    marker_size = 3;
    fontsize = 15;
    %---plot CN---
    subplot(4,1,1)
    set(gca,'FontSize',12);
    hold on
    set(gca,'YGrid','on','Box','on')
    %     axis ([-Inf Inf -0.05 7.05])
%     axis ([-Inf Inf -0.1 7.1])
    axis ([min_pos max_pos -0.1 7.1])
    set(gca,'YTick',[0:1:7],'XTick',[])
    set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
    ylabel('Copy number','FontSize',fontsize);
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN > 7
            CN = 7;
        end
        line_style = 'r-';
        plot([data_pos(indx(1)) data_pos(indx(end))],[CN CN], line_style, 'LineWidth',3.0);
    end
    barcode_m = barcode;
    tmp = strfind(barcode_m,'_');
    barcode_m(tmp) = '-';
%     title (['Chromosome ' num2str(i)])
    title (['Chromosome ' num2str(i) ', ' barcode_m])
    
    %---plot BAF---
    subplot(4,1,2)
    set(gca,'FontSize',12);
    hold on
    for j=reshape(indx1,1,[])
        CN = cn_seg(j);
        Muc = baf_seg(j);
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 1
            k = 1;
        end
        if CN == 1
            k = 2;
        end
        if CN > 1 && Muc < 1
            k = 3;
        end
        if CN == 2 && Muc == 1
            k = 4;
        end
        if CN > 2 && Muc == 1
            k = 5;
        end
        plot(data_pos(tv),data_baf(tv),'.','MarkerSize',marker_size,'Color',baf_colors(k,:));
    end
%     plot(data_pos,data_baf,'b.', 'MarkerSize',marker_size)
    for j = 0:0.25:1
        plot([data_pos(1) data_pos(end)],[j j],'k-','LineWidth',0.5)
    end
%     %plot expected BAF mean values
%     for j=reshape(indx1,1,[])
%         if state_C1_seg(j) == 21 %total deletion
%             baf_mean = 0.5;
%         else
%             temp1 = w(1)*CN_mapping(state_C1_seg(j))*Muc_mapping(state_C1_seg(j))+(1-w(1));
%             temp2 = w(1)*CN_mapping(state_C1_seg(j))+(1-w(1))*2;
%             baf_mean = temp1/temp2;
%         end
%         plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[baf_mean baf_mean],'r-','LineWidth',1.2);
%         plot([data_pos(start_seg(j)) data_pos(end_seg(j))],[1-baf_mean 1-baf_mean],'r-','LineWidth',1.2);
%     end
    ylabel('BAF','FontSize',fontsize);
    axis ([min_pos max_pos -0.03 1.03])
    set(gca,'YTick',[0 0.5 1],'Box','on')
    set(gca,'XTick',[])
    
    subplot(4,1,3);
    set(gca,'FontSize',12);
    hold on
    for j=reshape(indx1,1,[])
        CN = cn_seg(j);
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 1
            CN = 0;
        end
        k = CN+1;
        if k > 4
            k = 4;
        end
        plot(data_pos(tv),data_rd(tv),'.','MarkerSize',marker_size,'Color', rd_colors(k,:));
        
    end
    ylabel('RD','FontSize',fontsize);
    set(gca,'Box','on')
    set(gca,'XTick',[])
%     axis ([min_pos max_pos min_rd max_rd])
    axis ([min_pos max_pos -Inf Inf])
%     axis([-Inf Inf -Inf Inf])

    subplot(4,1,4);
    set(gca,'FontSize',12);
    hold on
%     score_all = [];
%     score_all_h = [];
    for j=reshape(indx1,1,[])
        line_style = 'r-';
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        plot([data_pos(indx(1)) data_pos(indx(end))],[score_seg(j) score_seg(j)], line_style, 'LineWidth',1.5);
    end
%     axis ([min_pos max_pos 0.001 0.005])
    axis ([min_pos max_pos -10 110])
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0:50:100],'Box','on')
    set(gca,'Box','on')
    ylabel(' Score','FontSize',fontsize);

    %save figure
%     figpath = [plotsdir '\Chr_' num2str(i) '_' barcode];
    figpath = [plotsdir '/Chr_' num2str(i) '_' barcode '.png'];
    %     eval(['print -djpeg -r600 ' figpath ])
    eval(['print -dpng -r400 ' figpath ])

end
