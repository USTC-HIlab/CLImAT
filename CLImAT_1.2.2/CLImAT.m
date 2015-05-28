function CLImAT(Datasource, Outputsource, genotypepath, configfile)
%--------------------------------------------------------------------%
%------------------>       version 1.2.2       <---------------------
%--------------------------------------------------------------------%
% 2014/12/09 by Zhenhua
% This is the version 1.2.2 of the CLImAT method

%------Input and output------%
% Datasource: a directory containing data files to be processed
% Outputsource: directory of output files
% genotypepath: directory of genotype files
% configfile: a file containing user-defined parameters

global current_version
global NoSolutionFlag
current_version = '1.2.2';
sourceflag = 1;

if nargin < 4
    error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m'] );
end

%parameters used in CLImAT
%===============================================
%---for cn up to 7---
depend_table = [...
    %w1
    1 1 0.01 0.5;...
    2 1 1 1.0;...
    3 1 2 0.5;...
    4 1 2 1.0;...
    5 1 3 0.67;...
    6 1 3 1.0;...
    7 1 4 0.75;...
    8 1 4 0.5;...
    9 1 4 1.0;...
    10 1 5 0.8;...
    11 1 5 0.6;...
    12 1 5 1.0;...
    13 1 6 5/6;...
    14 1 6 4/6;...
    15 1 6 0.5;...
    16 1 6 1.0;...
    17 1 7 6/7;...
    18 1 7 5/7;...
    19 1 7 4/7;...
    20 1 7 1.0;...
    %     21 1 0.01 0.5...
    ];

%ABBBBB,AABBBB,AAABBB,BBBBBB
%ABBBBBB,AABBBBB,AAABBBB,BBBBBBB
%===============================================
screening_method = 2;
thres_EM = 1e-4;
max_iter = 30;
verbose = 1;
init_CLImAT_paras = [{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]}]; % initial parameters will be assigned in the main function
                                                                   % parameters:pie,transmat,o,lambda,w,p

%initialization of global variable
global data_rd_sep
global data_baf_sep
global data_pos_sep
global gamma_sep
global condi_probs_sep
global condi_probs_fluct_sep
global tumor_range
global Chromosomes
global rd_median
global gender
global window

% read configuration file
CLImAT_read_config(configfile);

if sourceflag == 1 % reading a directory
    disp(['CLImAT (version ' current_version ') is loading...'])
    datafilelist = dir(Datasource);
    if length(datafilelist) < 3
        error(['No files in the directory ' Datasource]);
    else %now do batch annotation       
        filename = cell(1,(length(datafilelist)-2));
        for i = 3:length(datafilelist)
            filename{i-2} = datafilelist(i).name;
        end
        tumor_range_table = [0.1*ones(length(filename),1) ones(length(filename),1)];
        
        %--------------perform CLImAT--------------------      
        %record all the time used for batch annotation
        t_all = 0;
        if length(filename) == 1
            disp('-----CLImAT batch annotation starts now, ONE file is found-----');
        elseif length(filename)>1
            disp(['-----CLImAT batch annotation starts now, ' num2str(length(filename)) ' files are found-----']);
        end
        for fileindx = 1:length(filename)
            %clear global variables            
            gamma_sep = [];
            condi_probs_sep = [];
            condi_probs_fluct_sep = [];
            
            %record time cost
            tic
            tumor_range = tumor_range_table(fileindx,:);
            results = regexp(filename{fileindx}, '^(.+)\.+.+','tokens', 'once');
            if isempty(results)
                fn_nosuffix = filename{fileindx};
            else
                fn_nosuffix = results{1};
                if ~isempty(strfind(fn_nosuffix,'.'))
                    fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
                end
            end
            
            %--------open result files --------------
            o_fid = fopen([Outputsource '/' fn_nosuffix '.results'],'w');
            if o_fid == -1
                error(['Can not open result file for ' filename{fileindx}]);
            end
            
            %--------open genotype files --------------
            fpG = fopen([genotypepath '/' fn_nosuffix '.Gtype'],'w');
            if fpG == -1
                error(['Can not open genotype result file ' genotypepath '/' fn_nosuffix '.Gtype']);
            end
            
            %--------------load data--------------------
            [data_chr_all, data_pos_all, data_bc_all, data_tc_all, data_rd_all] = CLImAT_load_preprocessData([Datasource '/' filename{fileindx}]);

            rd_median = median(data_rd_all);
            Chromosomes = reshape(unique(data_chr_all),1,[]);
            
            nfilename = [Outputsource '/' fn_nosuffix '_normalized.mat'];
            
            eval(['save ' nfilename ' data_chr_all data_pos_all data_bc_all data_tc_all data_rd_all']);
            
            if size(data_pos_all,1)>size(data_pos_all,2) %make sure it's 1byN
                data_pos_all = data_pos_all';
            end
            if size(data_bc_all,1)>size(data_bc_all,2) %make sure it's 1byN
                data_bc_all = data_bc_all';
            end
            if size(data_tc_all,1)>size(data_tc_all,2) %make sure it's 1byN
                data_tc_all = data_tc_all';
            end
            if size(data_rd_all,1)>size(data_rd_all,2) %make sure it's 1byN
                data_rd_all = data_rd_all';
            end
            
            %use at least 30000 SNPs for screening
            stepsize_ds = max(floor(length(data_chr_all)/30000),1);
            
            %-------divide into different chromosomes----------
            chr_num = length(Chromosomes);
            data_pos_sep = cell(chr_num,1);
            data_bafA_sep = cell(chr_num,1);
            data_rd_sep = cell(chr_num,1);
            data_baf_sep = cell(chr_num,1);
            data_het_sep = cell(chr_num,1);
                        
            for i = 1:chr_num
                tv = ismember(data_chr_all,Chromosomes(i));
                % RD
                data_rd_sep{i} = data_rd_all(tv);
                % B alellic Count
                data_bc = data_bc_all(tv);
                data_tc = data_tc_all(tv);
                data_baf = data_bc./data_tc;               
                tvA = data_baf < 0.5;
                data_bafA_sep{i} = tvA;             
                data_bc(tvA) = data_tc(tvA)-data_bc(tvA);      
                data_baf_sep{i} = [data_bc; data_tc];
                clear tvA data_baf data_bc data_tc;
                % Position
                data_pos_sep{i} = data_pos_all(tv);
                % het
                data_het_sep{i} = zeros(1,length(data_rd_all(tv))) > 0;
                clear tv
            end
            clear data_rd_all data_bc_all data_tc_all data_pos_all data_chr_all;
            
            %------------------ call CLImAT --------------------
            [temp1,CLImAT_paras,best_indx] = ...
                CLImAT_main(screening_method,init_CLImAT_paras,depend_table,stepsize_ds,thres_EM,max_iter,verbose,fn_nosuffix);
            
            %-------------- process results ----------------------
            w = CLImAT_paras{5}{best_indx};
            lambda = CLImAT_paras{4}{best_indx};
            o = CLImAT_paras{3}{best_indx};
            p = CLImAT_paras{6}{best_indx};
           
            [p_states,num_SNP,aCN,p_s,states_pre] = CLImAT_process_results_new(depend_table);
            
            cn = [0 1 2 2 3 3 4 4 4 5 5 5 6 6 6 6 7 7 7 7 ];
            Bcn = [0 1 1 2 2 3 3 2 4 4 3 5 5 4 3 6 6 5 4 7 ];
            AI_mapping = [1 1 3 2 3 2 3 3 2 3 3 2 3 3 3 2 3 3 3 2 1];%DEL(1) LOH(2) HET(3)
            
            states_pre = states_pre{best_indx};
                       
            cn_segs_all = zeros(size(states_pre,1),1); % copy number
            mcn_segs_all = zeros(size(states_pre,1),1); % major copy number 
            AI_segs_all = zeros(size(states_pre,1),1); % allelic imblance status
            scores = zeros(size(states_pre,1),1); % reliability socre
            
            for i = 1:size(states_pre,1)
                Chr_i = states_pre(i,1);
                St_i = states_pre(i,2);
                Ed_i = states_pre(i,3);
                S_C1 = states_pre(i,4);
                data_rd = data_rd_sep{Chr_i}(St_i:Ed_i);
                data_baf = data_baf_sep{Chr_i}(:,St_i:Ed_i);

                % assign copy number and allelic imblance status to each segment
                if cn(S_C1) >= 6
                    rd_mean = mean(data_rd_sep{Chr_i}(St_i:Ed_i));
                    if gender == 1 && (Chromosomes(Chr_i) == 23 || Chromosomes(Chr_i) == 24)
                        real_cn = round(((rd_mean-o)*2/lambda-w)/(1-w));
                    else
                        real_cn = round(((rd_mean-o)/lambda-w)*2/(1-w));
                    end
                    if real_cn < 0
                        cn_segs_all(i) = 0;
                    else
                        cn_segs_all(i) = real_cn;
                    end
                    if gender == 1 && (Chromosomes(Chr_i) == 23 || Chromosomes(Chr_i) == 24)
                        mcn_segs_all(i) = cn_segs_all(i);
                        AI_segs_all(i) = -1;
                    else
                        if cn_segs_all(i) < 2
                            mcn_segs_all(i) = cn_segs_all(i);
                            AI_segs_all(i) = 1;
                        elseif cn_segs_all(i) == cn(S_C1)
                            mcn_segs_all(i) = Bcn(S_C1);
                            AI_segs_all(i) = AI_mapping(S_C1);
                        else
                            mcn_segs_all(i) = CLImAT_majorCN_estimate(data_baf, cn_segs_all(i), w);
                            if mcn_segs_all(i) == cn_segs_all(i)
                                AI_segs_all(i) = 2;
                            else
                                AI_segs_all(i) = 3;
                            end
                        end
                    end
%                     if cn_segs_all(i) < 2
%                         mcn_segs_all(i) = cn_segs_all(i);
%                         AI_segs_all(i) = 1;
%                     elseif cn_segs_all(i) == cn(S_C1)
%                         mcn_segs_all(i) = Bcn(S_C1);
%                         AI_segs_all(i) = AI_mapping(S_C1);
%                     else
%                         mcn_segs_all(i) = CLImAT_majorCN_estimate(data_baf, cn_segs_all(i), w);
%                         if mcn_segs_all(i) == cn_segs_all(i)
%                             AI_segs_all(i) = 2;
%                         else
%                             AI_segs_all(i) = 3;
%                         end
%                     end
                else
                    if gender == 1 && (Chromosomes(Chr_i) == 23 || Chromosomes(Chr_i) == 24)
                        cn_segs_all(i) = cn(S_C1);
                        mcn_segs_all(i) = cn(S_C1);
                        AI_segs_all(i) = -1;
                    else
                        cn_segs_all(i) = cn(S_C1);
                        mcn_segs_all(i) = Bcn(S_C1);
                        AI_segs_all(i) = AI_mapping(S_C1);
                    end
                end

                data_het_sep{Chr_i}(St_i:Ed_i) = CLImAT_het_estimate(Chr_i, data_baf, cn_segs_all(i), mcn_segs_all(i), w);
                
                scores(i) = CLImAT_get_reliabilityScore(Chr_i,data_baf,data_rd,data_het_sep{Chr_i}(St_i:Ed_i),o,lambda,w,p,cn_segs_all(i),mcn_segs_all(i));

            end
            
            len = sum(states_pre(:,3)-states_pre(:,2)+1);
            m_score = sum(scores.*(states_pre(:,3)-states_pre(:,2)))/len;
            scores = Scale_score(scores,m_score);
            
            pre_epos = 0;
            pre_chr = 0;
            bp_len = 0;
            cn_w = 0;
            for i = 1:size(states_pre,1)
                Chr_i = states_pre(i,1);
                St_i = states_pre(i,2);
                Ed_i = states_pre(i,3);
                St_pos = data_pos_sep{Chr_i}(St_i);
                Ed_pos = data_pos_sep{Chr_i}(Ed_i);
                if Chromosomes(Chr_i) ~= pre_chr
                    pre_chr = Chromosomes(Chr_i);
                    pre_epos = St_pos-1;
                end
                if St_pos ~= pre_epos+1
                    St_pos = pre_epos+1;
                end
                pre_epos = Ed_pos;
                bp_len = bp_len+(Ed_pos-St_pos+1);
                cn_w = cn_w+(Ed_pos-St_pos+1)*cn_segs_all(i);
            end
            
            ACN = cn_w/bp_len;
            
            listfile = 'LOG.txt';
            fp_list = fopen(listfile,'a+');
            if fp_list == -1
                warning('Can not open the file: "%s"!\n',listfile);
            else
                if fileindx == 1
                    fprintf(fp_list,'Version:\tDate\tTime\tSample\tTumor purity\tTumor ploidy\tCopy neutral RD\tShift of RD\n');
                end
                fprintf(fp_list,'CLImAT%s\t%s\t%s\t%s\t%f\t%f\t%6.0f\t%6.0f\n',...
                    current_version,datestr(clock,'mmm-dd-yyyy'),datestr(clock,'HH:MM:SS'),fn_nosuffix,1-w,ACN,lambda,o);
            end
            fclose(fp_list);
            
            % save genotype of each SNP               
            fprintf(fpG,'Chr\tposition\tCN\tB allele CN\tGenotype\n');
            for i = 1:length(Chromosomes)
                tv_chr = states_pre(:,1) == i;
                cn_segs = cn_segs_all(tv_chr);
                mcn_segs = mcn_segs_all(tv_chr);
                St_i = states_pre(tv_chr,2);
                Ed_i = states_pre(tv_chr,3);
                Bcn_status = zeros(length(data_pos_sep{i}),1);
                tv = (data_bafA_sep{i}==1 & data_het_sep{i});
                Bcn_status(tv) = 1;
                tv = (data_bafA_sep{i}==0 & ~data_het_sep{i});
                Bcn_status(tv) = 2;
                tv = (data_bafA_sep{i}==1 & ~data_het_sep{i});
                Bcn_status(tv) = 3;
                for j = 1:length(data_pos_sep{i})
                    tv = j >= St_i & j <= Ed_i;
                    if Bcn_status(j) == 0
                        b_cn = mcn_segs(tv);
                    elseif Bcn_status(j) == 1
                        b_cn = cn_segs(tv)-mcn_segs(tv);
                    elseif Bcn_status(j) == 2
                        b_cn = cn_segs(tv);
                    else
                        b_cn = 0;
                    end
                    fprintf(fpG,'%d\t%d\t%d\t%d',Chromosomes(i),data_pos_sep{i}(j),cn_segs(tv),b_cn);
                    if cn_segs(tv) == 0
                        fprintf(fpG,'\t%s\n','N/A');
                    else
                        fprintf(fpG,'\t%s\n',[repmat('A',1,cn_segs(tv)-b_cn) repmat('B',1,b_cn)]);
                    end 
                end
                clear tv_chr cn_seq mcn_segs St_i Ed_i
            end
            fclose(fpG);
            
            %-------------- output summary of the results--------------
            fprintf(o_fid,'---------------------------------------------------------------\n');
            fprintf(o_fid,['             Summary of CLImAT results (version ' current_version ')          \n']);
            if NoSolutionFlag
                fprintf(o_fid,'Warning: Prediction results may be inaccurate due to the failure\n');
                fprintf(o_fid,'in finding optimal initial global parameters!\n');
            end
            
            fprintf(o_fid,'General information of this cancer sample:                      \n');
            fprintf(o_fid,'   Copy neutral read depth: %6.0f\n',lambda);
            fprintf(o_fid,'   Shift of read depth: %6.0f\n',o);
            fprintf(o_fid,'   Proportion of abnormal cells in the sample: %6.4f\n',1-w);
            fprintf(o_fid,'   Average copy number: %1.2f\n',ACN);
            fprintf(o_fid,'   Parameter p of NB distributions:');
            for i = 1:length(p)
                fprintf(o_fid,' %6.4f',p(i));
            end
            fprintf(o_fid,'\n');
            fprintf(o_fid,'   Proportion of all abnormal chromosomal regions: %6.4f\n',1-p_s(4));
            fprintf(o_fid,'   Estimated average cancer DNA index: %6.4f\n',ACN/2);
            fprintf(o_fid,'---------------------------------------------------------------\n');
            fprintf(o_fid,'\n');            
            
            %--------------output state assignment in segments--------------
            fprintf(o_fid,'Chr\tStartPos\tEndPos\tCN\tBAF\tAI\tScore\tLength\n');
            for i = 1:size(states_pre,1)
                Chr_i = states_pre(i,1);
                St_i = states_pre(i,2);
                Ed_i = states_pre(i,3);
                St_pos = data_pos_sep{Chr_i}(St_i);
                Ed_pos = data_pos_sep{Chr_i}(Ed_i);

                if cn_segs_all(i) == 0
                    fprintf(o_fid,'%d\t%d\t%d\t%d\t%1.3f\t%d\t%3.1f\t%d\n',Chromosomes(Chr_i),...
                        St_pos,Ed_pos,cn_segs_all(i),0.5,AI_segs_all(i),scores(i),Ed_pos-St_pos+1);
                else
                    fprintf(o_fid,'%d\t%d\t%d\t%d\t%1.3f\t%d\t%3.1f\t%d\n',Chromosomes(Chr_i),...
                        St_pos,Ed_pos,cn_segs_all(i),mcn_segs_all(i)/cn_segs_all(i),AI_segs_all(i),scores(i),Ed_pos-St_pos+1);
                end

            end
            
            fclose(o_fid);
            
            clear CLImAT_paras best_indx data_id_sep
            %6.4finall display a brief report
            t = toc;
            t_all = t_all+t;
            disp ([num2str(fileindx) '. ' filename{fileindx} ' is done, time used: ' num2str(t)]);
            
        end %6.4for fileindx = 3:length(datafilelist)
    end % if length(datafilelist)<3
    
    disp(['-----CLImAT is finished, totally ' num2str(t_all/60) ' minites were used-----']);
    clear all
end %if sourceflag  == 1

clear all
close all
