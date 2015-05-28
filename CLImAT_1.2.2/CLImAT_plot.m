function CLImAT_plot(Datasource, resultssource, plotssource)
%this function is used to plots results of multiple samples.
%Datasource: the directory containing all data files
%resultssource: the directory containing all result files (one to one
%mapping assumed between data file and result file
%plotssource: the directory that plots will be saved
if isdeployed
    if nargin < 3
        error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.txt.'] );
    end
else
    if nargin < 3
        error(['Insufficient input parameters, Please check again! ' ...
            'More details in example.m.'] );
    end
end
current_version = '1.2.2';
disp(['CLImAT_plot (version ' current_version ') is loading...'])
tic
% if nargin == 4
% plot_mc = 0;
% end

datafilelist = dir(Datasource);

if length(datafilelist) < 3
    error('No data files in the directory!');
else %now do batch annotation
    if length(datafilelist) > 3
        disp(['-----Plot CLImAT results now, ' num2str(length(datafilelist)-2) ' files are found-----'])

    else
        disp('-----Plot CLImAT results now, ONE file is found-----')
    end
    filename = cell(1,(length(datafilelist)-2));
    for i = 3:length(datafilelist)
        filename{i-2} = datafilelist(i).name;
    end

    for fileindx = 1:length(filename)
        results = regexp(filename{fileindx},'^(.+)\.+.+','tokens','once');
        if isempty(results)
            fn_nosuffix = filename{fileindx};
        else
            fn_nosuffix = results{1};
            if ~isempty(strfind(fn_nosuffix,'.'))
                fn_nosuffix(strfind(fn_nosuffix,'.')) = '_';
            end
        end
        Datafile = [resultssource '/' fn_nosuffix '_normalized.mat'];
        resultsfile = [resultssource '/' fn_nosuffix '.results'];  
        plotsdir = [plotssource '/' fn_nosuffix];
        s = mkdir(plotsdir);
        if ~s %failed to make a directory
            error(['Can not make directory: ' plotsdir]);
        else
%             CLImAT_plot_normalized_results_new(Datafile,resultsfile,plotsdir,fn_nosuffix);
            CLImAT_plot_normalized_results(Datafile,resultsfile,plotsdir,fn_nosuffix);
            disp ([num2str(fileindx) '. ' filename{fileindx} ' is done'])
        end
    end %for fileindx = 3:length(datafilelist)
end %if length(datafilelist)<3

t_all = toc;
disp(['-----Batch plotting is finished, totally ' num2str(t_all/60) ' minites were used-----'])
clear all
close all