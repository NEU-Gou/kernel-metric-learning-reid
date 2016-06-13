% Script Showing Result
% clc
% clear
% close all;
patition_folder = [dropbox_folder '/Feature/'];
result_folder = pwd;
%%
% figure
% clf;
num_patch = [6];% 341 6 14 75
flag_preFilter = 0;
partition_name = 'Random';  %'Random'; %'SDALF'; %
algoname = {'PCCA','oLFDA','svmml','KISSME','rPCCA','LFDA','MFA'}; % 'PCCA','oLFDA','svmml','KISSME','rPCCA','LFDA','MFA'
kname={'linear' 'chi2' 'chi2-rbf'}; %  'linear','chi2' 'chi2-rbf'
color={ 'c','g', 'b','r','k','m','y'};
linesty ={'o-.','d--','*-'};
clear M Method
clear result_name
% cnt=1;
datasetname = {'iLIDS'}; %'VIPeR' 'iLIDS' 'CAVIAR' '3DPeS';
rrrr=[];
for kkk=1: length(datasetname)
    dataset_name = datasetname{kkk};
    load([patition_folder dataset_name '_Partition_' partition_name '.mat']);
    num_ref = sum(Partition(1).ix_test_gallery(1,:));
    h(kkk) = figure();
    cnt=1;
    
    for ixp =1:length(num_patch)
        %     figure;
        np =num_patch(ixp);
        
        for m = 1: length(algoname)
            for k =1:length(kname)
                kname{k};
                Featurename = [dataset_name '_HistLBP' num2str(np) 'Patch'];
                Featurename = ['Result_' dataset_name '_' algoname{m} '_' Featurename '_' partition_name '_' kname{k} '_'];
                fname = dir([result_folder '/' Featurename '*']);
                if size(fname,1) >1 
                    fname = fname(2);
                end
                if isempty(fname)
                    continue;
                end
                display(fname.name);
                fname = [result_folder '/' fname.name];
                fname(strfind(fname,' ')) =[];
                temp = strfind(fname,'.mat');
                load(fname(1:temp+3));
                M(ixp,m, k ,:,:) = Method;
                prob = zeros(size(Method,1), 10000);
                for i =1: size(Method,1)
                    r = Method{i,2}.Ranking(1:10,:); % 2nd col storing test result
                    [a, b] = hist(r' ,1: size(r,2));
                    a = a./repmat(ones(1,size(a,2))*size(r,2), size(a,1),1);
                    prob(i, 1:size(a,1)) =mean(a');
                end
                prob = mean(cumsum(prob,2),1);
                temp = [color{m} linesty{k}];
                
                if strcmp(algoname{m},'oLFDA')
                    result_name{cnt} = 'LFDA';
                elseif strcmp(algoname{m},'svmml')
                    result_name{cnt} = 'svmml';
                else
                    switch kname{k}
                        case {'linear'}
                            result_name{cnt} = [algoname{m} '-' 'L'];
                            if strcmp(algoname{m},'LFDA')
                                result_name{cnt} = ['k' algoname{m} '-' 'L'];
                            else
                                result_name{cnt} = [algoname{m} '-' 'L'];
                            end
                        case {'chi2'}
                            if strcmp(algoname{m},'LFDA')
                                result_name{cnt} = ['k' algoname{m} '-' '\chi^2'];
                            else
                                result_name{cnt} = [algoname{m} '-' '\chi^2'];
                            end
                        case {'chi2-rbf'}
                            if strcmp(algoname{m},'LFDA')
                                result_name{cnt} = ['k' algoname{m} '-' 'R{\chi^2}'];
                            else
                                result_name{cnt} = [algoname{m} '-' 'R{\chi^2}'];
                            end
                    end
                end
                if length(num_patch)>1
                    result_name{cnt} = [result_name{cnt} '-' num2str(np)];
                end
                plot(prob(1:25)*100,temp,'linewidth',1.2,'markersize',10); hold on;
                prob= ([ prob([1 5 10 20])  CalculatePUR(prob(1:size(Method{1,2}.Ranking,2)),num_ref)] )*100;
                rrr(cnt,:)=prob;
                display([dataset_name '_HistMoment' num2str(np) 'Patch' algoname{m} ' Rank 1 5 10 20 accuracy and PUR===>']);
                display(num2str(prob,'%.1f  '));
                cnt = cnt +1;
            end
        end
    end
 
    legend(result_name);
    grid on;
    xlabel('Rank Score','Fontsize', 18, 'fontweight', 'bold')
    ylabel('Matching Rate (%)','Fontsize', 18, 'fontweight', 'bold')
    title(['CMC curve on ' dataset_name ' with ' num2str(np) ' patch'],'Fontsize', 18, 'fontweight', 'bold');
%     rrrr =[rrrr; rrr([1:3 10 4:9], :)]
    rrrr =[rrrr; rrr];
end
M=squeeze(M);
matrix2latex_tab([rrrr]','%.1f  ');
