% Script Showing Result
clc
clear
close all;
dropbox_folder = pwd;
result_folder = pwd;
partition_folder = 'Feature/';
%%
% figure
% clf;
num_patch =[6 14 75];%14 75 341
flag_preFilter = 0;
partition_name = 'Random';  %'Random'; %'SDALF'; %
algoname = {'PCCA' 'oLFDA' 'svmml' 'KISSME' 'rPCCA' 'LFDA' 'MFA'}; % 'LFDA','rPCCA','MFA','PCCA','oLFDA','svmml'
kname={'linear' 'chi2' 'chi2-rbf'}; %  'k',
color={ 'c','g', 'b','r','k','m'};
linesty ={'o-.','d--','*-'};
clear M Method
clear result_name
% cnt=1;
datasetname = {'VIPeR','iLIDS','CAVIAR','3DPeS'}; % 'VIPeR','iLIDS','CAVIAR','3DPeS'
rrrr=[];
num_itr = 10;
for kkk=1: length(datasetname)
    clear M
    dataset_name = datasetname{kkk};
    load([partition_folder dataset_name '_Partition_' partition_name '.mat']);
    num_ref = sum(Partition(1).ix_test_gallery(1,:));
    h(kkk) = figure();
    cnt=1;
    for ixp =1:length(num_patch)
        %     figure;
        np =num_patch(ixp);
        
        for m = 1: length(algoname)
            for k =1:length(kname)
                kname{k};
                Featurename = [dataset_name '_HistMoment' num2str(np) 'Patch'];
                Featurename = ['Result_' dataset_name '_' algoname{m} '_' Featurename '_' partition_name '_' kname{k} '_'];
                fname = dir([result_folder '/' Featurename '*']);
%                 ixst =strfind(fname, 'Result');
%                 if length(ixst) >1
%                     fname = fname(ixst(2):end);
%                 end
                if size(fname,1) >1 
                    fname = fname(2,:);
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
                    r = Method{i,2}.Ranking(1:num_itr,:);
                    [a, b] = hist(r' ,1: size(r,2));
                    a = a./repmat(ones(1,size(a,2))*size(r,2), size(a,1),1);
                    if size(a,1) ~= 1
                        a = a';
                    end
                    prob(i, 1:size(a,2)) =mean(a,1);
                end
                prob = mean(cumsum(prob,2),1);
                %         [a, b] = hist(r',1:316);
                %         a = cumsum(a)./repmat(ones(1,size(a,2))*316, size(a,1),1);
                %         a = a';
                %         a =mean(a);
                temp = [color{mod(ceil(cnt/length(kname)), length(color))+1} linesty{k}];
                
                % form lengends.
                switch kname{k}
                    case {'linear'}
                        result_name{cnt} = [algoname{m} '-' 'L'];
                        if strcmp(algoname{m},'LFDA')
                            result_name{cnt} =[algoname{m}];
                            temp = ['k'  linesty{k}];
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
                if length(num_patch)>1
                    result_name{cnt} = [result_name{cnt} '-' num2str(np)];
                end
%                 plot(prob(1:25)*100,temp,'linewidth',1.2,'markersize',10); hold on;
                Prob(cnt,:)= prob;
                prob= ([ prob([1 5 10 20])  CalculatePUR(prob(1:size(Method{1,2}.Ranking,2)), num_ref)] )*100;
                rrr(cnt,:)=prob;
                
%                 display([dataset_name '_HistMoment' num2str(np) 'Patch' algoname{m} ' Rank 1 5 10 20 accuracy ===>']);
%                 display(num2str(prob,'%.1f  '));
                cnt = cnt +1;
            end
        end
    end
%     M=squeeze(M);
    %%
    [~,iix]= max(rrr(:,5));
    plot(Prob(iix, 1:25)*100,temp,'linewidth',1.2,'markersize',10); hold on;
    result_name = result_name(iix);
    rrr = rrr(iix,:);
    load([partition_folder dataset_name '_Partition_' partition_name '.mat']);
    load([dropbox_folder '/dataset/' dataset_name '_Images.mat'], 'gID', 'camID')
    [m] = Vote_MultiRankingAlgo(M, Partition,gID,1);
    m1 = cumsum(mean(reshape(m, [],size(m,3)),1));
    rrr= [rrr;( [m1([1 5 10 20])  CalculatePUR(m1, num_ref)] )*100];
    plot(m1(1:25)*100,'c+-','linewidth',1.2,'markersize',10);
    result_name = [result_name, {'Ensemble 1'}];
    
    [m] = Vote_MultiRankingAlgo(M, Partition,gID,0);
    m2 = cumsum(mean(reshape(m, [],size(m,3)),1));
    rrr= [rrr;( [m2([1 5 10 20])  CalculatePUR(m2, num_ref)] )*100];
    plot(m2(:,1:25)*100,'m+-','linewidth',1.2,'markersize',10);
    result_name = [result_name, {'Ensemble 2'}];
%     r2 = reshape(r, [],size(r,3));
    legend(result_name);
    grid on;
    xlabel('Rank Score','Fontsize', 18, 'fontweight', 'bold')
    ylabel('Matching Rate (%)','Fontsize', 18, 'fontweight', 'bold')
%     title(['CMC curve on ' dataset_name ' of ensemble' num2str(np) ' patch'],'Fontsize', 18, 'fontweight', 'bold');
    rrrr =[rrrr; rrr];
    
end

matrix2latex_tab([rrrr]','%.1f  ');


