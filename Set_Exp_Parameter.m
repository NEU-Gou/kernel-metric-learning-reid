% Set the experiment parameters for the test script
%% load feature
% flag_preFilter = 0;
% AlgoOption.isLDFV = 0;
partition_name = 'Random';  %'Random'; %'SDALF'; % cam3ref

% Featurename = [dataset_name '_HistMoment' num2str(num_patch) 'Patch'];
% if flag_preFilter
%     fname = [dataset_name '_HistMoment' num2str(num_patch) 'Patch_wPreFiltering.mat'];
% else
%     fname = [dataset_name '_HistMoment' num2str(num_patch) 'Patch_woPreFiltering.mat'];
% end
% switch featureType
%     case 'HistLBP'
%         fname = [dataset_name '_Hist' num2str(num_patch) 'Patch'];
%     case 'CLBP'
%     case 'LDFV'
%         fname = [dataset_name '_LDFV' num2str(num_patch) 'Patch'];
%         AlgoOption.isLDFV = 1;
%         FeatureAppearance = []; 
% end

fname = [dataset_name '_HistLBP' num2str(num_patch) 'Patch'];
try
    load([dropbox_folder '/feature/' fname]);
catch 
    error('No such type of feature precomputed.');
end
%% picking the useful features

% New appearance feature format
if num_patch == 341 && ~strcmp(algoname ,'oLFDA') &&...
        ~strcmp(algoname ,'svmml') && ~strcmp(algoname ,'KISSME')
    Feature = sparse(double(Feature));
else
    Feature = single(Feature);
end

AlgoOption.doPCA = 0;  % flag for PCA preprocessing
if strcmp(algoname ,'svmml')% apply PCA for svmml
    [COEFF,pc,latent,tsquare] = princomp(Feature,'econ');
    pcadim =  sum(cumsum(latent)/sum(latent)<0.95); %80;%
    Feature = pc(:, 1:pcadim);
    AlgoOption.doPCA = 1;
end
if strcmp(algoname ,'KISSME')
    [COEFF,pc,latent,tsquare] = princomp(Feature,'econ');
    %         pcadim =  sum(cumsum(latent)/sum(latent)<0.95); %80;%
    Feature = pc(:, 1:pcadim);
    AlgoOption.doPCA = 1;
end
% else
%     fname = [dataset_name '_Hist' num2str(num_patch) 'Patch'];
%     load([dropbox_folder '/Feature/' fname]);
%     AlgoOption.doPCA = 0; 
%     Feature = FeatureAppearence;
% end
% clear DataSet;
%% load dataset partition
load([dropbox_folder '/feature/' dataset_name '_Partition_' partition_name '.mat']);
load([dropbox_folder '/dataset/' dataset_name '_Images.mat'], 'gID', 'camID');

%%
% The number of test times with the same train/test partition.
% In each test, the gallery and prob set partion is randomly divided.
num_itr =10; 
np_ratio =10; % The ratio of number of negative and positive pairs. Used in PCCA
% default algorithm option setting
AlgoOption.name = algoname;
AlgoOption.func = algoname; % note 'rPCCA' use PCCA function also.
AlgoOption.npratio = np_ratio; % negative to positive pair ratio
AlgoOption.beta =3;  % different algorithm have different meaning, refer to PCCA and LFDA paper.
AlgoOption.d =40; % projection dimension
AlgoOption.epsilon =1e-4;
AlgoOption.lambda =0;
AlgoOption.w = [];
AlgoOption.dataname = fname;
AlgoOption.partitionname = partition_name;
AlgoOption.num_itr=num_itr;
% customize in different case
switch  algoname
    case {'LFDA'}
        AlgoOption.npratio =0; % npratio is not required.
        AlgoOption.beta =0.01;
        AlgoOption.d =40;
        AlgoOption.LocalScalingNeighbor =6; % local scaling affinity matrix parameter.
        AlgoOption.num_itr= 10;
    case {'oLFDA'}
        AlgoOption.npratio =0; % npratio is not required.
        AlgoOption.beta =0.15; % regularization parameter
        AlgoOption.d = 40;
        AlgoOption.LocalScalingNeighbor =6; % local scaling affinity matrix parameter.
        AlgoOption.num_itr= 10;
    case {'rPCCA'}
        AlgoOption.func = 'PCCA';
        AlgoOption.lambda =0.01;
    case {'svmml'}
        AlgoOption.p = []; % learn full rank projection matrix
        AlgoOption.lambda1 = 1e-8;
        AlgoOption.lambda2 = 1e-6;
        AlgoOption.maxit = 300;
        AlgoOption.verbose = 1;
    case {'MFA'}
        AlgoOption.Nw = 0; % 0--use all within class samples
        AlgoOption.Nb = 12;
        AlgoOption.d = 30;
        AlgoOption.beta = 0.01;
%     case {'PRDC'} % To be added in the future
%         AlgoOption.Maxloop = 100;
%         AlgoOption.Dimension = 1000;
%         AlgoOption.npratio = 0;
    case {'KISSME'}
        AlgoOption.PCAdim = pcadim;
        AlgoOption.npratio = 10;
        AlgoOption.nFold = 20;
end

% if strfind(fname, 'COV')
%     AlgoOption.isCOV = 1;
%     kname={'COV'};
% else
%     AlgoOption.isCOV = 0;
if strcmp(algoname ,'oLFDA')|| strcmp(algoname ,'svmml') || strcmp(algoname ,'KISSME')
    kname = {'linear'};
end
% end

