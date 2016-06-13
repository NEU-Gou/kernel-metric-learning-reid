clc
clear
warning off
dropbox_folder =  pwd;%'C:\Users\Goo\Dropbox\Fei_Mengran\Matlab Code\ReID';
addpath(genpath('Assistant Code'));
mkdir('feature');
%%
% name of metric learning algorithm 
algoname = 'PCCA'; %'oLFDA'; 'PCCA'; 'rPCCA'; 'LFDA'; 'MFA'; 'KISSME'; 'svmml' 
% dataset name
dataset_name = 'iLIDS'; %{'VIPeR' 'iLIDS' 'CAVIAR' '3DPeS'};
% kernel name
kname={'linear'}; % {'linear', 'chi2', 'chi2-rbf'};
% % % % % % Feature types
% % % % % featureType = 'LDFV';
% number of patches(stripes)
num_patch = 6; %6, 14, 75, 341
% PCA dimension, ONLY used in KISSME
pcadim = 65; %[77 45 65 70];
        
% 
if ~exist(['feature/' dataset_name '_Partition_Random.mat'],'file')
    Set_Partition(dataset_name);
end
if ~exist(['feature/' dataset_name '_HistLBP' num2str(num_patch) 'Patch.mat'],'file')
    steps =flipud([4, 4; 8 8; 16 16; 48 21]); % set the moving step size of the region.
    BBoxszs =flipud([8,8; 16 16; 32 32; 48 22]); % set the region size.
    numP = [6 14 75 341];
    
    load(fullfile('dataset',[dataset_name, '_Images.mat']));
    % image size check
    imgsizes = cellfun(@size, I, 'uni',0);
    imgsizes = cell2mat(imgsizes');
    if any(imgsizes(:,1)~=imgsizes(1,1)) || any(imgsizes(:,2)~=imgsizes(1,2))
        I = cellfun(@(x) imresize(x,[128 48]),I,'uni',0);
    end
    Feature = HistLBP(I,16,BBoxszs(numP==num_patch,:),steps(numP==num_patch,:));
    save(fullfile('feature',[dataset_name,'_HistLBP',num2str(num_patch),'Patch.mat']),'Feature','-v7.3');
end

test_Ranking
Script_demo_result;