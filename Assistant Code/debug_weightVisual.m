% visualization of the weight of features
clear 
clc
dropbox_folder = 'C:\Users\Mengran\Dropbox\Fei_Mengran\Matlab Code\ReID';
result_folder = 'C:\Users\Mengran\Research\code\ReID\EnsembleReID_v1.3\result\ECCV14\';
%%
partition_name = 'Random';
ALGO = {'PCCA','LFDA','MFA','rPCCA'}; %'LFDA','MFA','rPCCA'
DATA = {'VIPeR','iLIDS','CAVIAR','3DPeS'}; %'VIPeR','iLIDS','CAVIAR',
num_patch = 341;
[region_idx, BBox] =GenerateGridBBox([128 48], [8 8], [4 4]);


for aaa = 1:numel(ALGO)
    for ddd = 1:numel(DATA)
        algoname = ALGO{aaa};
        dataset_name = DATA{ddd};
        Set_Exp_Parameter;
        kname = 'linear';
        fname = ls([result_folder 'Result_' DATA{ddd} '_' ALGO{aaa} '_*' num2str(num_patch) 'Patch_' partition_name '_' kname '_*']);
        if isempty(fname)
            continue;
        end
        load([result_folder fname]);

        projM = Method{1,2}.P;
        idx_train = Partition(1).idx_train ;
        train = Feature(idx_train,:);
        weightP = projM*train;
        % weightP = weightP(1,:);
        weightP = sum(weightP.^2);
        weightP = weightP.^0.5;

        num_bin = size(weightP,2)/num_patch;
        wPatch = zeros(num_bin,num_patch);
        % back projection FEATURE-->PATCH
        fname = [DATA{ddd} '_HistMoment' num2str(num_patch) 'Patch_woPreFiltering.mat'];
        load(fname);
        featurename = fieldnames(DataSet.idxfeature);

        % color histogram
        cbin = 16;
        for i = 1:8 %length(featurename)
            temp = reshape(weightP(cbin*num_patch*(i-1)+1:cbin*num_patch*i),cbin,num_patch);    
            wPatch(cbin*(i-1)+1:cbin*i,:) = temp;
        end

        % LBP histogram
        Lbin = [59 243];
        for i = 1
            temp = reshape(weightP(cbin*num_patch*8+1:cbin*num_patch*8+Lbin(i)*num_patch),...
                Lbin(i),num_patch);
            wPatch(cbin*8+1:cbin*8+Lbin(i),:) = temp;
        end

        % Visualization the weight matrix
        wPatch = sum((wPatch).^2);
        wIM = reshape(wPatch,31,11);
        
        % Normalize the weight matrix
%         wIM = 50*wIM./sum(sum(wIM));
        wIM = wIM./max(max(wIM));
        h = figure;
        imagesc(wIM,[0,1]);
        axis equal
        xlim([1,11])
        ylim([1,31])
        title([ALGO{aaa}, ' on ', DATA{ddd}]);
        saveas(h,[ALGO{aaa}, ' on ', DATA{ddd} '.pdf']);
        saveas(h,[ALGO{aaa}, ' on ', DATA{ddd} '.png']);
        close(h);
    end
end


