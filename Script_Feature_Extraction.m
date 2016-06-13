%% set the parameter and load the dataset mat file
% dataset_name = 'CAVIAR'; % 'RPIReID' 'Cleveland'; 'VIPeR'; 'iLIDS'; % set the data set name
switch dataset_name
    case {'VIPeR'}
        load ([dropbox_folder '/dataset/VIPeR_Images.mat'], 'I','gID', 'img_path', 'camID');
    case {'iLIDS'}
        load ([dropbox_folder '/dataset/iLIDS_Images.mat'], 'I','gID', 'img_path', 'camID');
    case {'Cleveland'}
        load ([dropbox_folder '/dataset/Cleveland_Image.mat'], 'I','gID', 'img_path', 'camID');
    case {'RPIReID'} 
        load ([dropbox_folder '/dataset/RPIReID_Images.mat'], 'I','gID', 'img_path', 'camID');
    case {'CAVIAR'} 
        load ([dropbox_folder '/dataset/CAVIAR_Images.mat'], 'I','gID', 'img_path', 'camID');    
    case {'KTHZ'}
end
% setting the preprocessing filter
prefilter_gaussian = []; %fspecial('gaussian', 9, 1.5);

imsz = [128 48]; % set the image size.
step =flipud([4, 4; 8 8; 16 16;  21 48]); % set the moving step size of the region.
BBoxsz =flipud([8,8; 16 16; 32 32; 22 48]); % set the region size.
numP = [6 14 75 341];

step = step(numP == num_patch,:);
BBoxsz = BBoxsz(numP == num_patch,:);

%% extract features with different region division set up.
for kk =1: size(step,1)
    kk
    % Generate the grid boundingbox with the setup of region division setting
    [region_idx, BBox] =GenerateGridBBox(imsz, BBoxsz(kk, :), step(kk, :));
    DataSet.region = region_idx;
    DataSet.BBox =BBox;
    DataSet.num_HistBin =16; % set the number of bins for color histogram
    for i =1:length(I)
        i
        % if the preprocessing filter is set, then extract features from
        % the filtered images
        if ~isempty(prefilter_gaussian) 
            I_filtered{i} = single(imfilter(I{i}, h));
        else
            I_filtered{i} = single(I{i});
        end
        if size(I_filtered{i},1) ~= imsz(1) || size(I_filtered{i},2) ~= imsz(2) 
            [F(i,:), feature_list] = ComputePatchHistogramMoment(imresize(I_filtered{i},imsz),region_idx,BBox, DataSet.num_HistBin );
        else
            [F(i,:), feature_list] = ComputePatchHistogramMoment(I_filtered{i},region_idx,BBox, DataSet.num_HistBin );
        end
    end
    % save the feature and useful information in the Dataset struct.
    Featurename = [dataset_name '_HistMoment' num2str(length(region_idx)) 'Patch'];
    DataSet.data = single(F);
    DataSet.name = Featurename;
    DataSet.idxfeature = feature_list;
    DataSet.prefilter_gaussian = prefilter_gaussian;
    DataSet.filepath = img_path;
    DataSet.ID = gID;
    DataSet.camID = camID;
    if ~isempty(prefilter_gaussian)
        save(['Feature/' Featurename '_wPreFiltering.mat'], 'DataSet','-v7.3')
    else
        save(['Feature/' Featurename '_woPreFiltering.mat'], 'DataSet','-v7.3')
    end
    clear F DataSet Featurename feature_list
end