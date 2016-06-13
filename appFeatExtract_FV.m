function [ F ] = appFeatExtract_FV(I, train, camID, numPatch)
%APPFEATEXTRACT_FV Summary of this function goes here
%   Detailed explanation goes here

% addpath(genpath('~/Dropbox/Research/Code/piotr_toolbox/'));
% run('~/Dropbox/Research/Code/VLFeat/vlfeat-0.9.17/toolbox/vl_setup.m');
step =flipud([8 8; 16 16; 21 48; 64 24; 128 48]); % set the moving step size of the region.
BBoxsz =flipud([16 16; 32 32; 22 48; 64 24; 128 48]); % set the region size.
numP = [1 4 6 14 75];
imsz = [128 64];
numbin = 16;
numChn = 17;%14+2;
kk = find(numP == numPatch);

pPyramid = chnsPyramid();
pPyramid.pChns.shrink = 1;
pPyramid.nPerOct = 2;
pPyramid.nApprox = 0;
pPyramid.pChns.pColor.colorSpace = 'rgb';
pPyramid.pChns.pGradMag.enabled = 0;


[~, BBox, region_mask] = GenerateGridBBox(imsz, BBoxsz(kk,:), step(kk,:));

chnFeat = cell(1,numel(I));
fprintf('Begin to extract channel feature... \n');tic
for i = 1:numel(I) % per person      
    if mod(i,round(numel(I)/10))==0
            fprintf('.');
    end
    tmpI = I{i};       
    if iscell(tmpI)
        num_frame = numel(tmpI);
    else 
        num_frame = 1;
        tmpI = {tmpI};
    end
    
    if 0 %1--HOG FV;   0---LDFV
        tmpChnFeat = zeros(numChn,imsz(1)*imsz(2),num_frame, 'single');
        for f = 1:num_frame
            tmp_im = imresize(tmpI{f},imsz);
            tmp_data = {};
            % get the channel feature         
            pyramid = chnsPyramid( tmp_im, pPyramid );
            tmp_data_ch = pyramid.data{1};
            tmp_data{1} = tmp_data_ch(:,:,1:3);                 % RGB
            tmp_data{2} = rgbConvert(tmp_data{1},'luv',1);      % LUV
            tmp_data{3} = rgbConvert(tmp_data{1},'hsv',1);      % HSV
            tmp_data{4} = tmp_data_ch(:,:,4:9);                 % Gradient
            tmp_data = cat(3,tmp_data{:});
            tmp_data(:,:,6) = []; % remove duplicated V channel
            tmp_data(:,:,15) = repmat([1:imsz(1)]',1,imsz(2));   % x,y location
            tmp_data(:,:,15) = tmp_data(:,:,15)./imsz(1);
            tmp_data(:,:,16) = repmat([1:imsz(2)],imsz(1),1);
            tmp_data(:,:,16) = tmp_data(:,:,16)./imsz(2);
            
            tmp_data = permute(tmp_data, [3 1 2]);            
            tmp_data = reshape(tmp_data,size(tmp_data,1), []);
            tmpChnFeat(:,:,f) = tmp_data;
        end
    else 
        tmpChnFeat = zeros(numChn,imsz(1)*imsz(2),num_frame, 'single');
        for f = 1:num_frame
            tmp_im = imresize(tmpI{f},imsz);
            tmp_data = zeros(imsz(1),imsz(2),17);
            % get the channel feature
            tmp_data(:,:,1) = repmat([1:imsz(1)]',1,imsz(2));
            tmp_data(:,:,1) = tmp_data(:,:,1)./imsz(1);
            tmp_data(:,:,2) = repmat([1:imsz(2)],imsz(1),1);
            tmp_data(:,:,2) = tmp_data(:,:,2)./imsz(2);
            tmp_hsv = rgb2hsv(tmp_im);
            tmp_data(:,:,3:5) = tmp_hsv;
            for c = 1:3
                [tmpX,tmpY] = imgradientxy(reshape(tmp_hsv(:,:,c),imsz(1),imsz(2)));
                [tmpXX,~] = imgradientxy(tmpX);
                [~,tmpYY] = imgradientxy(tmpY);
                tmp_data(:,:,(c-1)*4+6:c*4+5) = cat(3,tmpX,tmpY,tmpXX,tmpYY);
            end
            tmp_data = permute(tmp_data, [3 1 2]);
            tmp_data = reshape(tmp_data,size(tmp_data,1),[]);
            tmpChnFeat(:,:,f) = tmp_data;
        end
    end
    chnFeat{i} = tmpChnFeat;
end
fprintf('Done!\n');toc
fprintf('Begin to build GMM model...\n');tic
% GMM encoding
idx = 1:numel(I);
idx_train = ismember(idx,train);
for s = 1:numPatch    
    tmpChnFeat = cellfun(@(x) x(:,logical(region_mask(:,s)),:), chnFeat(idx_train),'UniformOutput',0);
    tmpChnFeat = cellfun(@(x) reshape(x,numChn,[]), tmpChnFeat,'UniformOutput',0);
    tmpChnFeat = cell2mat(tmpChnFeat);
    % subsample
    MAXsample = 100000;
    idxsub = randsample(size(tmpChnFeat,2),min(MAXsample,size(tmpChnFeat,2)));  
    [means{s}, covariances{s}, priors{s}] = vl_gmm(tmpChnFeat(:,idxsub),numbin,'NumRepetitions',10);    
end
fprintf('Done!\n');toc
fprintf('Begin to extract fisher vectors...\n');tic
% FV encode
for i = 1:numel(chnFeat)
    if mod(i,round(numel(I)/10))==0
            fprintf('.');
    end
    tmpChnFeat = chnFeat{i};
    tmpF_perP = zeros(numChn*2*numbin*numPatch,size(tmpChnFeat,3),'single');
    for f = 1:size(tmpChnFeat,3)
        tmpFrame = tmpChnFeat(:,:,f);
        tmpF_perf = zeros(numChn*2*numbin,numPatch,'single');
        for s = 1:numPatch        
            tmpFrame_s = tmpFrame(:,logical(region_mask(:,s)));
            tmpF_perf(:,s) = vl_fisher(tmpFrame_s,means{s}, covariances{s}, priors{s});
        end
        tmpF_perP(:,f) = tmpF_perf(:);
    end
    % naive mean along the temporal dimension
    tmpF_mean = mean(tmpF_perP,2);
    % normalize on each FV
    tmpF_mean = sign(tmpF_mean).*(sqrt(abs(tmpF_mean)));    % power
    FVdim = numChn*2*numbin;
    for s = 1:numPatch        
        tmpF_mean((s-1)*FVdim+1:s*FVdim) = normc_safe(tmpF_mean((s-1)*FVdim+1:s*FVdim)); % L2
    end
    F(i,:) = tmpF_mean';
end
fprintf('Done!\n');
% F = sign(F).*(sqrt(abs(F)));
% F = normc_safe(F');
% F = F';




