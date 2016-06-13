clc
clear 
%% appearance feature extraction in X. Fei, et al ECCV 2014
load 'dataset\cuhk03_release\cuhk_03_detected.mat'
name = 'cuhk03_detected_Hist75Patch.mat';

imsz = [128 48];
step = [8,8]; %[4, 4; 8 8; 16 16;  21 48]
BBoxsz = [16,16]; % [8,8; 16 16; 32 32; 22 48]
n8LBP_Mapping = getmapping(8,'u2');
n16LBP_Mapping = getmapping(16,'u2');
num_colorChn = 8;
num_bin = 16;

[region_idx, BBox, region_mask] =GenerateGridBBox(imsz, BBoxsz, step);
tmpF = zeros(numel(I),numel(region_idx)*(num_bin*num_colorChn+n8LBP_Mapping.num+n16LBP_Mapping.num),'single');
dim_color = num_bin*num_colorChn;
dim_LBP = n8LBP_Mapping.num+n16LBP_Mapping.num;

bin = repmat(1e-5:1/num_bin:1,5,1);
bin = [bin; [16:((235-15)/num_bin):235] ./255];
bin = [bin; [16:((240-15)/num_bin):240] ./255];
bin = [bin; [16:((240-15)/num_bin):240] ./255];

for i = 1:numel(I)
    i
    tmpSeq = I{i};
    if ~iscell(tmpSeq)
        tmpSeq = {tmpSeq};
    end
    tmpF_color = zeros(numel(tmpSeq),num_bin*num_colorChn,numel(region_idx));
    tmpF_lbp = zeros(numel(tmpSeq),n8LBP_Mapping.num+n16LBP_Mapping.num,numel(region_idx));
    for f = 1:numel(tmpSeq)
        tmpSeq{f} = imresize(tmpSeq{f},imsz);
        tmpRGB = im2double(tmpSeq{f});
        tmpHSV = rgb2hsv(tmpRGB);
        tmpYUV = rgb2ycbcr(tmpRGB);
        tmpGray = rgb2gray(tmpRGB);
        tmpChn = cat(3,tmpRGB, tmpHSV(:,:,1:2), tmpYUV);
        for c = 1:size(tmpChn,3)
            tmp_ch = tmpChn(:,:,c);
            for bb = 1:size(region_mask,2)
                tmpHist = hist(tmp_ch(region_idx{bb}),bin(c,:));
                tmpF_color(f,(c-1)*num_bin+1:c*num_bin,bb) = tmpHist./sum(tmpHist);
            end
        end
        
        for bb = 1:numel(region_idx)
            % LBP--n8u2r1
            tmpF_lbp(f,1:n8LBP_Mapping.num,bb) = lbp(tmpGray(BBox(bb,2):BBox(bb,4), ...
                BBox(bb,1):BBox(bb,3),:),1,n8LBP_Mapping.samples,n8LBP_Mapping,'nh')';
            % LBP--n16u2r2
            tmpF_lbp(f,n8LBP_Mapping.num+1:end,bb) = lbp(tmpGray(BBox(bb,2):BBox(bb,4), ...
                BBox(bb,1):BBox(bb,3),:),2,n16LBP_Mapping.samples,n16LBP_Mapping,'nh')';
        end        
    end
    tmpF_color = reshape(tmpF_color,size(tmpF_color,1),[]);
%     tmp_feat = normc_safe(tmpF_color');
%     tmpF_color = tmp_feat';
    tmpF_lbp = reshape(tmpF_lbp,size(tmpF_lbp,1),[]);
%     tmp_feat = normc_safe(tmpF_lbp');
%     tmpF_lbp = tmp_feat';
    tmpF(i,:) = [mean(tmpF_color,1) mean(tmpF_lbp,1)];
end
FeatureAppearance = tmpF;
save(name,'FeatureAppearance','-v7.3');