% calculate RGB, LAB, LOG-RGB and Gabor feature mapping of the input image.
% By Fei Xiong, 
%    ECE Dept, 
%    Northeastern University 
%    2013-11-04
% Input: 
%       I is the inuput RGB image. H*W*3
%       option is the struct contains the flag indicating which feature
%           channel is required to be computed.
% Output:
%       F: the feature maps of the image I. H*W*Nf, Nf is the number of
%           feature channels.
%       feature_list: the name of feature channels.
% ATTENTION: REQUIRE COLORSPACE LIBARY
function [F, feature_list] = ComputeFeature(I, option)

ix_feat =0;
if option.RGB
    F(:,:,ix_feat+[1:3]) = single(I);
    feature_list(ix_feat+[1:3]) ={'RGB-R', 'RGB-G', 'RGB-B'} ;
    ix_feat = ix_feat+3;
end
if option.LAB
    F(:,:,ix_feat+[1:3])  = single(colorspace('lab<-rgb', I));
    feature_list(ix_feat+[1:3])={'LAB-L', 'LAB-A', 'LAB-B'} ;
    ix_feat = ix_feat+3;
end
if option.LOG
    F(:,:,ix_feat+[1:3])  = single(log2(single(I)+1));
    feature_list(ix_feat+[1:3])={'LOG-L', 'LOG-A', 'LOG-B'} ;
    ix_feat = ix_feat+3;
end

Igray = rgb2gray(uint8(I));
if ~isempty(option.GaborFilter);
    for i =1:length(option.GaborFilter)
        F(:,:,ix_feat+1) = ...
        sqrt(single(imfilter(Igray,squeeze(option.GaborFilter{i}(1,:,:)), 'symmetric').^2 +...
        imfilter(Igray,squeeze(option.GaborFilter{i}(2,:,:)), 'symmetric').^2));
        feature_list(ix_feat+1) ={'Gabor'};
        ix_feat = ix_feat +1;
%         subplot(3, 4, i), imshow(squeeze(F(:,:,ix_feat)) , []);
    end
end

end