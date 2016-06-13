% I           --- 1XN cell variable of all person images
% img_path    --- 1XN cell variable of all image's path (can be empty)
% camID       --- 1XN double variable of the camera number of each person
% gID         --- 1XN double variable of the ground truth ID of each person

%%
% load 3DPeS* dataset
% *Baltieri, D., Vezzani, R., Cucchiara, R.: 3dpes: 3d people dataset for surveillance
% and forensics. In: Proceedings of the 1st International ACM Workshop on Mul-
% timedia access to 3D Human Objects. pp. 59{64. Scottsdale, Arizona, USA (Nov
% 2011)
dataset_path = pwd;
fid = fopen('3dpes_train.al');
I = {};
count = 1;
gID = [];
viewAng  =[];
tmp = fgetl(fid);
impath = {};
while tmp ~= -1 % read until the end of the file 
    tmp = fgetl(fid);
    if strcmp(tmp, '<annotation>')
        tmp = fgetl(fid);
        impath{count} = tmp(15:end-15);        
        tmpIM = imread([dataset_path,'/', impath{count}]);
        [maxY,maxX,~] = size(tmpIM);
        [IDstr dummy] = strtok(impath{count}, '_');
        [dummy IDstr] = strtok(IDstr, '\');
        gID(count) = str2num(IDstr(2:end));
        tmp = fgetl(fid);
        x1 = str2num(tmp(regexp(tmp,'<x1>')+4:regexp(tmp,'</x1>')-1));
        x2 = str2num(tmp(regexp(tmp,'<x2>')+4:regexp(tmp,'</x2>')-1));
        y1 = str2num(tmp(regexp(tmp,'<y1>')+4:regexp(tmp,'</y1>')-1));
        y2 = str2num(tmp(regexp(tmp,'<y2>')+4:regexp(tmp,'</y2>')-1));
        
        x1 = max(1,x1);
        x2 = min(x2, maxX);        
        y1 = max(1,y1);        
        y2 = min(y2,maxY);
        I{count} = imresize(tmpIM(y1:y2,x1:x2,:),[128,48]); % normlized to 128*48
        viewAng(count) = str2num(tmp(regexp(tmp,'<id>')+4:regexp(tmp,'</id>')-1));
        count = count+1
    end
end
fclose(fid);

% remove the individuals with only one image
[numh,xh] = hist(gID,unique(gID));
iso_idv = xh(numh == 1);
for i = 1:numel(iso_idv)
    rm_idx = find(gID == iso_idv(i));
    gID(rm_idx) = [];
    I(rm_idx) = [];
    viewAng(rm_idx) = [];
    impath(rm_idx) = [];
end
save('3DPeS_Images.mat', 'I', 'gID', 'viewAng', 'impath');
