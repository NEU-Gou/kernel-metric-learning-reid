% I           --- 1XN cell variable of all person images
% img_path    --- 1XN cell variable of all image's path (can be empty)
% camID       --- 1XN double variable of the camera number of each person
% gID         --- 1XN double variable of the ground truth ID of each person

%%
% load *VIPeR dataset
% *Gray, D., Tao, H.: Viewpoint invariant pedestrian recognition with an ensemble of
% localized features. In: Computer Vision{ECCV 2008, pp. 262{275. Springer (2008)

dataset_path = pwd;

dirc = dir([dataset_path 'cam_a/']);
cnt =1;
for i =1: length(dirc)
    if ~dirc(i).isdir && strcmp(dirc(i).name(end-2:end), 'bmp')
        I{cnt} = imread([dataset_path 'cam_a/' dirc(i).name]);
        img_path{cnt} = [dataset_path 'cam_a/' dirc(i).name];
        gID(cnt) = str2num(dirc(i).name(1:3));
        camID(cnt) = 1;
        cnt = cnt +1;
    end
end
dirc = dir([dataset_path 'cam_b/']);
for i =1: length(dirc)
    if ~dirc(i).isdir && strcmp(dirc(i).name(end-2:end), 'bmp')
        I{cnt} = imread([dataset_path 'cam_b/' dirc(i).name]);
        img_path{cnt} = [dataset_path 'cam_b/' dirc(i).name];
        gID(cnt) = str2num(dirc(i).name(1:3));
        camID(cnt) = 2;
        cnt = cnt +1;
    end
end

save('VIPeR_Images.mat','I','img_path','camID','gID','-v7.3');
