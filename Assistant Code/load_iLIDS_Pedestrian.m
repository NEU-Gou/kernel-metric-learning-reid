% I           --- 1XN cell variable of all person images
% img_path    --- 1XN cell variable of all image's path (can be empty)
% camID       --- 1XN double variable of the camera number of each person
% gID         --- 1XN double variable of the ground truth ID of each person

%%
% load *iLIDS dataset
% *Zheng, W.S., Gong, S., Xiang, T.: Associating groups of people. In: BMVC (2009)

clear
dataset_path = pwd;

dirc = dir(dataset_path );
cnt =1;
for i =1: length(dirc)
    if ~dirc(i).isdir && strcmp(dirc(i).name(end-2:end), 'jpg')
        I{cnt} = imread([dataset_path  dirc(i).name]);
        img_path{cnt} = [dataset_path  dirc(i).name];
        gID(cnt) = str2num(dirc(i).name(1:4));
        camID(cnt) = str2num(dirc(i).name(5:8));
        cnt = cnt +1;
    end
end
save('iLIDS_Images.mat','I','img_path','camID','gID','-v7.3');
