% I           --- 1XN cell variable of all person images
% img_path    --- 1XN cell variable of all image's path (can be empty)
% camID       --- 1XN double variable of the camera number of each person
% gID         --- 1XN double variable of the ground truth ID of each person

%%
% Load *CAVIAR4REID dataset
% *Cheng, D.S., Cristani, M., Stoppa, M., Bazzani, L., Murino, V.: Custom pictorial
% structures for re-identification. In: British Machine Vision Conference (BMVC).
% pp. 68.1{68.11 (2011)

img_name = ls('*.jpg');
datapath = pwd;

num_img = size(img_name,1);
I = cell(1,num_img);
img_path = cell(1,num_img);
gID = zeros(1,num_img);
camID = zeros(1,num_img);
for i = 1:size(img_name,1)
    i
    imfile = img_name(i,:);
    img_path{i} = fullfile(datapath, imfile);
    tmpImg = imread(img_path{i});
    I{i} = imresize(tmpImg, [128,48]);
    gID(i) = str2num(imfile(1:4));
end
save('CAVIAR_Images.mat', 'I', 'img_path', 'gID', 'camID');