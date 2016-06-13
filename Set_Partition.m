function Set_Partition(dataset, num_gal, num_test)


% load partition
% num_gal = 36;

% Default setting for each data set 
if nargin == 1    
    num_gal = 1;
    switch dataset
        case {'VIPeR'}
            num_test = 316;
        case {'iLIDS'}
            num_test = 60;
        case {'CAVIAR'}
            num_test = 36;
        case {'3DPeS'}
            num_test = 95;
        otherwise
            warning('Needs more parameters for unknown dataset');
    end
end
if num_gal > 1
    partition_name = 'multishot';
else
    partition_name = 'Random';
end


% num_gal = 3;
% if strcmp(partition_name,'Random') % make sure the single shot has correct set up
%     num_gal = 1;
% end

% generate 10 partitions with random selection of training identities 10
% trails for each partition with random selection of gallary image

num_partition =10;
num_trail = 10;

load(['dataset/' dataset '_Images.mat']);

ID = unique(gID);
ID = sort(ID);
for idx_partition=1:num_partition
    idxtemp = randperm(length(ID));
    idx_test = idxtemp(1:num_test);
    idx_train = idxtemp(num_test+1:end);
    gIDtemp = gID;
    while 1
        [C,ia,ib] = intersect(gIDtemp,ID(idx_train));
        if isempty(ia)
            break;
        end
        gIDtemp(ia) = -1;
    end
    Partition(idx_partition).idx_train =uint16(find(gIDtemp==-1));
    Partition(idx_partition).idx_test =uint16(find(gIDtemp>0));
    Partition(idx_partition).num_trainPerson = length(idx_train);
end


% generate gallery part and prob part for both training and testing data.
for idx_partition=1:num_partition
    ID_train = gID(Partition(idx_partition).idx_train);
    uID_train = unique(ID_train);
    idx_train_gallery = zeros(num_trail, length(Partition(idx_partition).idx_train));
    for m =1: length(unique(ID_train)) % for each person random choosing gallery sample
        iix_tID = find(ID_train == uID_train(m) );
        iix_temp = zeros(num_trail,num_gal);
        for t = 1:num_trail
            iix_temp(t,:) = randperm(length(iix_tID),num_gal);
        end
        for n = 1:num_gal
            temp = sub2ind(size(idx_train_gallery), [1:num_trail], iix_tID(iix_temp(:,n)));
            idx_train_gallery(temp)= 1; % gallery index in the training set for the individual of m.
        end
    end    
    
    ID_test = gID(Partition(idx_partition).idx_test);
    uID_test = unique(ID_test);
    idx_test_gallery = zeros(num_trail, length(Partition(idx_partition).idx_test));
    for m =1: length(unique(ID_test)) % for each person random choosing gallery sample
        iix_tID = find(ID_test == uID_test(m) );
        iix_temp = zeros(num_trail,num_gal);
        for t = 1:num_trail
            iix_temp(t,:) = randperm(length(iix_tID),num_gal);
        end
        for n = 1:num_gal
            temp = sub2ind(size(idx_test_gallery), [1:num_trail], iix_tID(iix_temp(:,n)));
            idx_test_gallery(temp)= 1; % gallery index in the training set for the individual of m.
        end        
    end
    Partition(idx_partition).ix_train_gallery = idx_train_gallery>0;
    Partition(idx_partition).ix_test_gallery = idx_test_gallery>0;
    
    [ix_pos_pair, ix_neg_pair]=GeneratePair(ID_train);
    Partition(idx_partition).idx_train_pos_pair = uint16(ix_pos_pair); % positive pairs for train set
    Partition(idx_partition).idx_train_neg_pair = uint16(ix_neg_pair); % negative pairs for train set
%     for k =1: num_trail
%         % random permutated negative pairs for train set
%         ix_radp(k,:) = randperm(length(ix_neg_pair));
%     end
%     Partition(idx_partition).idx_train_radp = uint32(ix_radp);
    [ix_pos_pair, ix_neg_pair]=GeneratePair(ID_test);
    Partition(idx_partition).idx_test_pos_pair = uint16(ix_pos_pair); % positive pairs for test set
    Partition(idx_partition).idx_test_neg_pair = uint16(ix_neg_pair); % negative pairs for test set
%     for k =1: num_trail
%         % random permutated negative pairs for test set
%         ix_radp(k,:) = randperm(length(ix_neg_pair));
%     end
%     Partition(idx_partition).idx_train_radp = uint32(ix_radp);
end
if num_gal > 1 % keep N for multishot scenario
    partition_name = [partition_name, '_' num2str(num_gal)];
end
save(['feature/' dataset '_Partition_' partition_name '.mat'], 'Partition', '-v7.3');

return