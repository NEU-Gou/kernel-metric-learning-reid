% given the ID set, find the positive and negative pair of samples.
% NOTE that the negative pair samples are random permutated, so that taking
% first N pairs is the same as randomly pick the negative pairs.
function [ix_pos_pair, ix_neg_pair]=GeneratePair(ID, varargin)
ID = double(ID);
R = repmat(ID.^2, length(ID), 1);
R = R + R' -2*ID'*ID;
idx_triu = find(triu(ones(size(R)) - eye(size(R)))>0);
idx_pos = idx_triu(R(idx_triu(:)) ==0);
idx_neg = idx_triu(R(idx_triu(:)) ~=0);
rndp = randperm(length(idx_neg));
idx_neg = idx_neg(rndp);
[ix_neg_pair(:,1), ix_neg_pair(:,2) ]= ind2sub(size(R), idx_neg);
[ix_pos_pair(:,1), ix_pos_pair(:,2) ]= ind2sub(size(R), idx_pos);
if ~isempty(varargin)
    npratio = varargin{1}; % negative to positive pair ratio
    num_neg_pair = min(npratio*size(ix_pos_pair,1), size(ix_neg_pair,1));
    ix_neg_pair =ix_neg_pair(1:num_neg_pair, :);
end
return;