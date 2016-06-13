% calculate the matching characteristics for original LFDA algorithm.
% the major difference between this and compute_rank2.m is no kernel
% technique is used. So the projection matrix is d'-by-d
% By Fei Xiong, 
%    ECE Dept, 
%    Northeastern University 
%    2013-11-04
% Input: 
%       Method: the distance learning algorithm struct. In this function
%       two field are used. 
%               P is the projection matrix. d'-by-d
%               kernel is the name of the kernel function. 
%       train: The data used to learn the projection matric. Each row is a
%               sample vector. Ntr-by-d (not used in this function)
%       test: The data used to test and calculate the CMC for the
%               algorithm. Each row is a sample vector. Nts-by-d
%       ix_partition: the randomly generated partition of the test set.
%               Each row is a randomly generated partition. 1 represents
%               this test sample is used as reference sample, while 0
%               represents such sample is used as probe sample. Nit-by-Nts
%       IDs: The identity of the samples in the test set. Nts-by-1, where
%               Nts is the size of test set. Nts-by-1

function [R ,Alldist , ixx] = compute_rank_oLFDA(Method, train, test, ix_partition, IDs)

A = Method.P;

for k =1:size(ix_partition,1) % calculate the CMC for each random partition.
    % set the reference and prob data set.
    ix_ref = ix_partition(k,:) ==1;
    ix_prob = ix_partition(k,:) ==0;
    X_ref = test(ix_ref, : );
    X_prob = test(ix_prob, : );
    ref_ID = IDs(ix_ref);
    prob_ID = IDs(ix_prob);
    % calculate the distance and ranking for each prob sample
    for i =1: size(X_prob,1)
        diff = bsxfun(@minus, X_ref,X_prob(i, :));
        diff = A*diff';
        dis(i, :) = sum(diff.^2,1);
        [tmp, ix] = sort(dis(i, :));
        r(i) =  find(ref_ID(ix) == prob_ID(i));
        ixx(i,:)=ix;
    end
    R(k, :) = r;
    Alldist{k} = dis;
end
return;