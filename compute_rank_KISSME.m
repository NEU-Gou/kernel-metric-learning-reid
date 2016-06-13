function [R,Alldist,ixx] = compute_rank_KISSME(Method,test,ix_partition,IDs)

for k = 1:size(ix_partition,1)
    ix_ref = ix_partition(k,:) == 1;
    if min(min(double(ix_partition))) < 0
        ix_prob = ix_partition(k,:) ==-1; 
    else
        ix_prob = ix_partition(k,:) ==0;
    end
    ref_ID = IDs(ix_ref);
    prob_ID = IDs(ix_prob);
    for c = 1:numel(test)
        M = Method{c}.ds.kissme.M;
        feat_ref = test{c}(ix_ref,:);
        feat_prob = test{c}(ix_prob,:);
        dist = sqdist(feat_prob', feat_ref' ,M);
        for p = 1:size(feat_prob,1)            
            [tmp, ix] = sort(dist(p, :));
            r(p) =  find(ref_ID(ix) == prob_ID(p));
            ixx(p,:)=ix;
        end
    end
    R(k, :) = r;
    Alldist{k} = dist; % distance matrix
end
