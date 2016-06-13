% voting with multiple ranking algorithm

function [MC] = Vote_MultiRankingAlgo(Method, Partition, gID, flag)
num_itr = 10;%size(Partition(1).ix_test_gallery,1);
Method = shiftdim(Method, length(size(Method))-2); % shift the algorithm matrix so that the 2-by-partition-by....
Method = reshape(Method, size(Method,1), size(Method, 2), []); 
for ix_part =1: size(Method,1) % partition loop
    idx_test = Partition(ix_part).idx_test ;
    ix_test_gallery =Partition(ix_part).ix_test_gallery;
    for m =1: num_itr % iteration loop
        p = zeros(size(Method{ix_part,2, 1}.Dist{m}));
        IDs= gID(idx_test);
        ix_ref = ix_test_gallery(m,:) ==1;
        ix_prob = ix_test_gallery(m,:) ==0;
        ref_ID = IDs(ix_ref);
        prob_ID = IDs(ix_prob);
        for k =1: size(Method, 3)
            if isempty(Method{ix_part,2,k})
                continue;
            end
            prob = Method{ix_part,1, k}.Prob(1:num_itr,1:numel(ref_ID));
            prob =[prob(:,1) prob(:,2:end) - prob(:, 1:end-1)];
            prob = prob(:,1:numel(ref_ID));
            dis = Method{ix_part,2, k}.Dist{m}; % test-by-reference
            prob = repmat(mean(prob,1), sum(ix_prob),1 );
            [~, idx] = sort(dis');
            idx = idx'; % test-by-reference
            temp = repmat([1:size(idx,1)]', 1, size(idx,2));
            idxtmp= sub2ind(size(idx), temp(:), idx(:));
            if flag>0
                p(idxtmp)=p(idxtmp)+(prob(:));
            else
                p(idxtmp)=p(idxtmp)+log(prob(:));
            end
        end
        for i =1: size(p,1)
            [tmp, ix] = sort(p(i, :),'descend');
            r(i) =  find(ref_ID(ix) == prob_ID(i));
            ixx(i,:)=ix;
        end
        [a, b]=hist(r, [1:sum(ix_ref)]);
        MC(ix_part,m,1:sum(ix_ref)) = a/sum(a);
%         if sum(a)~=sum(ix_prob)
%             display('wrong')
%         end
        clear r;
    end
end
return;