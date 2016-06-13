% test script 
% Set up computer environment
t1 = clock;
timstr =[ num2str(int32(t1(1))) '_'  num2str(int32(t1(2))) '_'  num2str(int32(t1(3)))...
    '_'  num2str(int32(t1(4))) '_'  num2str(int32(t1(5))) '_'  num2str(int32(t1(6)))];
%% set up Experiment parameters
tic
Set_Exp_Parameter
disp('Data loading time is:')
eltime_dataLoading = toc
%%
rng('default'); % Fix the random number generation for repeat results
for k =1:length(kname) % kernel for loop, select different kernel
    AlgoOption.kernel =kname{k};    
    savename = ['Result_' dataset_name '_' AlgoOption.name '_' fname '_' partition_name '_' kname{k} '_' timstr '.mat'];
    display(savename);
    if ~iscell(Feature)
        tmpF{1} = Feature;
        Feature = tmpF;
    end
    
    for c = 1:numel(Feature)
        if numel(gID)~=size(Feature{c},1) % make sure the feature is N-by-d
            feature_set{c} = double(Feature{c})';
        else
            feature_set{c} = double(Feature{c});
        end

        if ~isempty(strfind(AlgoOption.func, 'PCCA')) && strcmp(AlgoOption.kernel, 'linear') 
            % linear kernel for PCCA using sqrt
            feature_set{c} =sqrt(feature_set{c});
        end
        if AlgoOption.doPCA == 1
            if ~strcmp(AlgoOption.func,'KISSME')
%                 For svmml, L2 normalization is applied, for KISSME do
%                 nothing
                feature_set{c} = normc_safe(feature_set{c}'); % L2 normalizaiton
                feature_set{c} = feature_set{c}';
            end
        else
            % Others apply L1 normalization                
            feature_set{c} = bsxfun(@times, feature_set{c}, 1./sum(feature_set{c},2));
            if ~strcmp(AlgoOption.kernel, 'chi2-rbf')
                feature_set{c} =feature_set{c}*100;
            end

        end
    end
    %%
    for idx_partition=1:length(Partition) % partition loop        
        display(['==============' kname{k} ' ' num2str(idx_partition),' set ============='])
        for iset = 1:2
            %% iset: 1. validation step, 2. testing step
            if iset == 1
                %% load validating data
                % inverse the training and testing set to validate the prior                     
                idx_train = Partition(idx_partition).idx_test ;
                idx_test = Partition(idx_partition).idx_train ;
                ix_train_neg_pair = Partition(idx_partition).idx_test_neg_pair; % negative pair index
                ix_train_pos_pair = Partition(idx_partition).idx_test_pos_pair; % positive pair index
                ix_test_gallery =Partition(idx_partition).ix_train_gallery;
            else
                %% papre training and testing data
                idx_train = Partition(idx_partition).idx_train ;
                idx_test = Partition(idx_partition).idx_test ;
                ix_train_neg_pair = Partition(idx_partition).idx_train_neg_pair;
                ix_train_pos_pair = Partition(idx_partition).idx_train_pos_pair;
                ix_test_gallery =Partition(idx_partition).ix_test_gallery;
            end
            
            
            for c = 1:numel(feature_set)                    
                train{c} = feature_set{c}(idx_train,:); % training set
                test{c} = single(feature_set{c}(idx_test,:)); % test set                    
            end
            
            %% training step
            if ~strcmp(AlgoOption.name, 'LFDA') % setting the training pair samples
                Nneg = min(AlgoOption.npratio* length(ix_train_pos_pair), length(ix_train_neg_pair));
                ix_pair = [ix_train_pos_pair ; ix_train_neg_pair(1:Nneg,:) ]; % both positive and negative pair index
                y = [ones(size(ix_train_pos_pair,1), 1); -ones(Nneg,1)]; % annotation of positive and negative pair
            end
            if strcmp(AlgoOption.name, 'rPCCA')
                AlgoOption.w = ones(1, size(train{1},2));
                switch AlgoOption.kernel
                    case {'linear'}
                        AlgoOption.lambda =0.01; %0.3;
                    case {'chi2'}
                        AlgoOption.lambda =0.01; %0.3;
                    case {'chi2-rbf'}
                        AlgoOption.lambda =0.01;
                end
            end
            tic
            for c = 1:numel(train)
                switch AlgoOption.func
                    case {'MFA'}
                        [algo{c}, V] = MFA(single(train{c}), gID(idx_train)', AlgoOption);
                    case {'PCCA'}
                        [algo{c}, L, AKA] = PCCA(single(train{c}), ix_pair, y, AlgoOption);
                    case {'LFDA'}
                        [algo{c}, V]= LFDA(single(train{c}),gID(idx_train)' ,AlgoOption);
                    case {'oLFDA'}
                        [algo{c}, V]= oLFDA(single(train{c}),gID(idx_train)' ,AlgoOption);
                    case {'svmml'}
                        [algo{c}] = svmml_learn_full_final(single(train{c}),gID(idx_train)',AlgoOption);
                    case {'KISSME'}                        
                        [algo{c}] = kissme(train{c}',ix_pair,y,AlgoOption);
                end
            end
            disp('The training time is:')
            eltime_training = toc
            %% testing step
            tic
            switch algoname
                case {'oLFDA'}
                    [r, dis] = compute_rank_oLFDA(algo, train, test, ix_test_gallery, gID(idx_test));
                case {'svmml'}
                    [r, dis] = compute_rank_svmml(algo,train,test,ix_test_gallery, gID(idx_test));
                case {'KISSME'}
                    [r, dis] = compute_rank_KISSME(algo,test,ix_test_gallery,gID(idx_test));
                otherwise
                    [r, dis] = compute_rank2(algo, train, test, ix_test_gallery, gID(idx_test));
            end
            disp('The testing time is:')
            eltime_testing = toc
            [a, b] = hist(r',1:316);
            % when the probe set is not the same as test gallery set, it will be
            % labeled as "-1"
            if min(min(double(ix_test_gallery)))<0
                a = cumsum(a)./repmat(sum(ix_test_gallery==-1,2)', size(a,1),1);
            else
                a = cumsum(a)./repmat(sum(ix_test_gallery==0,2)', size(a,1),1);
            end
            a = a';
            for itr =1: size(a,1)
                rr(itr,:)= [a(itr,1) a(itr,5) a(itr,10) a(itr,20) a(itr,25) a(itr,50)];
                display(['itration ' num2str(itr) ' Rank 1 5 10 20 25 50 accuracy ===>' num2str(rr(itr,:)) ' ====']);
            end
            display(num2str([mean(rr,1)]));
            
            elTime.dataLoading = eltime_dataLoading;
            elTime.testing = eltime_testing;
            elTime.training = eltime_training;
            Method{idx_partition, iset} = algo{1};
            if strcmp(algoname, 'svmml')
                for c = 1:numel(algo)
                    tmpA{c} = algo{c}.A;
                    tmpB{c} = algo{c}.B;
                    tmpb{c} = algo{c}.b;
                end
                Method{idx_partition, iset}.A = tmpA;
                Method{idx_partition, iset}.B = tmpB;
                Method{idx_partition, iset}.b = tmpb;
            else
                for c = 1:numel(algo)
                    tmpP = algo{c}.P;
                end
                Method{idx_partition, iset}.P = tmpP;
            end
            Method{idx_partition, iset} = algo{1};
            Method{idx_partition, iset}.Prob = a;
            Method{idx_partition, iset}.Ranking = r;
            Method{idx_partition, iset}.Dist = dis;
            Method{idx_partition, iset}.time = elTime;
            
            save(savename, 'AlgoOption','Method');
        end % validation/test loop
    end % partition loop    
    clear Method
end % kernel loop
