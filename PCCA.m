% PCCA distance learning algorithm proposed in "PCCA: A New Approach for
% Distance Learning from sparse pairwise constraints" cvpr 2012
% By Fei Xiong, 
%    ECE Dept, 
%    Northeastern University 
%    2013-11-04
% INPUT
%   X: N-by-d data matrix. Each row is a sample vector.
%   ix_pair: the index for pairwise constraints.
%   y: the annotation for pairwise constraints. {+1, -1}
%   option: algorithm options
%       beta: the parameter in generalized logistic loss function
%       d: the dimensionality of the projected feature
%       eps: the tolerence value.
% OUTPUT
%   Method: the structure contains the learned projection matrix and
%       algorithm parameters setup.
%   l_old: the objective value
%   AKA: the regularizer value
% UPDATE LOG:
%   05/26/2014 speeding up the initialize part
function [Method, l_old, AKA]= PCCA(X, ix_pair, y, option)
ix_pair= double(ix_pair);
if length(size(X)) >2 % input feature are covariance matrix
    [Method, l_old, AKA]= PCCA_Cov(X, ix_pair, y, option);
    return;
end
A =[];
l_old = [];
AKA = [];
display(['begin PCCA ' option.kernel]);
beta= option.beta;
d = option.d;
eps = option.epsilon;
eta = 0.1;
if option.lambda>0
    W= eye(size(X,1))*option.lambda;
end

% compute the kernel matrix
Method = struct('rbf_sigma',0);
[K, Method] = ComputeKernel(X, option.kernel, Method);
K= K*size(K,1)/trace(K); % scale the kernel matrix

% identity basis
I = sparse(eye(size(X,1)));
% K*J_n in PCCA paper equation(10)
KJ= sparse([], [], [], size(ix_pair, 1), size(X,1)^2, 2*size(X,1)*size(ix_pair, 1));
tic
chop_num = 2000;  % chop into samll piece for speed up (2000 needs 17GB for CAVIAR)
% compute K*J_n as in equation 10
for ns = 1:chop_num:size(ix_pair,1)%-mod(size(ix_pair,1),chop_num)
    chop_mat = zeros(chop_num,size(X,1)^2);
    n = 1;
    for i = ns:min(chop_num+ns-1,size(ix_pair,1))
%         chop_row = zeros(1,size(X,1)^2);
        ix_row1= sub2ind(size(K), 1:size(X,1), ones(1,size(X,1))*ix_pair(i,1));
%         chop_row(1,ix_row1) = K(:,ix_pair(ns+i-1,1))- K(:,ix_pair(ns+i-1,2));
        chop_mat(n,ix_row1) = K(:,ix_pair(i,1))- K(:,ix_pair(i,2));
        ix_row2= sub2ind(size(K), 1:size(X,1), ones(1,size(X,1))*ix_pair(i,2));
%         chop_row(1,ix_row2) = -chop_row(1,ix_row1);
        chop_mat(n,ix_row2) = -chop_mat(n,ix_row1);
%         chop_mat(i,:) = chop_row;        
        n = n+1;
    end
    if chop_num+ns-1 < size(ix_pair,1)
        KJ(ns:chop_num+ns-1,:) = sparse(chop_mat);
    else 
        KJ(ns:size(ix_pair,1),:) = sparse(chop_mat(1:mod(size(ix_pair,1),chop_num),:));
    end
end
% for i = size(ix_pair,1)-mod(size(ix_pair,1),chop_num):size(ix_pair,1)
%     ix_row1= sub2ind(size(K), 1:size(X,1), ones(1,size(X,1))*ix_pair(i,1));
%     KJ(i,ix_row1)= K(:,ix_pair(i,1))- K(:,ix_pair(i,2));
%     ix_row2= sub2ind(size(K), 1:size(X,1), ones(1,size(X,1))*ix_pair(i,2));
%     KJ(i,ix_row2)= -KJ(i,ix_row1);
% end
% parfor i =1:size(ix_pair,1)
%     KJ_row = zeros(1,size(X,1)^2);
% %     compute K*J_n as in equation 10
%     ix_row1= sub2ind(size(K), 1:size(X,1), ones(1,size(X,1))*ix_pair(i,1));
% %     KJ(i,ix_row1)= K(:,ix_pair(i,1))- K(:,ix_pair(i,2));
%     KJ_row(1,ix_row1) = K(:,ix_pair(i,1))- K(:,ix_pair(i,2));
%     ix_row2= sub2ind(size(K), 1:size(X,1), ones(1,size(X,1))*ix_pair(i,2));
% %     KJ(i,ix_row2)= -KJ(i,ix_row1);
%     KJ_row(1,ix_row2) = -KJ_row(1,ix_row1);
%     KJ(i,:) = sparse(KJ_row);
% end
toc
%%
%initialization
A =randn(d, size(X,1))/1e4;
temp = K(ix_pair(:,1),:) - K(ix_pair(:,2),:);
temp = A*temp';
D = sum(temp.^2, 1);
AKA =trace(A*K*A');
l_old = sum(logistic_loss(D,y,option)) + option.lambda* AKA; %
% gradient search
cnt =0;
while 1
    L = -y'.*(1 -D); % there is a sign mistake in equation (3) and equation for L_n^t
    L = 1./ (1+exp(-beta*L));
    L = double(y'.*L);
    if option.lambda ==0
        Y = reshape(L*KJ,[size(X,1) size(X,1)]);
    else
        Y = reshape(L*KJ,[size(X,1) size(X,1)]) + W;
    end
%     % optimization eta does not work
%     temp =reshape(sum(KJ),[size(X,1) size(X,1)]);
%     temp = temp*K*Y*A'*A;
%     eta = trace(temp)/(2*trace(temp*Y));
    
    A_new = A - 2*A *eta *Y;%*K
    temp = K(ix_pair(:,1),:) - K(ix_pair(:,2),:);
    temp = A_new*temp';
    D = sum(temp.^2, 1);
    AKA_new =trace(A_new*K*A_new');
    l_new = sum(logistic_loss(D,y,option))+ option.lambda * AKA_new; % + option.lambda* trace(A*K*A')
    % adjust learning rate
    if l_new >  l_old
        eta = eta*0.9;
        if eta <1e-50
            break;
        else
            continue;
        end
    else
        eta = eta*1.1;
    end
    if mod(cnt,100)==0
        display(num2str([ cnt l_new  AKA  norm(A-A_new,'fro')/norm(A, 'fro') eta ]))
%         plot(D); drawnow; pause(0.1)
    end
    
    if l_old - l_new < eps && norm(A-A_new, 'fro')/norm(A, 'fro')<eps
        break;
    end
    l_old = l_new;
    A = A_new;
    AKA = AKA_new;
    cnt =cnt+1;
end
display(num2str([ cnt l_old  AKA norm(A-A_new,'fro')/norm(A, 'fro') eta ]));
%% save the algorithm information and trained projection matrix.
if option.lambda>0
    Method.name = 'rPCCA';
else
    Method.name = 'PCCA';
end
Method.P=A;
Method.kernel=option.kernel;
Method.Prob = [];
Method.Dataname = option.dataname;
Method.Ranking = [];
Method.Dist = [];
Method.Trainoption=option;
return;

% Computing the generalized logistic loss value, the objective function 
% value eq(2). 
function L =logistic_loss(D,y,option)
beta = option.beta;
L =1/beta*log(1+ exp(beta*(y'.*(D-1))));
return;
