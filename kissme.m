% Wrapping data and modified the structure to fullfill the kissme* algorithm.
%    *Proposed by Martin Koestinger (2011), please refer the following paper 
%    if you use this code please cite the following paper:
%    Kostinger, M., Hirzer, M.,Wohlhart, P., Roth, P.M., Bischof, H.: 
%    Large scale metric learning from equivalence constraints. In CVPR 2012
    
%%
function [Method] = kissme(X,ix_pair,y,option)

params.numCoeffs = option.PCAdim; %dimensionality reduction by PCA to 34 dimension
params.N = size(ix_pair,2); %number of image pairs, 316 to train 316 to test
params.numFolds = option.nFold; %number of random train/test splits
params.saveDir = '';
params.pmetric = 0;

cHandle =  LearnAlgoKISSME(params);

s = learnPairwise(cHandle,X,ix_pair(:,1),ix_pair(:,2),y>0);
ds.(cHandle.type) = s;

Method.name = 'KISSME';
Method.ds = ds;
Method.P = ds.kissme.M;
Method.kernel = option.kernel;
Method.Prob = [];
Method.Dataname = option.dataname;
Method.Ranking = [];
Method.Dist = [];
Method.Trainoption = option;
return