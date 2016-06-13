INFO:
This is the matlab code of the paper 
X. Fei, M. Gou, O. Camps and M. Sznaier: "Person Re-Identification using Kernel-based Metric Learning Methods". In ECCV 2014.

STRUCTURE:
-All the main script and function codes are located under root path;
-All parsed datesets are located under "dataset/";
-All additional libs and data parsing codes are located under "Assistant Code/".

HOW TO RUN THE DEMO:
- Data preparation 
    Parsing images to appropriate format by running "load_DATASET.m" in folder "Assistant Code";
    If you want to use your own data set, please refer any load function for the detail of data structure 
- To run the example, please start with "run_demo.m";
- To get the final average results and PUR value, please run "Script_demo_result.m"
- To use ensemble fashion, one needs to run all the algorithms with all different features first and then run "Script_demo_Ensemble_result.m"
- Most parameters are set in "Set_Exp_Parameter.m", please modify them based on particular experiment. 
- Right now, this version support seven different algorithms. Among them, LFDA, rPCCA and MFA are proposed by ourselves; oLFDA and PCCA are re-implemented based on others' work; KISSME and svmml are modified from authors' code. PLEASE cite the corresponding paper properly. 

THIRD PARTY CODE:
- KISSME    
- svmml    
    PLEASE NOTE THAT THIS PART CODE COULD ONLY BE USED FOR RESEARCH PURPOSE.
- Local Binary Pattern (Assistant Code/LBP) 
    http://www.scholarpedia.org/article/Local_Binary_Patterns

DATASET:
Parsed iLIDS dataset is included in this package for quick testing. Please refer the following paper for more details about this dataset.
Zheng, W.S., Gong, S., Xiang, T.: Associating groups of people. In: BMVC (2009)

LICENSE:
Copyright (c) 2013, Fei Xiong and Mengran Gou @ Robust Systems Lab of Northeastern University
All rights reserved.
This package (EXCEPT svmml part) is licensed under BSD 3-Clause Licence. 

CITATION:
If you use this code please cite the following paper:
X. Fei, M. Gou, O. Camps and M. Sznaier: "Person Re-Identification using Kernel-based Metric Learning Methods". In ECCV 2014.
If you use oLFDA, PCCA, KISSME or svmml in this package, please refer the original paper properly.
- oLFDA
    Pedagadi, S., Orwell, J., Velastin, S., Boghossian, B.: Local sher discriminant analysis for pedestrian re-identification. In CVPR 2013
- PCCA
    Mignon, A., Jurie, F.: Pcca: A new approach for distance learning from sparse pairwise constraints. In CVPR 2012
- KISSME
    Kostinger, M., Hirzer, M.,Wohlhart, P., Roth, P.M., Bischof, H.: Large scale metric learning from equivalence constraints. In CVPR 2012
- svmml
    Li, Z., Chang, S., Liang, F., Huang, T.S., Cao, L., Smith, J.R.: Learning locally-adaptive decision functions for person verication. In CVPR 2013


CHANGELOG:
v0.0: the original version. 
v0.1: Fix LFDA, Set_Exp_Parameter 
v0.3: Add two more algorithms: MFA and SVMML(UIUC). Change the structure for multi layer projection
v0.4 - 04/09/2014: Update oLFDA, LFDA, moment feature extraction. Beta version released.
v1.0 - 08/27/2014: Wrapped up and release the public version
v1.1 - 07/23/2015: Update LBP feature extraction for Matlab 8.0 (2012b) or higher; 
v1.2 - Move to github
