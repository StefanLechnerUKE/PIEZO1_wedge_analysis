fns = fieldnames(IndivTrimersRAW);

WT_CTL = [13 22 77 92 94];
WT_yoda1 = [60 65 81 83 87];


A2094W_cytoD = [3 13 26 59 62];
A2094W_cytoD_Rot = [0 0 0 0 0];

A2094W_Yoda1 = [8 40 48 56 79];
A2094W_Yoda1_Rot = [0 -140 140 145 0];

VF_AA_CTL = [];
VF_AA_Yoda1 = [];

VF_AA = [2 19 23 35 95];
VF_AA_Rot =[0 135 0 0 0];

VF_AA_Yoda1 = [15 24 30 42 55];
VF_AA_Yoda1_Rot = [0 0 0 0 0];

All_rot_1 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -140 140 145 0 0 0 0 0 0 0 0 0 0 0];

for k = VF_AA_Yoda1
Fig4_examples.(strcat(fns{k},'_VF_AA_Yoda1_',string(k))) = IndivTrimersRAW.(fns{k})
end
