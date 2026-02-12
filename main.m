clc
clear

%% load data
interactionA = xlsread('C:\Users\ahmet\Desktop\NFMCLDA-main\data\05-lncRNA-disease');
interactionB = xlsread('C:\Users\ahmet\Desktop\NFMCLDA-main\data\08-lncRNA-miRNA');
interactionC = xlsread('C:\Users\ahmet\Desktop\NFMCLDA-main\data\06-miRNA-disease');

sl1 = xlsread('C:\Users\ahmet\Desktop\NFMCLDA-main\data\04-lncRNA-lncRNA');
sd1 = xlsread('C:\Users\ahmet\Desktop\NFMCLDA-main\data\07-disease-disease');

%% gkernel + transition matrix
x = gKernel(240,interactionA);
sl = (x+sl1)*0.5;
L = transition_matrix_network_propagation1(sl);

y = gKernel(412,interactionA');
sd = (y+sd1)*0.5;
D = transition_matrix_network_propagation1(sd);

sm = gKernel(495,interactionC);
M = transition_matrix_network_propagation1(sm);

%% thrrw
A_ori = interactionA;

score = thrrw(L,D,interactionA,interactionB,interactionC,M, 1,3, 1,0.3,0.1,0.9);

%% MC
T = score;

trIndex = double(T ~= 0);

[WW,iter] = MC(0.9,0.6, 100, T, trIndex, 0.001, 300, 0,1);