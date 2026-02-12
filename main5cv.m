clc
clear

%% load data
interactionA = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\06-miRNA-disease');

interactionB = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\08-lncRNA-miRNA');
interactionB = interactionB';

interactionC = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\05-lncRNA-disease');

sl1 = importdata('C:\Users\ahmet\Desktop\NFMCLDA-main\data\miRNASimilarity.txt');
sd1 = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\07-disease-disease');
sm1 = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\04-lncRNA-lncRNA');

%%
x = gKernel(495,interactionA);
sl = (x+sl1)*0.5;
L = transition_matrix_network_propagation1(sl);

y = gKernel(412,interactionA');
sd = (y+sd1)*0.5;
D = transition_matrix_network_propagation1(sd);

z = gKernel(240,interactionC);
sm = (z+sm1)*0.5;
M = transition_matrix_network_propagation1(sm);

%%
A_ori = interactionA;

score = thrrw(L, D, interactionA, interactionB, interactionC, M, 1, 3, 1, 0.3, 0.1, 0.9);

T = score;
trIndex = double(T ~= 0);
[WW,iter] = MC(0.9, 0.6, 100, T, trIndex, 0.001, 300, 0, 1);

%%
CV = 5; 

[row,col] = size(A_ori);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row_index,col_index] = find(A_ori==1); %%find all the elments that its value is 1
link_num = sum(A_ori(:));  %% caculate the numbers of the interaction

random_index = randperm(link_num);
size_of_CV = round(link_num/CV);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
result = zeros(1,7);
for k = 1:CV
    fprintf('cross validation: round= %d/%d\n', k, CV);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (k~=CV) %% We allocate all the interaction elements into 5 parts
        test_row_index = row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index = col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index = row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index = col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    train_set = A_ori;
    test_link_num = size(test_row_index,1);
    for i = 1:test_link_num %% let interaction elements' value equal to 0 in the test matrix
        train_set(test_row_index(i),test_col_index(i)) = 0;  %
    end
    
    score = thrrw(L, D, train_set, interactionB, interactionC, M, 1, 3, 1, 0.3, 0.1, 0.9);
    T = score;
    trIndex = double(T ~= 0);
    [WW,iter] = MC(0.9, 0.6, 100, T, trIndex, 0.001, 300, 0, 1);
    
    result = result+model_evaluate(A_ori,WW,train_set);
end

disp('aupr,auc,sen,spec,precision,accuracy,f1')
result = result/CV

%%
hold on
title('5-fold CV','FontWeight','normal')