clc
clear

%% load data
interactionA = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\06-miRNA-disease');

interactionB = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\08-lncRNA-miRNA');
interactionB = interactionB';

interactionC = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\05-lncRNA-disease');

sl1 = importdata('C:\Users\ahmet\Desktop\NFMCLDA-main\data\miRNASimilarity.txt');
sd1 = readmatrix('C:\Users\ahmet\Desktop\NFMCLDA-main\data\07-disease-disease');

%%
x = gKernel(495,interactionA);
y = gKernel(412,interactionA');

sl = (x+sl1)*0.5;
L = transition_matrix_network_propagation1(sl);

sd = (y+sd1)*0.5;
D = transition_matrix_network_propagation1(sd);


sm = gKernel(240,interactionC);
M = transition_matrix_network_propagation1(sm);


A_ori = interactionA;

score = thrrw(L, D, interactionA, interactionB, interactionC, M, 1, 3, 1, 0.3, 0.1, 0.9);
T = score;
trIndex = double(T ~= 0);
[WW,iter] = MC(0.9, 0.6, 100, T, trIndex, 0.001, 300, 0, 1);

%% ====================== LOOCV AYARI ======================
[row, col] = size(A_ori);

% Değeri 1 olan tüm etkileşimlerin indeksleri
[row_index, col_index] = find(A_ori == 1);
link_num = sum(A_ori(:));      % toplam etkileşim sayısı

CV = link_num;                 % LOOCV: her fold için 1 etkileşim
random_index = randperm(link_num);

size_of_CV = 1;                % her fold'ta 1 test örneği

%% ====================== LOOCV DÖNGÜSÜ ======================
result = zeros(1,7);

for k = 1:CV
    fprintf('LOOCV: round = %d/%d\n', k, CV);

    % Test için seçilecek tek etkileşimin index’i
    if (k ~= CV)
        test_row_index = row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index = col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index = row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index = col_index(random_index((size_of_CV*(k-1)+1):end));
    end

    % Eğitim matrisi: test etkileşimini 0 yap
    train_set = A_ori;
    test_link_num = size(test_row_index, 1);
    for i = 1:test_link_num
        train_set(test_row_index(i), test_col_index(i)) = 0;
    end

    % Modeli yeniden hesapla
    score = thrrw(L, D, train_set, interactionB, interactionC, M, 1, 3, 1, 0.3, 0.1, 0.9);
    T = score;
    trIndex = double(T ~= 0);
    [WW, iter] = MC(0.9, 0.6, 100, T, trIndex, 0.001, 300, 0, 1);

    % Performans metriklerini topla
    result = result + model_evaluate(A_ori, WW, train_set);
end

disp('aupr,auc,sen,spec,precision,accuracy,f1')
result = result / CV;   % LOOCV ortalaması

save result_loocv.mat