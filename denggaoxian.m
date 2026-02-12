%Draw contour lines for the correlation matrix
LD=xlsread('...\data\05-lncRNA-disease');
contour(LD);
xlabel('diseases');
ylabel('lncRNAs');