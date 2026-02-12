function [result_mi]=gKernel(nm,interaction)  %Gaussian similarity
for i=1:nm
        sm(i)=norm(interaction(i,:))^2;  %The function of finding norm
    end
    gamam=nm/sum(sm')*1;   %Kernel bandwidth
    for i=1:nm
        for j=1:nm
            pkm(i,j)=exp(-gamam*(norm(interaction(i,:)-interaction(j,:)))^2);   
        end
    end   
    
    result_mi=pkm;
end