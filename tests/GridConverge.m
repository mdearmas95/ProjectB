clc;clear;close all
N=[10 20 40 60 80 100 120];

NormU=[174 250 413 583 756 930 1156];

n=length(N);
for i=2:n-1
    NormDiff(i) = NormU(i+1) - NormU(i);
end

plot(N(1:6),NormDiff,'*')
figure(1)
title('Grid Convergence')
xlabel('Nodes')
ylabel('Norm of difference in U')
for i=1:6
Error(i)=abs(NormDiff(i)-NormDiff(5))/(NormDiff(5));
end

figure(2)
plot(N(2:6),Error(2:6),'*')
title('Error')
xlabel('Nodes')
ylabel('Error')