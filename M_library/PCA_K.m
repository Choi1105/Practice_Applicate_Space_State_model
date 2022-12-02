function [PCm, eigen_val, Vm ] = PCA_K(Y)
% Y = T by k
% PCm = T by k

% 1 단계: 상관계수행렬 계산하기
CORM = corrcoef(Y);

% 2 단계: 특성근, 특성벡터 계산하기
[V, D] = eig(CORM);

% V = 특성벡터, D의 대각값들 = 특성근
   
% 3 단계
eigen_val = diag(D); % 3 by 1
[eigen_val, index] = sort(eigen_val, 'descend'); % 특성근 정렬하기
Vm = V(:, index);% 특성벡터 정렬하기

% 4 단계: PC 계산하기
PCm = Y*Vm; % T by 3

disp(['특성근 = ', num2str(eigen_val'/cols(Y))]);
end