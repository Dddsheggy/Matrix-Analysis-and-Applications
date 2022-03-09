clc
clear all
close all
%A simple example
% A = load('toy_matrix.txt');
A = rand(18,15);
k = 13;

% traversal algorithm
R = traversal_CSSP(A,k)
C = A(:,R);
loss_t = norm(A-C*pinv(C)*A, 'fro')^2


% our algorithm
R = approx_CSSP(A, k,1)
C = A(:,R);
loss_a = norm(A-C*pinv(C)*A, 'fro')^2


R = iter_fs(A, k);
C = A(:,R);
loss_i = norm(A-C*pinv(C)*A, 'fro')^2

% greedy algorithm
% R = RndGreedyCSS(A, k, 0)
% C = A(:,R);
% loss = norm(A-C*pinv(C)*A, 'fro')^2