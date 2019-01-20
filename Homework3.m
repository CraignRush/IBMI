%% IBMI Homework 4
clear all;
close all;
%% Assignment 1 Newton's Method
f = @(x) ((x.^2)./2) - sin(x);
f_der = @(x) x - cos(x);
x_prev = 0.5;
e = 1;
iterations = 1;

while e >= 1e-5
    x_next = x_prev - (f(x_prev)/f_der(x_prev));
    e = abs(x_next - x_prev);
    x_prev = x_next;
    iterations = iterations + 1;
end

%% Assignment 2 An imaging problem

%% a
s.orig.D = [-1,1;1,1;1,-1;-1,-1];
s.orig.E = [-0.5,0.5;0,0;0.5,0.5];
s.orig.I_e = [4,3,1];
s.orig.A = zeros(length(s.orig.D),length(s.orig.E));

for i = 1:length(s.orig.D)
    for j = 1:length(s.orig.E)
        s.orig.A(i,j) = 1./sqrt(sum((s.orig.D(i,:)-s.orig.E(j,:)).^2));
    end
end

%% b
s.orig.I_d = s.orig.A * s.orig.I_e';

%% c
tol = 1e-5;
[s.orig.moore.I_e, s.orig.moore.A_inv] = SolPseudoInvMoore(s.orig.A,s.orig.I_d,tol);
%% d
[s.orig.SVD.I_e, s.orig.SVD.A_inv] = SolPseudoInvSVD(s.orig.A,s.orig.I_d,0);
%% e
s.orig.LSQR.I_e = lsqr(s.orig.A,s.orig.I_d);

%% f
interval = [-1,1];
noise = interval(1) + (interval(2) - interval(1)) * rand(length(s.orig.I_d),1);

s.noise.I_d = s.orig.I_d + noise;

[s.noise.moore.I_e,s.noise.moore.A_inv] = SolPseudoInvMoore(s.orig.A,s.noise.I_d,tol);
[s.noise.SVD.I_e, s.noise.SVD.A_inv] = SolPseudoInvSVD(s.orig.A,s.noise.I_d,0);


% Compare Results with detectors manifold times farer away
factor = 3;
s.orig.moved.D = s.orig.D * factor;
for i = 1:length(s.orig.D)
    for j = 1:length(s.orig.E)
        s.orig.moved.A(i,j) = 1./sqrt(sum((s.orig.moved.D(i,:)-s.orig.E(j,:)).^2));
    end
end
s.orig.moved.I_d = s.orig.moved.A * s.orig.I_e';

% Conditioning
cond_num = cond(s.orig.A); %7.7273 --> ill-conditioned
cond_num_hat = cond(s.orig.moved.A); % 25 --> ill-conditioned

%% g
distance = [0.1, 0.1];
s.orig.added.E = [s.orig.E;s.orig.E(3,:) + distance;s.orig.E(3,:) - distance];
for i = 1:length(s.orig.D)
    for j = 1:length(s.orig.added.E)
        s.orig.added.A(i,j) = 1./sqrt(sum((s.orig.D(i,:)-s.orig.added.E(j,:)).^2));
    end
end

%% h
s.orig.added.I_e = s.orig.I_e';
s.orig.added.I_e = [s.orig.added.I_e; 1; 1];

s.orig.added.I_d = s.orig.added.A * s.orig.added.I_e;

noise = interval(1) + (interval(2) - interval(1)) * rand(length(s.orig.added.I_e),1);
s.noise.added.I_e = s.orig.added.I_e + noise;
s.noise.added.I_d = s.orig.added.A * s.noise.added.I_e;

%% i
s.orig.added.SVD.truncation = 0:1:2;
s.orig.added.SVD.I_e = cell(1,length(s.orig.added.SVD.truncation));
s.orig.added.SVD.A = cell(1,length(s.orig.added.SVD.truncation));

for i = s.orig.added.SVD.truncation
  [s.orig.added.SVD.I_e{i+1},s.orig.added.SVD.A{i+1}] = SolPseudoInvSVD(s.orig.added.A,s.orig.added.I_d,i);
  [s.noise.added.SVD.I_e{i+1},s.noise.added.SVD.A{i+1}] = SolPseudoInvSVD(s.orig.added.A,s.noise.added.I_d,i);  
end

%% Tikhonov Regularization
M = s.orig.added.A;
L = diag([2,1,0.5,0.1]);
s_tik = s.orig.added.I_d;

q_lam = @(lam) (M'*M+lam*L'*L)\M'*s_tik;

mq_norm = @(lam) norm(M*q_lam(lam)-s_tik,2);
lq_norm = @(lam) norm(L*q_lam(lam),2);
q_norm = @(lam) norm(q_lam(lam),2);

range = 1e-4 * 2.^(0:1:14);

x = arrayfun(mq_norm,range);
y = arrayfun(lq_norm,range);

figure;
plot(x,y,'x-');





