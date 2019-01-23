%% IBMI Homework 4
clear all;
close all;
rng default;
%% Assignment 1 Newton's Method
syms f(x)
f(x) =  ((x.^2)./2) - sin(x);
df = diff(f,x);
ddf = diff(df,x);
x_prev = 0.5;
e = 1;
iterations = 1;
x = -0.5:0.001:1.5;

figure;
plot(x,double(f(x)),'LineWidth',1.2,'DisplayName','f(x)');
hold on
plot(x,double(df(x)),'LineWidth',1.2,'DisplayName','f''(x)');
l = refline([0 0]);
l.Color = [0,0,0];
l.LineStyle = '-';
l.LineWidth = 0.1;
%legend;
while e >= 1e-5
    x_next = x_prev - (double(df(x_prev))/double(ddf(x_prev)));
    e = abs(x_next - x_prev);    
    %l = refline(df(x_prev),df(x_prev)-ddf(x_prev)*x_prev);
    plot(x_prev,double(f(x_prev)),'kx');
    plot(x_prev,0,'rx');
    x_prev = x_next;
    iterations = iterations + 1;
end

xlabel('x');
ylabel('y');
print('fig/NewtonMethod.pdf','-dpdf');
%% Assignment 2 An imaging problem
% a
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
interval = [-0.1,0.1];
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

s.noise.moved.I_d = (s.orig.moved.A * s.orig.I_e') + noise;

[s.noise.moved.moore.I_e,s.noise.moved.moore.A_inv] = ...
    SolPseudoInvMoore(s.orig.moved.A,s.noise.moved.I_d,tol);
[s.noise.moved.SVD.I_e, s.noise.moved.SVD.A_inv] = ...
    SolPseudoInvSVD(s.orig.moved.A,s.noise.moved.I_d,0);

% Conditioning
cond_num = cond(s.orig.A); %7.7273 --> ill-conditioned
cond_num_hat = cond(s.orig.moved.A); % 25.01 --> ill-conditioned

%% g
distance = [0.1, 0.1];
s.orig.added.E = [s.orig.E;s.orig.E(3,:) + distance;s.orig.E(3,:) - distance];
for i = 1:length(s.orig.D)
    for j = 1:length(s.orig.added.E)
        s.orig.added.A(i,j) = 1./sqrt(sum((s.orig.D(i,:)-s.orig.added.E(j,:)).^2));
    end
end

%% h
s.orig.added.I_e = [s.orig.I_e'; 1; 1];

s.orig.added.I_d = s.orig.added.A * s.orig.added.I_e;

noise = interval(1) + (interval(2) - interval(1)) * rand(length(s.orig.added.I_e),1);
s.noise.added.I_e = s.orig.added.I_e + noise;
s.noise.added.I_d = s.orig.added.A * s.noise.added.I_e;

%% i
s.orig.added.SVD.truncation = 0:1:2;
s.noise.added.SVD.I_e = cell(1,length(s.orig.added.SVD.truncation));
s.noise.added.SVD.A = cell(1,length(s.orig.added.SVD.truncation));

Level = noise*[0,1,10];

for i = s.orig.added.SVD.truncation
    for j = 1:size(Level,2)
     s.noise.added.SVD.I_d = s.orig.added.A * (s.orig.added.I_e + Level(:,j));
  [s.noise.added.SVD.I_e{i+1},s.noise.added.SVD.A{i+1}] = SolPseudoInvSVD(s.orig.added.A,s.noise.added.SVD.I_d,i);  
    digits(4);
    sprintf('Truncation:%d, Noise Level:%d\n',i,Level(1,j)./noise(1))
%    latex(vpa(sym(s.noise.added.SVD.I_e{i+1})));
    end
end
%% Tikhonov Regularization
M = s.orig.added.A;
s_tik = s.noise.added.I_d;
range = 1e-4 * 2.^(0:1:14);
x = zeros(length(range));
y = zeros(length(range));

L = eye(5) + 0.5*randn(5);

q_lam = @(lam) (M'*M+lam*(L'*L))\M'*s_tik;

mq_norm = @(lam) norm(M*q_lam(lam)-s_tik,2);
lq_norm = @(lam) norm(L*q_lam(lam),2);
q_norm = @(lam) norm(q_lam(lam),2);

x = arrayfun(mq_norm,range);
y = arrayfun(lq_norm,range);

dy = [0,diff(y)];
ddy = [0, diff(dy)];
dddy = [0, diff(ddy)];
max_dist = max(dddy);

fig =figure;
plot(x,y,'x-','LineWidth',1.2,'DisplayName','L-curve');
hold on
legend;
xlabel('||Mq(\lambda) - s||');
ylabel('||Lq(\lambda)||');
pos_opt = find(max(dddy)== dddy)+2;
plot(x(pos_opt),y(pos_opt),'ro','MarkerSize',10,'DisplayName',...
'Optimal regularization parameter');
q_lam(x(pos_opt))
print(fig,'fig/Tikhonov.eps','-dpdf');




