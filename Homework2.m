%% Homework 2
close all;
clear all;

x_cont = linspace(0,1,1000);
x_interval = [0,1];
%% Plot ODE f'(x) = cos(pi x)^2 and its solution
syms f(x)
ODE = diff(f,x) == (cos(pi.*x)).^2;
bc = f(0) == 0;
sol = dsolve(ODE,bc);

figure;
yyaxis right
plot(x_cont,(cos(pi.*x_cont)).^2,'DisplayName',['ODE (f''(x) ='...
    'cos(x \pi)^2)']); hold on
xlabel('x');
ylabel('f''(x)');
yyaxis left
fplot(sol,x_interval,'DisplayName',['Analytical Solution '...
    '(f(0)=0, ' newline 'f(x)=^{x}/_{2}+sin(2\pix)/4\pi )']);
legend;
xlabel('x');
ylabel('f(x)');
%% Finite Difference Method

% Initialize Matrices
h = [0.2 0.1 0.01];
mat_size = 1./h;
A = cell(length(mat_size),1);
b = cell(length(mat_size),1);
f = cell(length(mat_size),1);
x = cell(length(mat_size),1);
fin_diff = cell(length(mat_size),1);

figure;
fplot(sol,x_interval,'DisplayName','Analytical Solution');hold on
% Populate A and b
for i = 1:length(mat_size)
    x{i} = linspace(0,1,mat_size(i))';
    A{i} = eye(mat_size(i));
    for j = 2:(mat_size(i))
        A{i}(j,j-1) = -1;
    end
    
    b{i} = double(h(i) .* ones(mat_size(i),1) .* (cos(pi .* x{i})).^2);
    
    % Evaluate Af = b
    f{i} = linsolve(double(A{i}),b{i});
    
    % Plot solution
    plot(x{i},[0; f{i}(1:end-1)],'DisplayName',['Finite Difference h = ' num2str(h(i))]);
end
legend({},'Location','northwest');
xlabel('x_i');
ylabel('f(x_i)');

%% 3D Finite Element Method

% Local x
x_local = 0:0.01:1;
L = [1,0,0; 1,0.5,0.25;1,1,1];
L_inv = inv(sym(L));

syms x
X = [1,x,x^2];
N = X*L_inv;

figure; %visualize cubic shape functions
for i = 1:length(N)
    plot(x_local, subs(N(i),x_local),'DisplayName',['Square Shape Function ' num2str(i)]); hold on;
end
xlabel('Local Coordinates x');
ylabel('N(x)');
%% Stiffness Matrix
k_stiff = zeros(length(N));

for i = 1:length(N)
    for j = 1:length(N)
        k_stiff(i,j) = double(int( diff(N(i), x).*diff(N(j), x), x, 0, 1 ));
    end
end

nodes = 3; % number of elements
Ke = double(k_stiff);
RHSe = int(N, x, 0, 1)';


L = 1; % domain length
k = 1; %for k = [100,10,1,0.1] %; % thermal conductivity,1
S = 10; %    for S = [100,10,1,0.1]% % source,100
l = L/nodes; % node spacing assuming equidistant nodes
dim = nodes*(nodes - 1) +1; %dimension of stiffness matrix
K = zeros(dim); %initialize stiffness matrix with zeros
RHS = zeros([dim, 1]); %initialize load vector with zeros

for i = 1:nodes %cycle through elements
    
    %compute indices of entries for a given element
    ind_start = 1+(i-1)*(nodes-1);
    ind_end = ind_start + (nodes-1);
    
    % update stiffness matrix with the stiffness matrix of a single element
    K(ind_start:ind_end, ind_start:ind_end) = K(ind_start:ind_end, ...
        ind_start:ind_end) + Ke;
    
    %update load vector
    RHS(ind_start:ind_end) = RHS(ind_start:ind_end)+RHSe;
    
end

K = (k/l)*K; %multiply by constants
F = S*l*RHS;
large_number = 1e6; %a trick not to exclude T(0)
K(1,1) = K(1,1)+large_number;
T_FEM = linsolve(K, F)% solve for T
%
x_coord_FEM = 0:l/(nodes-1):1; %coordinates of nodes
figure;
plot(x_coord_FEM, T_FEM,'ks','DisplayName',['FEM, S=',num2str(S),' k=',...
    num2str(k)]); hold on %plot FEM solution
x_coord_solution = 0:0.01:1; %a grid to compute direct solution on
solution = -S/(2*k)*x_coord_solution.^2+S/k*x_coord_solution; %values of the direct
plot(x_coord_solution, solution,'DisplayName', ['Solution, S=',num2str(S),' k=',...
    num2str(k)]); hold on %visualize the direct solution
xlabel('x');
ylabel('T(x)');
xlim([0 1.1]);
legend({},'Location','northwest');
%% Error Computations
error_FEM = T_FEM - (-S/(2*k)*x_coord_FEM.^2+S/k*x_coord_FEM)'; %compute FEM error
error_FEM_norm = norm(error_FEM) %norm of the error vector
figure;
plot(error_FEM,'DisplayName',['Error of the FEM, S=',num2str(S),' k=',...
    num2str(k)]');
legend({},'Location','northwest');
xlabel('Elements');
ylabel('Error between Solution and FEM');
%end
%end


