%% Assignment 2
%Assume the 2x2 discretization of the spatial distribution of absorption coefficient shown in Fig. 1 and
% 2 point detectors with spacing L.
clear all;
close all;
%% 1) Design a model matrix A that relates X-ray measurements 
% p=(P1,P2,P3,P4,P5,P6)T at 3 shown angles to 
% the (unknown) absorption mu=(mu1,mu2,mu3,mu4)T as in: A*mu = p. Show A.

% scaling for 45° beams
a = sqrt(2) - 1; 

syms L;

% setup model matrix
A = [L, 0, L, 0;
0, L, 0, L;
a*L, 0, L, a*L;
a*L, L, 0, a*L;
0, 0, L, L;
L, L, 0, 0];

pretty(A);
%% 2) Assume a specific distribution (values) of mu_test and a specific value of L
% Simulate the corresponding measurements mu_test for this distribution
% using the model matrix A. Show p_test.

L_spec = 4;
A_subs = double(subs(A,L,L_spec));
mu_test = rand([4,1]); % assume (random) values of absorption mu

p_test = A_subs*mu_test % simulate corresponding measurements b

%% 3) Using the simulated measurements p_test, reconstruct absorption mu_rec,
% i.e. solve A*mu_rec=p_test. Show both the assumed (mu_test) and
% the reconstructed (mu_rec) absorption distributions.

% now we try to recover mu back from our measurements:
% A*mu = p => mu = inv(A)*p - we try to solve for mu that is assumed unknown
% inv(A)*p % doesn't work, rank of A is 3!

mu_rec = lsqr(A_subs, p_test); % min||A*mu-p||.^2 - we try to find a solution using minmization procedure
rel_error_perc = abs(mu_test - mu_rec)./mu_test*100 % looking at the relative error in percent, we're pretty far off from the real values


disp('4) Shortly interpret your result:');
disp('a. Does the reconstruction mu_rec correspond to the assumed distribution mu_test of the absorption coefficient?');
disp('According to var rel_error_perc, mu_rec matches the array mu_test at a relative error');
disp('bigger than 1E-10 percent. This means, the arrays can considered to be the same.');
disp('b. Did you use inv() for inverting the model? Why?');
disp('No, because the model matrix is no square matrix.');