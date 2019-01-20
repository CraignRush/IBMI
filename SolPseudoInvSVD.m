function [I_e_SVD,A_inv_SVD] = SolPseudoInvSVD(A,b,truncate)
[U,S,V] = svd(A);

[m,n] = size(U);
U = U(1:m,1:n-truncate);

[m,n] = size(V);
V = V(1:m,1:n-truncate);

[m,n] = size(S);
S = S(1:m-truncate,1:n-truncate);

S_degger = S;
S_degger(S_degger > 0) = 1./S_degger(S_degger > 0);
A_inv_SVD = V * S_degger' * U';
I_e_SVD = A_inv_SVD * b;
end

