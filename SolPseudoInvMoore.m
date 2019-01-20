function [I_e_moore,A_inv_moore] = SolPseudoInvMoore(A,b,tol)
A_inv_moore = pinv(A,tol);
I_e_moore = A_inv_moore * b;
end

