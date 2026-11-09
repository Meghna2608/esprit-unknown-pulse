function [arr_tau_hat, arr_eigenvalues] = Esprit(Y_hat, idx1, idx2, K, Gamma)
% ESPRIT estimate of K pulse locations.
%   Y_hat      : |Omega| x L noisy Fourier measurements
%   idx1, idx2 : rows of Y_hat forming the two sub-arrays (Omega_2 = Omega_1 + Gamma)
%   K          : number of pulses
%   Gamma      : shift between the two sub-arrays
% Returns the estimated locations in [0,1) and the eigenvalues of Psi.

L = size(Y_hat,2);
R_Y_hat = (1/L)*(Y_hat*Y_hat');            % empirical covariance matrix
[U_hat,D_hat] = eig(R_Y_hat);
[~,idx] = sort(diag(D_hat),'descend');
U_hat = U_hat(:,idx(1:K));                 % signal subspace

U1_hat = U_hat(idx1,:);
U2_hat = U_hat(idx2,:);
Psi = mldivide(U1_hat,U2_hat);             % least-squares ESPRIT

[~, D_psi] = eig(Psi);
arr_eigenvalues = diag(D_psi);
arr_tau_hat = -(1/(2*pi*Gamma))*atan2(imag(arr_eigenvalues),real(arr_eigenvalues));
arr_tau_hat = arr_tau_hat - floor(arr_tau_hat);

end
