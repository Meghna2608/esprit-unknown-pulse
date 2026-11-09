%% ESPRIT for pulses of unknown shape from multiple snapshots
% M. Kalra and K. Lee, IEEE Trans. Signal Process., vol. 72, pp. 2637-2648, 2024.
% Model: Y_hat = G*Phi*X + Z, with G = diag(G(omega)) (unknown pulse),
%        Fourier matrix Phi = exp(-1j*2*pi*omega*tau), X = S x L amplitudes.
%
% Data are generated on the fine grid Omega_tilde = (0:2*M_tilde-1)/(2T):
% odd indices are the grid l/T, even indices are l/T + 1/(2T).
%   mode 1, max overlap (Sec. III):   Omega_1 = {0..M-1}/T, Omega_2 = Omega_1 + 1/T,
%                                     M = 2*M_half-1, Gamma = 1/T
%   mode 2, random doublets (Sec. IV): M_half random points l/T,
%                                     Omega_2 = Omega_1 + 1/(2T), Gamma = 1/(2T)
% Both designs use 2*M_half measurements. T = 1 is assumed.

clear; clc; rng(1);
T = 1;  L = 50;

% (a) Max overlap: S = 5, Delta = 0.05, 40 measurements, SNR = 25 dB
S = 5;  Delta = 0.05;  M_half = 20;  M_tilde = 40;  SNR = 25;  a = 0.75;
tau = [0, cumsum(Delta + (Delta/10)*rand(1, S-1))];
Y_hat_tilde = generate_data(M_tilde, T, tau, L, SNR, a);
[Omega_idx, Pi1, Pi2] = design_subarray(M_half, M_tilde, 1, 1);
tau_hat = Esprit(Y_hat_tilde(Omega_idx,:), Pi1, Pi2, S, 1/T);
fprintf('(a) max overlap:     md = %.2e\n', md(tau, tau_hat, T));

% (b) Random doublets: S = 7, Delta = 0.008, 40 measurements, SNR = 28 dB
S = 7;  Delta = 0.008;  M_half = 20;  M_tilde = 128;  SNR = 28;  a = 1;
tau = [0, cumsum(Delta + (Delta/10)*rand(1, S-1))];
Y_hat_tilde = generate_data(M_tilde, T, tau, L, SNR, a);
[Omega_idx, Pi1, Pi2] = design_subarray(M_half, M_tilde, 1, 2);
tau_hat = Esprit(Y_hat_tilde(Omega_idx,:), Pi1, Pi2, S, 1/(2*T));
fprintf('(b) random doublets: md = %.2e\n', md(tau, tau_hat, T));


%% ---------------------------- Functions --------------------------------

function tau_hat = Esprit(Y_hat, Pi1, Pi2, S, Gamma)
% ESPRIT. Y_hat: |Omega| x L data; Pi1, Pi2: rows of Omega_1, Omega_2
% (Omega_2 = Omega_1 + Gamma); S: number of pulses.
    L = size(Y_hat, 2);
    R_Y_hat = (1/L) * (Y_hat * Y_hat');               % empirical covariance
    [U_hat, D_hat] = eig(R_Y_hat);
    [~, order] = sort(real(diag(D_hat)), 'descend');
    U_hat = U_hat(:, order(1:S));                     % signal subspace
    Psi_hat = U_hat(Pi1,:) \ U_hat(Pi2,:);            % pinv(U1_hat) * U2_hat
    tau_hat = -angle(eig(Psi_hat).') / (2*pi*Gamma);  % -arg(lambda)/(2*pi*Gamma)
    tau_hat = tau_hat - floor(tau_hat);               % wrap to [0, 1)
end

function [Omega_idx, Pi1, Pi2] = design_subarray(M_half, M_tilde, mc_idx, mode)
% Omega_idx: indices of Omega in Omega_tilde; Pi1, Pi2: indices of Omega_1,
% Omega_2 in Omega. 
    switch mode
        case 1  % max overlap: shift 1/T
            Omega1_idx = 2*(1:2*M_half-1) - 1;
            Omega2_idx = Omega1_idx + 2;
        case 2  % random doublets: shift 1/(2T)
            % Omega_1: M_half distinct points sampled without replacement.
            rng(M_half * mc_idx);
            tmpv = randperm(M_tilde);
            Omega1_idx = 2*tmpv(1:M_half) - 1;
            Omega2_idx = Omega1_idx + 1;
        otherwise
            error('mode must be 1 or 2');
    end
    Omega_idx = union(Omega1_idx, Omega2_idx);
    [~, Pi1] = intersect(Omega_idx, Omega1_idx);
    [~, Pi2] = intersect(Omega_idx, Omega2_idx);
end

function err = md(tau, tau_hat, T)
% Matching distance: min over permutations of max_k |tau_k - tau_hat_pi(k)|_T.
    tau = tau(:).';  tau_hat = tau_hat(:).';
    P = perms(tau_hat);
    tmperr = abs(repmat(tau, size(P,1), 1) - P);
    tmperr = min(tmperr, T - tmperr);
    err = min(max(tmperr, [], 2));
end

function Y_hat_tilde = generate_data(M_tilde, T, tau, L, SNR, a)
% Noisy data on Omega_tilde. G is the Fourier transform of
% g(t) = 4*rect(t/w)*cos^2(40*pi*a*t), w = pi/(40*a).
    S = numel(tau);
    Omega_tilde = (0:2*M_tilde-1).' / (2*T);
    c1 = pi/(40*a);
    G_tilde = c1*( sinc((Omega_tilde - 40*a)*c1) + sinc((Omega_tilde + 40*a)*c1) ...
                 + 2*sinc(Omega_tilde*c1) );
    Phi_tilde = exp(-1j*2*pi*Omega_tilde*tau(:).');
    X = (randn(S,L) + 1j*randn(S,L)) / sqrt(2);
    Y_tilde = diag(G_tilde) * Phi_tilde * X;
    sigma = sqrt(norm(Y_tilde,'fro')^2 * 10^(-SNR/10) / (numel(Y_tilde)/2));
    Z = sigma * (randn(size(Y_tilde)) + 1j*randn(size(Y_tilde))) / sqrt(2);
    Y_hat_tilde = Y_tilde + Z;
end