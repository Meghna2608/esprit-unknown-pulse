%% Fig. 2(d): estimation error vs. condition number of the amplitude covariance R_X.
%  Non-Dirac pulse (a = 0.75), maximum overlapping sub-arrays,
%  S = 5, L = 40, Delta = 0.05, M = 40, SNR = 25 dB, 500 Monte Carlo trials.
close all; clc; clear;

T = 1;
L = 40;
num_reps = 500;
K = 5;
min_sep = 0.05;
M = 1/min_sep;
SNR = 25;
arr_tau = tau_k_gen(min_sep,K);
M_tilde = 1/min_sep;
size_Omega_tilde = 2*M_tilde;
Omega_tilde = (1/T)*(0:1:(size_Omega_tilde-1)).';
condition_num_array = linspace(1,100,10);
X = (1/sqrt(2)) * (randn(K, L) + 1i*randn(K, L));
[U, S, V] = svd(X);
rng(123); % fixes the seeds matrix
seeds = randi(2^31,num_reps,numel(condition_num_array));

[g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,0.75);
error_locations_array_mo = zeros(num_reps,numel(condition_num_array));
for mc_idx = 1:num_reps
    for c = 1:numel(condition_num_array)
        % amplitudes with kappa(R_X) = condition_num_array(c)
        cond_num = condition_num_array(c);
        temp = linspace(1,cond_num,K);
        S_diag = sqrt(temp);
        S_new = [diag(S_diag) zeros(K,L-K)];
        X_new = U * S_new * V';

        rng(seeds(mc_idx,c));
        Y_tilde = diag(g_tilde)*Phi_tilde*X_new;
        sigma_square = (norm(Y_tilde,'fro')^2*10^(-SNR/10))/(M*L);
        sigma = (sigma_square)^0.5;
        Z = sigma.*((1/sqrt(2))*(randn(size(Y_tilde))+1i*randn(size(Y_tilde))));
        Y_tilde_noisy = Y_tilde + Z;

        [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,1);
        Y_hat = Y_tilde_noisy(sel,:);
        arr_tau_hat = Esprit(Y_hat,idx1,idx2,K,1/T);
        error_locations_array_mo(mc_idx,c) = location_error(arr_tau,arr_tau_hat,T);
    end
end
median_error = median(error_locations_array_mo,1);

%% Plot
figure(1);
plot(condition_num_array,median_error,'LineWidth',1.5);
xlabel('$$\kappa(\mathbf{R}_{\mathbf{X}})$$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('estimation error','interpreter','LaTex','fontsize',22);
xlim([1,100])
grid on;
