%% Fig. 2(b): estimation error vs. number of Fourier measurements M.
%  Non-Dirac pulse (a = 0.75), maximum overlapping sub-arrays,
%  S = 5, L = 50, Delta = 0.05, SNR = 25 dB, 500 Monte Carlo trials.
close all; clc; clear;

T = 1;
num_reps = 500;
K = 5;
min_sep = 0.05;
snr = 25;
L = 50;
arr_tau = tau_k_gen(min_sep,K);
M_array = 10:10:100;
M_tilde_array = M_array;
rng(123); % fixes the amplitudes and the seeds matrix
X = ((1/sqrt(2))*(randn(K,L)+1i*randn(K,L)));
seeds = randi(2^31, length(M_array),num_reps);

error_locations_array_mo = zeros(numel(M_array),num_reps);
for mc_idx = 1:num_reps
    for m = 1:numel(M_array)
        M = M_array(m);
        M_tilde = M_tilde_array(m);
        size_Omega_tilde = 2*M_tilde;
        Omega_tilde = (1/T)*(0:1:(size_Omega_tilde-1)).';
        [g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,0.75);
        rng(seeds(m,mc_idx));

        Y_tilde = diag(g_tilde)*Phi_tilde*X;
        sigma_square = (norm(Y_tilde,'fro')^2*10^(-snr/10))/(M*L);
        sigma = (sigma_square)^0.5;
        Z = sigma.*((1/sqrt(2))*(randn(size(Y_tilde))+1i*randn(size(Y_tilde))));
        Y_tilde_noisy = Y_tilde + Z;

        [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,1);
        Y_hat = Y_tilde_noisy(sel,:);
        arr_tau_hat = Esprit(Y_hat,idx1,idx2,K,1/T);
        error_locations_array_mo(m,mc_idx) = location_error(arr_tau,arr_tau_hat,T);
    end
end
median_error = median(error_locations_array_mo,2);

%% Plot
figure(1);
plot(2*M_array,median_error,'LineWidth',1.5)
grid on;
xlabel('$$M$$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('estimation error','interpreter','LaTex','fontsize',15);
