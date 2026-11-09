%% Fig. 2(a): estimation error vs. number of snapshots L.
%  Non-Dirac pulse (a = 0.75), maximum overlapping sub-arrays,
%  S = 5, Delta = 0.05, M = 40, 500 Monte Carlo trials.
close all; clc; clear;

use_saved = true;    % true: plot the saved results in data/; false: rerun the simulation

T = 1;
num_reps = 500;
K = 5;
min_sep = 0.05;
M = 1/min_sep;
L_array = 6:6:60;
SNR_array = 10:10:30;

if use_saved
    load('data/Fig2a_error_vs_L.mat','median_error_location_mo');
else
    arr_tau = tau_k_gen(min_sep,K);
    M_tilde = 1/min_sep;
    size_Omega_tilde = 2*M_tilde;
    Omega_tilde = (1/T)*(0:1:(size_Omega_tilde-1)).';
    rng(123); % fixes the seeds matrix
    seeds = randi(2^31, length(L_array),length(SNR_array), num_reps);

    [g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,0.75);
    error_locations_array_mo = zeros(numel(L_array),numel(SNR_array),num_reps);

    for l = 1:numel(L_array)
        L = L_array(l);
        for s = 1:numel(SNR_array)
            snr = SNR_array(s);
            for mc_idx = 1:num_reps
                rng(seeds(l,s,mc_idx));
                X = ((1/sqrt(2))*(randn(K,L)+1i*randn(K,L)));
                Y_tilde = diag(g_tilde)*Phi_tilde*X;
                sigma_square = (norm(Y_tilde,'fro')^2*10^(-snr/10))/(M*L);
                sigma = (sigma_square)^0.5;
                Z = sigma.*((1/sqrt(2))*(randn(size(Y_tilde))+1i*randn(size(Y_tilde))));
                Y_tilde_noisy = Y_tilde + Z;

                [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,1);
                Y_hat = Y_tilde_noisy(sel,:);
                arr_tau_hat = Esprit(Y_hat,idx1,idx2,K,1/T);
                error_locations_array_mo(l,s,mc_idx) = location_error(arr_tau,arr_tau_hat,T);
            end
        end
    end
    median_error_location_mo = median(error_locations_array_mo,3);
end

%% Plot
median_error_snr_30 = median_error_location_mo(:,3);
median_error_snr_20 = median_error_location_mo(:,2);
median_error_snr_10 = median_error_location_mo(:,1);
figure(1);
plot(L_array,median_error_snr_30,'LineWidth',1.5);
hold on;
plot(L_array,median_error_snr_20,'LineWidth',1.5);
plot(L_array,median_error_snr_10,'LineWidth',1.5);
hold off;
xlabel('$$L$$', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('estimation error','interpreter','LaTex','fontsize',15);
legend('SNR = 30','SNR = 20','SNR = 10','fontsize',9)
xlim([6,60])
grid on;
