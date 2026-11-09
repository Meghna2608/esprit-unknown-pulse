%% Fig. 1 (a): phase transition of the estimation error, L vs. 10^(-SNR/10),
%  Dirac pulse, maximum overlapping sub-arrays.
%  S = 5, M = 240, Delta = 0.0125, SNR = 10..20 dB, 500 Monte Carlo trials.
close all; clc; clear;

use_saved = true;    % true: plot the saved results in data/; false: rerun the simulation

T = 1;
num_reps = 500;
K = 5;
min_sep = 0.0125;
M = 1.5/min_sep;
L_array = 6:3:78;
SNR_array = 20:-1:10;

if use_saved
    load('data/Fig1a_phase_transition_Dirac.mat','median_error_location_mo');
else
    arr_tau = tau_k_gen(min_sep,K);
    M_tilde = 1.5/min_sep;
    size_Omega_tilde = 2*M_tilde;
    Omega_tilde = (1/T)*(0:1:(size_Omega_tilde-1)).';
    rng(123); % fixes the seeds matrix
    seeds = randi(2^31, length(L_array),length(SNR_array), num_reps);

    [~,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,0.75);
    g_tilde = ones(size(Omega_tilde));   % Dirac pulse: G(omega) = 1
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
figure(1);
imshow(kron(flip(median_error_location_mo,1),ones(20,30)),[])
set(gca, 'visible','on')
yticks(2*(5:40:250))
yticklabels((78:-12:6))
xticks(3*(5:30:130))
xticklabels([0.01,0.02,0.04,0.08])
colorbar
colormap('jet');
set(gca,'fontsize',15);
xlabel('$10^{-SNR/10}$','interpreter','LaTex','fontsize',15);
ylabel('$L$','interpreter','LaTex','fontsize',15);
