%% Fig. 3(b): phase transition of the estimation error, M vs. S, random doublets.
%  Non-Dirac pulse (a = 1), L = 50, Delta = 0.008, M_tilde = 130 (fine grid of
%  2*M_tilde = 260 points), SNR = 26 dB, 500 Monte Carlo trials.
close all; clc; clear;

use_saved = true;    % true: plot the saved results in data/; false: rerun the simulation

L = 50;
T = 1;
snr = 26;
num_reps = 500;
min_sep = 0.008;
M_array = 8:4:72;    % number of random doublets (2*M measurements)
K_array = 2:1:8;
M_tilde = 130;

if use_saved
    load('data/Fig3b_phase_transition_M_vs_S.mat','median_error_location_rd');
else
    arr_tau_k = tau_k_gen(min_sep,8);
    rng(123); % fixes the seeds matrix
    seeds = randi(2^31, length(M_array),length(K_array), num_reps);
    Omega_tilde = (0:1:2*M_tilde-1).'/(2*T);

    error_locations_array_rd = zeros(numel(M_array),numel(K_array),num_reps);
    for m = 1:numel(M_array)
        M = M_array(m);
        for k = 1:numel(K_array)
            K = K_array(k);
            arr_tau = arr_tau_k(1:K);       % first K of the 8 locations
            [g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,1);
            for mc_idx = 1:num_reps
                rng(seeds(m,k,mc_idx));
                X = ((1/sqrt(2))*(randn(K,L)+1i*randn(K,L)));
                Y_tilde = diag(g_tilde)*Phi_tilde*X;
                sigma_square = (norm(Y_tilde,'fro')^2*10^(-snr/10))/(M_tilde*L);
                sigma = (sigma_square)^0.5;
                Z = sigma.*((1/sqrt(2))*(randn(size(Y_tilde))+1i*randn(size(Y_tilde))));
                Y_tilde_noisy = Y_tilde + Z;

                [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,2);
                Y_hat = Y_tilde_noisy(sel,:);
                arr_tau_hat = Esprit(Y_hat,idx1,idx2,K,1/(2*T));
                error_locations_array_rd(m,k,mc_idx) = location_error(arr_tau,arr_tau_hat,T);
            end
        end
    end
    median_error_location_rd = median(error_locations_array_rd,3);
end

%% Plot
figure(1);
imshow(kron(flip(median_error_location_rd,1),ones(15,30)),[])
set(gca, 'visible','on')
yticks(1.5*(5:20:170))
yticklabels(2*(72:-8:8))
xticks(3*(5:10:70))
xticklabels(2:1:8)
colorbar
colormap('jet');
set(gca,'fontsize',20);
xlabel('$S$','interpreter','LaTex','fontsize',15);
ylabel('$M$','interpreter','LaTex','fontsize',15);
