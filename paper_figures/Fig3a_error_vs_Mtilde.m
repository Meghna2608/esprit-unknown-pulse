%% Fig. 3(a): estimation error vs. full sampling size M_tilde, random doublets.
%  Dirac and non-Dirac (a = 1) pulses, S = 7, L = 50, M = 40, Delta = 0.008,
%  SNR = 28 dB, 1000 Monte Carlo trials.
close all; clc; clear;

use_saved = true;    % true: plot the saved results in data/; false: rerun the simulation

L = 50;
T = 1;
snr = 28;
num_reps = 1000;
M = 20;              % number of random doublets (2*M = 40 measurements)
K = 7;
min_sep = 0.008;
M_tilde_critical = ceil(2+1/min_sep);
M_tilde_array = 2*ceil(M_tilde_critical*(0.5:0.1:3.9)/2);

if use_saved
    d = load('data/Fig3a_error_vs_Mtilde_nonDirac.mat','median_error_location_rd');
    median_error_nonDirac = d.median_error_location_rd;
    d = load('data/Fig3a_error_vs_Mtilde_Dirac.mat','median_error_location_rd');
    median_error_Dirac = d.median_error_location_rd;
else
    arr_tau = tau_k_gen(min_sep,K);
    pulses = {'nonDirac','Dirac'};
    for p = 1:2
        error_locations_array_rd = zeros(numel(M_tilde_array),num_reps);
        for n = 1:numel(M_tilde_array)
            M_tilde = M_tilde_array(n);
            Omega_tilde = (0:1:2*M_tilde-1).'/(2*T);
            [g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,1);
            if strcmp(pulses{p},'Dirac')
                g_tilde = ones(size(Omega_tilde));   % Dirac pulse: G(omega) = 1
            end
            for mc_idx = 1:num_reps
                rng(100*M*K*L*mc_idx);
                X = ((1/sqrt(2))*(randn(K,L)+1i*randn(K,L)));
                Y_tilde = diag(g_tilde)*Phi_tilde*X;
                sigma_square = (norm(Y_tilde,'fro')^2*10^(-snr/10))/(M_tilde*L);
                sigma = (sigma_square)^0.5;
                Z = sigma.*((1/sqrt(2))*(randn(size(Y_tilde))+1i*randn(size(Y_tilde))));
                Y_tilde_noisy = Y_tilde + Z;

                [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,2);
                Y_hat = Y_tilde_noisy(sel,:);
                arr_tau_hat = Esprit(Y_hat,idx1,idx2,K,1/(2*T));
                error_locations_array_rd(n,mc_idx) = location_error(arr_tau,arr_tau_hat',T);
            end
        end
        if p == 1
            median_error_nonDirac = median(error_locations_array_rd,2);
        else
            median_error_Dirac = median(error_locations_array_rd,2);
        end
    end
end

%% Plot
figure(1);
semilogy(M_tilde_array',median_error_nonDirac,'LineWidth',1.5);
hold on;
semilogy(M_tilde_array',median_error_Dirac,'LineWidth',1.5);
set(gca,'fontsize',20);
axis([64 250 0 0.35])
xlabel('$\tilde{M}$','interpreter','LaTex','fontsize',15);
ylabel('estimation error','interpreter','LaTex','fontsize',15);
legend('non-Dirac','Dirac','interpreter','LaTex','fontsize',15);
grid on
