%% Fig. 2(c): estimation error vs. pulse width w = pi/(40*a).
%  Maximum overlapping sub-arrays, S = 5, L = 50, Delta = 0.05, M = 40,
%  SNR = 25 dB, 500 Monte Carlo trials.
clear; clc; close all;

L = 50;
T = 1;
K = 5;
min_sep = 0.05;
snr = 25;
Gamma = 1/T;

M_tilde = 100;
M = 1/min_sep;
rng(2608)
size_Omega_tilde = 2*M_tilde;
Omega_tilde = Gamma*(0:1:(size_Omega_tilde-1)).';
iter = 500;
arr_tau = tau_k_gen_spread(T,min_sep,K,123);
a_arr = linspace(0.2618,31.4159,10);
width_time_domain_arr = pi./(40*a_arr);
error_arr = zeros(numel(a_arr),iter);

for mc_idx = 1:iter
    for i = 1:numel(a_arr)
        a = a_arr(i);
        [g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,a);

        rng(mc_idx*K*L);
        X = ((1/sqrt(2))*(randn(K,L)+1i*randn(K,L)));
        Y_tilde = diag(g_tilde)*Phi_tilde*X;
        sigma_square = (norm(Y_tilde,'fro')^2*10^(-snr/10))/(M*L);
        sigma = (sigma_square)^0.5;
        rng(M*L*mc_idx*sigma);
        Z = sigma.*((1/sqrt(2))*(randn(size(Y_tilde))+1i*randn(size(Y_tilde))));
        Y_tilde_noisy = Y_tilde + Z;

        [sel,idx1,idx2] = design_subarray(M,M_tilde,mc_idx,1);
        Y_hat = Y_tilde_noisy(sel,:);
        arr_tau_hat = Esprit(Y_hat,idx1,idx2,K,Gamma);
        error_arr(i,mc_idx) = location_error(arr_tau,arr_tau_hat.',T);
    end
end
median_err = median(error_arr,2);

%% Plot
figure(1);
plot(width_time_domain_arr',median_err,'LineWidth',1.5)
axis([0.002 0.2 0 0.0144])
grid on;
xlabel('pulse width','interpreter','LaTex','fontsize',15);
ylabel('estimation error','interpreter','LaTex','fontsize',15);


%% ------------------------------------------------------------------------
function [seq,err] = tau_k_gen_spread(T,min_sep,K,mc_idx)
% K pulse locations with separations between min_sep and min_sep + (T-min_sep)/(K-1).
rng(10*min_sep*K*mc_idx);

err = 1;
max_iter = 50;
iter = 0;

while err == 1 && iter < max_iter

    iter = iter + 1;
    seperation = min_sep + [0 rand(1,K-1)*(T-min_sep)/(K-1)];
    seperation = seperation(randperm(K-1));
    seq = [0 cumsum(seperation)];

    tmpmat = abs(repmat(seq',[1 K]) - repmat(seq,[K 1]));
    tmpmat = min(tmpmat,1-tmpmat) + diag(ones(K,1));

    if min(tmpmat(:)) <= min_sep
        err = 0;
    end
end

if iter == max_iter
    error('tau_k_gen failed')
end

end
