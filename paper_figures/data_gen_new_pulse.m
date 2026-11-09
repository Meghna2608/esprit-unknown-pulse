function [g_tilde,Phi_tilde] = data_gen_new_pulse(Omega_tilde,arr_tau,a)
% Pulse spectrum g_tilde = G(Omega_tilde) and Phi_tilde = exp(-1j*2*pi*Omega_tilde*tau).
% G is the Fourier transform of g(t) = 4*rect(t/w)*cos^2(40*pi*a*t), w = pi/(40*a).

c1 = (pi/(40*a));
g_tilde = c1*(sinc((Omega_tilde-(40*a))*c1)+sinc((Omega_tilde+(40*a))*c1)+(2*sinc(Omega_tilde*c1)));
Phi_tilde = exp(-1i*2*pi*Omega_tilde*arr_tau);

end
