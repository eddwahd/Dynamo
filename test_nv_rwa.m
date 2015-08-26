function [Hp, Cp, Hp_slow, Cp_slow] = test_nv_rwa(theta)
% Basic NV RWA test.
%  [Hp, Cp] = test_nv_rwa(theta)
%
%  Transforms the electrons-only NV Hamiltonian into a rotating frame.
%  Control field is applied in the x direction.
%  theta is the angle between B0 and the NV symmetry axis.
%  Returns H'(t) and C'(t).

% Ville Bergholm 2014

dim = [3, 3];
ops;

TU = 1e-6; % time unit, in s

Delta = 2*pi * 2.87e9 * TU;
omega = -2*pi * 100e6 * TU;
omega_carrier = Delta +omega;

% Hamiltonians
H = Delta * Z^2 +cos(theta) * omega * Z +sin(theta) * omega * X;
freqs = eig(H)/2/pi
% Rotating frame Hamiltonian
%H0 = omega_carrier * Z;
H0 = omega_carrier * Z^2;

% Drift Hamiltonian
[Hp, Hp_slow, a_min, a_max, omega_rot] = RWA(H0, H-H0, 0, 1, true, 300)

% Control Hamiltonian
[Cp, Cp_slow, a_min, a_max, omega_rot] = RWA(H0, X, omega_carrier, 1, true, 300)
