function [conc_AB_nM, fraction_A_bound, fraction_B_bound] = bimolecular_reaction_solver(kD_nM, conc_A_0_nM, conc_B_0_nM)
%% calculate equilibrium concentrations of a bimolecular reaction A+B->AB
%   INPUTS
%   kD_nM: dissociation constant in nanomolar units
%   conc_A_0_nM: Total concentration of A species in nanomolar units
%   conc_B_0_nM: Total concentration of B species in nanomolar units

% change units from nanomolar to molar
kD = kD_nM / 10^9;
conc_A_0 = conc_A_0_nM / 10^9;
conc_B_0 = conc_B_0_nM / 10^9;

% calculate association constant from dissocation constant
kA = 1/kD;

% solve quadratic equation
a = kA;
b = -(1 + conc_A_0 * kA + conc_B_0 * kA);
c = conc_A_0 * conc_B_0 * kA;

conc_AB = (-b - sqrt(b^2 - 4 * a * c)) / (2 * a);
conc_AB_nM = conc_AB * 10^9;

fraction_A_bound = conc_AB_nM / conc_A_0_nM;
fraction_B_bound = conc_AB_nM / conc_B_0_nM;
end