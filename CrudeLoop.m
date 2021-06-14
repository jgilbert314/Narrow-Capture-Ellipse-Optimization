clear;

N_0 = 10;

eps_inds = 1:N_0; % Trap indices (ind > N_0-1 are larger traps)
for itr = 1:N_0
    TestingDiffTraps;
    eps_inds = [eps_inds(end), eps_inds(1:end-1)]; % Change position of larger traps
end