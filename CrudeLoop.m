clear;

N_0 = 10;

eps_inds = 1:N_0;
for itr = 1:N_0
    TestingDiffTraps;
    eps_inds = [eps_inds(end), eps_inds(1:end-1)];
end