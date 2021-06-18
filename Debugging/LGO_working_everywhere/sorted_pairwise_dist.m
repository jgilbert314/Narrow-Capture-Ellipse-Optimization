function Ret= sorted_pairwise_dist(theta_phi1) %input: LGO needs a ROW [all theta, all phi]. This proc does not care.

    Ntraps2=numel(theta_phi1);
    Ntraps=Ntraps2/2;

    theta_phi=reshape(theta_phi1,1,Ntraps2);

    thetas=theta_phi(1:Ntraps);
    phis=theta_phi(Ntraps+1:Ntraps2);
    
    theta_comb=nchoosek(thetas,2); %Note: combinations do not repeat, so we don't have to divide energy by 2
    phi_comb=nchoosek(phis,2);
    N_comb=size(theta_comb,1);

    ct=cos(theta_comb);
    st=sin(theta_comb);

    cp=cos(phi_comb);
    sp=sin(phi_comb);

    %under roots: q = sqrt(2)*    (   1 - sin(theta)*cos(phi)*sin(theta[i])*cos(phi[i])-sin(theta)*sin(phi)*sin(theta[i])*sin(phi[i]) -cos(theta)*cos(theta[i]) )   )  

    q =  sqrt(2)*sqrt(abs(1 - st(:,1).*st(:,2).*( cp(:,1).*cp(:,2) + sp(:,1).*sp(:,2)) - ct(:,1).*ct(:,2)));

    Ret=reshape(sort(q,'descend'),numel(q),1);
end


