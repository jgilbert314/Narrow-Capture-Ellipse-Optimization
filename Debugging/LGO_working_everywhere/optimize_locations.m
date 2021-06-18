function [Opt_angles] = optimize_locations(theta_phi,way,proc)%input: a matrix with 1 column: [all theta; all phi]
%Optimizes locations of given traps using selected alg

	Ntraps2=size(theta_phi,1);
	Ntraps=Ntraps2/2;

    theta=theta_phi(1:Ntraps);
    phi=theta_phi(Ntraps+1:Ntraps2);

    Ntraps=size(theta,1);

    %initial condition for solution
    opt_x0 = [theta;phi];

%    options = optimset('Display','iter');
    options = optimset('Display','off');

    % possibility 1: use fsolve to get local min; not implemented here
    if(way==1)
        Opt_angles = opt_x0;
    elseif(way==2)%LGO
        nvars = 2*Ntraps;
        ncons = 0;
        ctypes = [];
        min_all=zeros(nvars,1);
        max_all=[pi*ones(Ntraps,1);2*pi*ones(Ntraps,1)];
        
        add_vec=2*(pi-eps)*(opt_x0<0);
        opt_x0=opt_x0+add_vec; %correct phi-entries that got negative

        subtr_vec=(pi-eps)*(opt_x0>max_all);
        opt_x0=opt_x0-[subtr_vec(1:Ntraps) ; 2*subtr_vec(Ntraps+1:2*Ntraps)]; %correct phi-entries that got negative

        bnds=[min_all, opt_x0, max_all];
        min_fun_val=10000;
        %lgo_Options=[3, 10000, 1000, 1000, 1.0, min_fun_val, min_fun_val, 0.000001, 0.000001, 0.0000001, 0, 300];%50->101: 6s, LGO_Coul=4559...
        
        %[lg_Optimal, lg_X, lg_CViol, lg_AdtlInfo] = matlablgo(proc, nvars, bnds, ncons ,ctypes, lgo_Options);
        [lg_Optimal, lg_X, lg_CViol, lg_AdtlInfo] = matlablgo(proc, nvars, bnds, ncons ,ctypes);
        
        Opt_angles = reshape(lg_X, Ntraps2,1);
    end
    
end