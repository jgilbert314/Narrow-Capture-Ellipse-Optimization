function [h1]= plot_traps_sph(X_spher,h, use_h, color, size) %X_spher=row(theta phi)
% plots traps at 3 columns of X_Cart_coord.
%If figure with handle h is to be used, use_h=1.
%color is trap color.
%size  is trap size.

	Ntraps2=numel(X_spher);
    Ntraps=Ntraps2/2;

    theta_phi=reshape(X_spher,Ntraps2,1);

    theta=theta_phi(1:Ntraps);
    phi=theta_phi(Ntraps+1:Ntraps2);
    
    X_Cart_coord=[sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta)];
    
    if(use_h==1)
        h1=figure(h);
    else
        h1=figure;
    end;

    if(use_h~=1)
        [X,Y,Z] = sphere;
        mesh(X,Y,Z,'facealpha',0.8);
        hold on;
        %colormap([0 0 0]);
        colormap('jet');
        axis equal;
        axis off;
        zoom(1.8);
    end;
    hold on;
    plot3(X_Cart_coord(:,1),X_Cart_coord(:,2),X_Cart_coord(:,3),'.','Markersize',size,'markeredgecolor',color);
    hold off;
    
end