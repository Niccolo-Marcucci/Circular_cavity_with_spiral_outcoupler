clear,close all
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
folder="SIM02_no_cavity_spiral_outcoupler/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
top_charge=2;
for n_g=[5]

load(strcat(folder,"far_field_data_positive_charge2_N5_finer_mesh"));
% load(strcat(folder,"far_field_data","_charge",string(top_charge),"_N",string(n_g)));
% load(strcat(folder,"far_field_data","_charge",string(top_charge)));
% 
% useful_ux = abs(ux)<0.2;
% useful_uy = abs(uy)<0.2;

[ux,uy]=meshgrid(ux,uy);
% ux = ux(useful_ux,useful_uy);
% uy = uy(useful_ux,useful_uy);
% E_phi = E_phi(useful_ux,useful_uy);
% E_theta = E_theta(useful_ux,useful_uy);


% cos(thneta) = uz
theta = real( acos( sqrt(1 - ux.^2 - uy.^2)));
cos_phi = ux./sin(theta);
sin_phi = uy./sin(theta);

Ex = E_theta.*cos_phi - E_phi.*sin_phi;
Ey = E_theta.*sin_phi + E_phi.*cos_phi;

EL = +sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);
ER = -sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);

S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
chi = 0.5*asin( real(S3)./(abs(Ex).^2+abs(Ey).^2));

% figure
% plot_surf(ux,uy,sin_phi,'hot',"something about e_phi");
% plot_surf(ux,uy,theta,'hot',"something about e_theta");

figure
subplot(1,2,1)
plot_surf(ux,uy,abs(Ey).^2   +abs(Ex).^2,'hot','Intensity from Ex and Ey');
subplot(1,2,2)
plot_surf(ux,uy,abs(E_phi).^2+abs(E_theta).^2,'hot','Intensity from E_\phi and E_\theta');
%%
figure
subplot(2,3,1)
plot_surf(ux,uy,abs(ER).^2,'hot',"Right circular polarization intensity");
subplot(2,3,2)
plot_surf(ux,uy,abs(EL).^2,'hot',"Left circular polarization intensity");
subplot(2,3,4)
plot_surf(ux,uy,angle(ER),'hsv',"Right circular polarization phase");
subplot(2,3,5)
plot_surf(ux,uy,angle(EL),'hsv',"Left circular polarization phase");

subplot(2,3,3)
plot_surf(ux,uy,real(S3),'jet',"S3 Stokes parameter");
subplot(2,3,6)
plot_surf(ux,uy,abs(tan(chi)),'hot',"eccentricity |tan\chi|");
end



function plot_surf(ux,uy,Quantity,map,picture_title)
    s=surf(ux,uy,Quantity);
    s.EdgeColor='none';
    ax=gca;
    colormap(ax,map)
    colorbar
    view(2)
    xlabel("ux");
    ylabel('uy');
    axis('square')
    title(picture_title)
end