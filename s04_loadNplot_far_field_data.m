clear,
close all
folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
top_charge=2;
for n_g=[5]

% load(strcat(folder,"far_field_data_design_gd3_onSiO2_2_N35_positive"));
load(strcat(folder,"far_field_data","_charge",string(top_charge),"_N",string(n_g)));
% load(strcat(folder,"far_field_data","_charge",string(top_charge)));
% load(strcat(folder,"far_field_dataRing_positive"));
% 
% useful_ux = abs(ux)<0.2;
% useful_uy = abs(uy)<0.2;

[ux,uy]=meshgrid(ux,uy);
% ux = ux(useful_ux,useful_uy);
% uy = uy(useful_ux,useful_uy);
% E_phi = E_phi(useful_ux,useful_uy);
% E_theta = E_theta(useful_ux,useful_uy);

E_theta2 = E_theta;
E_theta = E_phi;
E_phi = E_theta2;

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
plot_surf(ux,uy,angle(ER),'hot',"Right circular polarization phase");
subplot(2,3,5)
plot_surf(ux,uy,angle(EL),'hot',"Left circular polarization phase");

subplot(2,3,3)
plot_surf(ux,uy,real(S3),'hot',"S3 Stokes parameter",1);
subplot(2,3,6)
plot_surf(ux,uy,abs(tan(chi)),'hot',"eccentricity |tan\chi|");

% sgtitle(['Topological charge ',num2str(top_charge)],'fontsize',18,'fontweight','bold') 
% saveas(figure(2),strcat(folder,"far_field_PLOT","_charge",string(top_charge)),'png')
% close all
end



function plot_surf(ux,uy,quantity,map,picture_title,symmetric)
%     s=surf(ux,uy,Quantity);
%     s.EdgeColor='none';
%     view(2)
    imagesc(ux(1,:),uy(:,1),quantity);
    ax=gca;
    colormap(ax,map)
    if nargin > 5 && symmetric
        c = max(abs([min(min(quantity)),max(max(quantity))]));
        caxis([-c c])
    end
    colorbar
    xlabel("ux");
    ylabel('uy');
    axis('square')
    title(picture_title)
end