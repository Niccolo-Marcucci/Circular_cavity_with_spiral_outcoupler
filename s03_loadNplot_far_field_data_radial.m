
clear,
% close all
addpath('Lumerical-Objects/multilayer_design/functions');


% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
folder="SIM04_complex_outcouplers/far_field_data/";%sweep_scatter_size/far_field_data/"; %
% folder="SIM04_complex_outcouplers/sweep_start_radius/far_field_data/"; %

folder="SIM05_metasurface_outcoupler/far_field_data/";%";%
n_g=5;
top_charge=0;

name = "_TM_AlOTiO2_N10positive_filled_Dphi60_N12_sigma-1_charge_1";

fig = figure('units','normalized','outerposition',[0 0 1 1]);
i=0;

E_diff = zeros(1,20);
ER_max = zeros(1,20);
EL_max = zeros(1,20);
% for radial_FF  = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
% for lateral_FF = [0.8, 1]%, 1]
for j = 1%start_radius = linspace(1, 4, 20)*1e3
    i = i+1;
% 	details =strcat('_design_gd3_onSiO2_positive_filled_Dphi60_N12_sigma-1_m1_radialFF',string(radial_FF),'_lateralFF',string(lateral_FF), '.mat');
% 	details =strcat('_design_gd3_onSiO2_positive_filled_Dphi60_N12_sigma-1_m1_diameter',string(round(start_radius*2)), '.mat');
    load(strcat(folder,"far_field_data",name));
    % convert to matlab reference frame
    Ex=transpose(Ex);
    Ey=transpose(Ey);
    
    [Ux,Uy]=meshgrid(ux,uy);
    Ux=Ux';
    Uy=Uy';
%     E_phi = transpose(E_phi);
%     E_theta = transpose(E_tetha);
%     
%     % since cos(theta) = uz
%     theta = real( acos( sqrt(1 - ux.^2 - uy.^2)));
%     cos_phi = ux./sin(theta);
%     sin_phi = uy./sin(theta);
%     
%     % compute Ex Ey from Etheta and Ephi
%     Ex = E_theta.*cos_phi- E_phi.*sin_phi;
%     Ey = E_theta.*sin_phi+ E_phi.*cos_phi;
    
%     % add Ez if it not negligible
%     Ex = Ex + Ez./cos_phi;
%     Ey = Ey + Ez./sin_phi;
% 
%     Extmp=
%     Ex=Ex+Ey*exp(-1i*pi/2);
%     Ey=Ey-Extmp*exEx;
%     Ex=Ex+Ey;
%     Ey=Ey-Extmp;
%     Ex=Ex+Ey*exp(1i*pi/2);
%     Ey=Ey-Extmp*exp(1i*pi/2);
    
    ER = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(-1i*pi/2);
    EL = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);
    E_min = min( min(min(abs(ER).^2)), min(min(abs(EL).^2)));
    E_max = max( max(max(abs(ER).^2)), max(max(abs(EL).^2)));
    S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
    S0 = (abs(Ex).^2+abs(Ey).^2);
    chi = 0.5*asin( real(S3)./S0);
    
    R = sqrt(Ux.^2 +Uy.^2);
    THETA = atan2(Uy,Ux);
    
    N = 200;
    r_grid = linspace(0,1,N+1); % max(max(R)),N+1);
    theta_grid = linspace(0,2*pi,N+1);
    beta = r_grid(1:end-1);
    ER_r = zeros(1,N);
    EL_r = zeros(1,N);
    counter = zeros(1,N);
    for j = 1:length(R(:))
        r = R(j);
        idx = find((r <= r_grid(2:end)) .* (r > r_grid(1:end-1)));
        if ~isempty(idx)
            ER_r(idx) = ER_r(idx) + abs(ER(j)).^2;
            EL_r(idx) = EL_r(idx) + abs(EL(j)).^2;
            counter(idx) = counter(idx) + 1;
        end
        
    end
    ER_r = ER_r ./ counter.*(2*pi*beta);
    EL_r = EL_r ./ counter.*(2*pi*beta);
%     [R_grid, THETA_grid] = meshgrid(r_grid, theta_grid);
%     
%     E_fit = scatteredInterpolant(R(:), THETA(:), abs(ER(:)).^2,'nearest','nearest');
%     ER2 = E_fit(R_grid, THETA_grid);
%     
%     E_fit = scatteredInterpolant(R(:), THETA(:), abs(EL(:)).^2,'nearest','nearest');
%     EL2 = E_fit(R_grid, THETA_grid);
    
    %%
%     sgtitle('Period: 212nm - \Delta\Phi: \pi/3','fontsize',18,'fontweight','bold');
%     period = 212;
    subplot(1,1,i)
    plot(beta, ER_r, beta, EL_r);
    nicePlot
%     title(strcat(string(round(radial_FF * period)),'nmX',string(round(lateral_FF * period)),'nm'),'fontsize',12)
    % title(strcat("diameter ", string(round(start_radius*2/1000,1))))
    xlabel('u_T') 
    ylabel('Average Energy on a circle')
    nicePlot;
    drawnow;
    
    ROI = beta<0.3;
    E_diff(i) = max(ER_r(ROI)) / (max(ER_r(ROI)) + max(EL_r(ROI))) * 100;
    ER_max(i) = max(ER_r(ROI));
    EL_max(i) = max(EL_r(ROI));
%     end
end
saveas(fig,strcat(folder,"radial_far_field_PLOT",name),'png')   
saveas(fig,strcat(folder,"radial_far_field_PLOT",name),'fig')   
stop;
%%
fig = figure
subplot(1,2,1)
d = linspace(1,4,20)*2;
plot(d, E_diff, d, 100-E_diff)
legend("RHC","LHC")
title('Energy ripartition')
xlabel('diameter [um]')
ylabel('%')
nicePlot
subplot(1,2,2)
plot(d, ER_max, d, EL_max)
legend("RHC","LHC")
xlabel('diameter [um]')
ylabel('max(E) [a.u.]')
nicePlot
saveas(fig,strcat(folder,"Energy_ripartition",details,'.png')   )
saveas(fig,strcat(folder,"Energy_ripartition",details,'.fig')   )


function plot_surf(ux,uy,quantity,map,picture_title,symmetry,massimo)
%     s=surface(ux,uy,quantity);
%     s.EdgeColor='none';
    imagesc(ux,uy,quantity);
    set(gca,'YDir','normal') 
    ax=gca;
    if nargin > 3
        colormap(ax,map)
    if nargin > 4
        title(picture_title)
    if nargin > 5 
        if symmetry == "symmetric"
            c = max(abs([min(min(quantity)),max(max(quantity))]));
            caxis([-c c])
        elseif symmetry == "unitary"
            caxis([0 1])
        else
            error('Wrong "symmetry" signment')
        end
    if nargin > 6
        if symmetry == "symmetric"
            caxis([-1 1]*massimo)
        elseif symmetry == "unitary"
            caxis([0 1]*massimo)
        else
            error('Wrong "symmetry" signment')
        end
    end;end;end;end
    colorbar
    xlabel("ux");
    ylabel('uy');
    axis('square')
end