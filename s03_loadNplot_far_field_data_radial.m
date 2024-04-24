
clear,
% close all
addpath('Lumerical-Objects/multilayer_design/functions');


% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
folder="SIM04_complex_outcouplers/far_field_data/";%sweep_scatter_size/far_field_data/"; %
% folder="SIM04_complex_outcouplers/sweep_start_radius/far_field_data/"; %
folder="SIM05_metasurface_outcoupler/scatterTests_PMMA_topped_positive/far_field_data/";%";%

names = [];
for dphi = -60%-60:120:60
    % for sigma = 1%-1:2:1
    %     for charge = 0%-1:1
    i = 0;
    sigma = 1;
    charge = 0;
    for sc_width  = [25, 50, 75, 100, 125, 150]
        for sc_length = [250, 275, 300, 325]
            details =['_TM_AlOTiO2_N10positive_filled_scShapeI_Dphi',num2str(dphi),'_N12_sigma',num2str(sigma),'_charge', num2str(charge), '_scWidth', num2str(sc_width), '_scLength', num2str(sc_length)];
            % details = ['_TM_AlOTiO2_N10positive_filled_scShapeI_Dphi',num2str(dphi),'_N24_sigma',num2str(sigma),'_charge', num2str(charge)];
            names = [names, string(details)];
        end
    end
end
for name = names
    load(strcat(folder,"far_field_data",name))

    fig = figure();%'units','normalized','outerposition',[0 0 1 1]);
    i=0;
    
    E_diff = zeros(1,20);
    ER_max = zeros(1,20);
    EL_max = zeros(1,20);

    i = i+1;

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
    % S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));  %% equivalent to -2*imag(Ex*conj(Ey))
    S3 = -2*imag(Ex.*conj(Ey));             %% equivalent to abs(Er)^2-abs(EL)^2
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
    S3_r = zeros(1,N);
    counter = zeros(1,N);
    for j = 1:length(R(:))
        r = R(j);
        idx = find((r <= r_grid(2:end)) .* (r > r_grid(1:end-1)));
        if ~isempty(idx)
            ER_r(idx) = ER_r(idx) + abs(ER(j)).^2;
            EL_r(idx) = EL_r(idx) + abs(EL(j)).^2;
            S3_r(idx) = S3_r(idx) + abs(S3(j));
            counter(idx) = counter(idx) + 1;
        end
        
    end
    ER_r(counter>0) = ER_r(counter>0) ./ counter(counter>0);%.*(2*pi*beta); %if one would want to integrate in the azimutal direction
    EL_r(counter>0) = EL_r(counter>0) ./ counter(counter>0);%.*(2*pi*beta); %instead of averaging only, it should moltyply by (2*pi*beta)
    S3_r(counter>0) = S3_r(counter>0) ./ counter(counter>0);%.*(2*pi*beta); %(counter>0) is required to avoid 0/0 division  
    % chi_r(counter>0) = chi_r(counter>0) ./ counter(counter>0);
    E_r = ER_r+EL_r;
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
    % subplot(6,4,i)
    plot(beta, S3_r, beta, E_r);
    legend("S3", "|E|^2 ")
    nicePlot
%     title(strcat(string(round(radial_FF * period)),'nmX',string(round(lateral_FF * period)),'nm'),'fontsize',12)
    % title(strcat("diameter ", string(round(start_radius*2/1000,1))))
    xlabel('u_T') 
    ylabel('Average azimutal Energy')
    nicePlot;
    drawnow;
   
    saveas(fig,strcat(folder,"radial_far_field_PLOT",name),'jpg')   
    saveas(fig,strcat(folder,"radial_far_field_PLOT",name),'fig')  
    ROI = beta<0.3;
    E_diff(i) = max(ER_r(ROI)) / (max(ER_r(ROI)) + max(EL_r(ROI))) * 100;
    ER_max(i) = max(ER_r(ROI));
    EL_max(i) = max(EL_r(ROI)); 
    close(fig)
end
%     end 

%%
% fig = figure
% subplot(1,2,1)
% d = linspace(1,4,20)*2;
% plot(d, E_diff, d, 100-E_diff)
% legend("RHC","LHC")
% title('Energy ripartition')
% xlabel('diameter [um]')
% ylabel('%')
% nicePlot
% subplot(1,2,2)
% plot(d, ER_max, d, EL_max)
% legend("RHC","LHC")
% xlabel('diameter [um]')
% ylabel('max(E) [a.u.]')
% nicePlot
% saveas(fig,strcat(folder,"Energy_ripartition",details,'.png')   )
% saveas(fig,strcat(folder,"Energy_ripartition",details,'.fig')   )


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