clear,
% close all
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
folder="SIM04_complex_outcouplers/far_field_data/";%";%
n_g=5;
top_charge=0;
names = ["_design_gd3_onSiO2_positive_filled_Ishape_Dphi-60_N12_sigma1_charge2"];
% 
% i = 0;
% for  DD_phi = [30 45 60 90]
%     i = i+1;
%     names(i) = strcat("_TM_Descrovi_multilayered_outcoupler_negative_filled_Dphi60_N12_sigma-1_DDphi", string(DD_phi));
% end
% % 
% 	details =strcat('_design_gd3_onSiO2_positive_filled_Dphi60_N21_sigma1_m1_radialFF',string(radial_FF),'_lateralFF',string(lateral_FF), '.mat');
%     details = "_TM_Descrovi_negative_filled_Dphi60_N12_sigma-1_scwidth24";


% load(strcat(folder,"far_field_data","_TM_Descrovi_empty"));

% convert to matlab reference frame
% Ex=transpose(Ex);
% Ey=transpose(Ey);
% 
% ER0 = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(-1i*pi/2);
% EL0 = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);

for name = names
    

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
    E2 = abs(Ex).^2+abs(Ey).^2;
    ER = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(-1i*pi/2);
    EL = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);
    E_min = min( min(min(abs(ER).^2)), min(min(abs(EL).^2)));
    E_max = max( max(max(abs(ER).^2)), max(max(abs(EL).^2)));
    S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
    S0 = (abs(Ex).^2+abs(Ey).^2);
    chi = 0.5*asin( real(S3)./S0);
    
    E = sqrt(real(Ex).^2+real(Ey).^2)+1i*sqrt(imag(Ex).^2+imag(Ey).^2);

    % color maps
    col_vec=linspace(0,1,256);
    map_wave=[col_vec' ones(256,1) col_vec' ;
              ones(255,1) col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'];
    map_wave_dark=[zeros(255,1) zeros(255,1) col_vec(end-1:-1:1)';
                   col_vec' zeros(256,1) zeros(256,1)];
    map_yell_dark=[col_vec(end-1:-1:1)' col_vec(end-1:-1:1)' zeros(255,1) ;
                   col_vec' zeros(256,1) zeros(256,1)];
    map_yell_lght=[ones(256,1) ones(256,1) col_vec';
                   ones(255,1) col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'];
    map_oran_dark=[col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'.^2 zeros(255,1) ;
                   col_vec' col_vec([1:end/2,end/2:-1:1])' zeros(256,1)];
    map_intensity=[ones(256,1) ones(256,1) col_vec(end:-1:1)';
                   ones(255,1) col_vec(end-1:-1:1)' zeros(255,1)];
    %%
    xy_lim = [-1 1]*.2;
%     fig = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,2,1)
    plot_surf(ux,uy,abs(ER).^2,'hot',"Right circular polarization intensity");
    caxis([E_min E_max])
    xlim(xy_lim)
    ylim(xy_lim)
    subplot(2,2,2)
    plot_surf(ux,uy,abs(EL).^2,'hot',"Left circular polarization intensity");
    caxis([E_min E_max])
    xlim(xy_lim)
    ylim(xy_lim)
    subplot(2,2,3)
    plot_surf(ux,uy,angle(ER),'hot',"Right circular polarization phase");
    xlim(xy_lim)
    ylim(xy_lim)
    subplot(2,2,4)
    plot_surf(ux,uy,angle(EL),'hot',"Left circular polarization phase");
    xlim(xy_lim)
    ylim(xy_lim)
    subplot(2,3,3)
    plot_surf(ux,uy,E2,'parula','Total Intensity');
%     plot_surf(ux,uy,real(S3./max(max(S0))),'jet',"S3/max(S0) Stokes parameter",1);
%     plot_surf(ux,uy,real(S3),map_yell_dark,"S3 Stokes parameter",'symmetric');
    xlim(xy_lim)
    ylim(xy_lim)
    subplot(2,3,6)
%     plot_surf(ux,uy,angle(E),'hot',"E=\surd(Ex^2+Ey^2) phase");
    
    plot_masked(ux,uy,real(tan(chi)),E2,map_wave,"eccentricity tan\chi",'symmetric',1);
%     xlim(xy_lim)
%     ylim(xy_lim)
    
%     sgtitle({strcat('{\fontsize{8} ',strrep(details,'_','\_'),'}');...
%         ['Topological charge ',num2str(top_charge)]},'fontsize',18,'fontweight','bold');
    
    saveas(fig,strcat(folder,"far_field_PLOT",name),'png')   
end
% 
function rgbImage = getRGB(data,map)
    dmin = min(data(:));
    dmax = max(data(:));
    data = round((data - dmin)/(dmax - dmin)*( length(map)-1 ) ) +1 ;
    rgbImage = ind2rgb(data,map);
end

function plot_masked(ux,uy,quantity,mask,map,picture_title,zsymmetry,massimo)
    rgbImage = getRGB(quantity,map);
    mask = (mask - min(mask(:))) / (max(mask(:)) - min(mask(:)));
    binary_mask = repmat(mask < 0.05,1,1,3);
    
    grayImage =  repmat(mean(rgbImage,3),1,1,3);
    rgbImage(binary_mask) = grayImage(binary_mask);
    
    imshow(rgbImage);
    
    ax = gca;
    set(ax,'YDir','normal') 
    if nargin > 4
        colormap(ax,map)
    if nargin > 5
        title(picture_title)
    if nargin > 6 
        if zsymmetry == "symmetric"
            c = max(abs([min(min(quantity)),max(max(quantity))]));
            caxis([-c c])
        elseif zsymmetry == "unitary"
            caxis([0 1])
        else
            error('Wrong "symmetry" signment')
        end
    if nargin > 7
        if zsymmetry == "symmetric"
            caxis([-1 1]*massimo)
        elseif zsymmetry == "unitary"
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

function plot_surf(ux,uy,quantity,map,picture_title,zsymmetry,massimo)
%     s=surface(ux,uy,quantity);
%     s.EdgeColor='none';
    imagesc(ux,uy,quantity);
    ax = gca;
    set(ax,'YDir','normal') 
    colorbar
    xlabel("ux");
    ylabel('uy');
    axis('square')
    if nargin > 3
        colormap(ax,map)
    if nargin > 4
        title(picture_title)
    if nargin > 5 
        if zsymmetry == "symmetric"
            c = max(abs([min(min(quantity)),max(max(quantity))]));
            caxis([-c c])
        elseif zsymmetry == "unitary"
            caxis([0 1])
        else
            error('Wrong "symmetry" signment')
        end
    if nargin > 6
        if zsymmetry == "symmetric"
            caxis([-1 1]*massimo)
        elseif zsymmetry == "unitary"
            caxis([0 1]*massimo)
        else
            error('Wrong "symmetry" signment')
        end
    end;end;end;end
end