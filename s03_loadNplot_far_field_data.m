clear,
% close all
clear all
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM04_complex_outcouplers/far_field_data/";
folder="SIM05_metasurface_outcoupler/a/far_field_data/";%scatterTests_Gold_topped_negative/far_field_data/";%";%

names = [];
for dphi = 60%-60:120:60
    % for sigma = 1%-1:2:1
        % for charge = 0%-1:1
    i = 0;
    sigma = -1;
    charge = -2;
    for sc_width  = [75] % [25, 50, 75, 100, 125, 150]
        for sc_length = [250] %[250, 275, 300, 325]
            details = ['_TM_SiO2TiO2_532_N9negative_GoldPallik_filled_scShapeI_Dphi',num2str(dphi),'_N12_sigma',num2str(sigma),'_charge', num2str(charge), '_scWidth', num2str(sc_width), '_scLength', num2str(sc_length)];
            % details = ['_TM_SiO2TiO2_532_N9positive_filled_scShapeI_Dphi',num2str(dphi),'_N12_sigma',num2str(sigma),'_charge', num2str(charge)];
            % details = ['_negative_charge', num2str(charge)];
            % details = ['_design_gd3_onSiO2_positive_filled_Ishape_Dphi-60_N12_sigma1_charge2'];
            names = [names, string(details)];
        end
    end
end
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
    load(strcat(folder,"far_field_data",name))
    % load(strcat(folder,"near_field_index_data",name))
    % fig1=figure();
    % imagesc(x,y,real(transpose(index)));
    % ax = gca;
    % set(ax,'YDir','normal') 
    % colorbar
    % xlabel("x");
    % ylabel('y');
    % axis('equal')
    % saveas(fig1,strcat(folder,"grating_PLOT",name),'png');
    % close(fig1);

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
    % S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));  %% equivalent to -2*imag(Ex*conj(Ey))
    S3 = -2*imag(Ex.*conj(Ey));             %% equivalent to abs(Er)^2-abs(EL)^2
    S0 = (abs(Ex).^2+abs(Ey).^2);
    chi = 0.5*asin( real(S3)./S0);
    
    E = sqrt(real(Ex).^2+real(Ey).^2)+1i*sqrt(imag(Ex).^2+imag(Ey).^2);

    % color maps
    col_vec=linspace(0,1,256);
    map_wave=[col_vec' ones(256,1) col_vec'  ;
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
    map_wave_yellow = hsv2rgb([linspace(1/3,0,256*2)' ones(256*2,1) ones(256*2,1)]);
    %%
    xy_lim = .2;
    ER = ER(abs(ux)<xy_lim,abs(uy)<xy_lim);
    EL = EL(abs(ux)<xy_lim,abs(uy)<xy_lim);
    E2 = E2(abs(ux)<xy_lim,abs(uy)<xy_lim);
    S0 = S0(abs(ux)<xy_lim,abs(uy)<xy_lim);
    S3 = S3(abs(ux)<xy_lim,abs(uy)<xy_lim);
    chi= chi(abs(ux)<xy_lim,abs(uy)<xy_lim);
    ux = ux(abs(ux)<xy_lim);
    uy = uy(abs(uy)<xy_lim);
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,3,1)
    ax = plot_surf(ux,uy,abs(ER).^2,'hot',"Right circular polarization intensity");
    caxis([E_min E_max])
    % xlim(xy_lim)
    % ylim(xy_lim)
    subplot(2,3,2)
    plot_surf(ux,uy,abs(EL).^2,'hot',"Left circular polarization intensity");
    caxis([E_min E_max])
    % xlim(xy_lim)
    % ylim(xy_lim)
    subplot(2,3,4)
    plot_surf(ux,uy,angle(ER),'hot',"Right circular polarization phase");
    % plot_masked(ux,uy,angle(ER),E2,hot,"S3/S0 & Total Intensity",'symmetric', pi , 1)
    % xlim(xy_lim)
    % ylim(xy_lim)
    subplot(2,3,5)
    plot_surf(ux,uy,angle(EL),'hot',"Left circular polarization phase");
    % xlim(xy_lim)
    % ylim(xy_lim)
    subplot(2,3,3)
    plot_surf(x,y,real(transpose(index)),'cool','Refractive index map');
    xlim([-4,4]*1e-6)
    ylim([-4,4]*1e-6)
    xlabel("x");
    ylabel('y');
    % plot_surf(ux,uy,E2,'parula','Total Intensity');
%     plot_surf(ux,uy,real(S3./max(max(S0))),'jet',"S3/max(S0) Stokes parameter",1);
%     plot_surf(ux,uy,real(S3),map_yell_dark,"S3 Stokes parameter",'symmetric');
    % xlim(xy_lim)
    % ylim(xy_lim)
    subplot(2,3,6)
%     plot_surf(ux,uy,angle(E),'hot',"E=\surd(Ex^2+Ey^2) phase");
    
    % plot_masked(ux,uy,real(tan(chi)),E2,map_wave,"eccentricity tan\chi",'symmetric',1, 0.1);
    ax2 = plot_masked(ux,uy,S3./S0,E2,map_wave_yellow,"far field",'symmetric', 1 , 1);
    set(ax2,'XTickLabel', ax.XTickLabel, 'XTick', linspace(1, length(E2(1,:)),length(ax.XTickLabel)),...
            'YTickLabel', ax.YTickLabel, 'YTick', linspace(1, length(E2(:,1)),length(ax.YTickLabel)));
    % plot_surf(ux,uy,S3,map_wave,"S3",'symmetric');
    % xlim(xy_lim)
    % ylim(xy_lim)
    
%     sgtitle({strcat('{\fontsize{8} ',strrep(details,'_','\_'),'}');...
%         ['Topological charge ',num2str(top_charge)]},'fontsize',18,'fontweight','bold');
    
    % saveas(fig,strcat(folder,"far_field_PLOT",name),'png')
    % close(fig);
end
% 
function rgbImage = getRGB(data,map,range)
    dmin = min(data(:));
    dmax = max(data(:));
    if nargin > 2 
        data =   round( (data - min(range)) / (max(range)-min(range)) *( length(map)-1 ) ) +1 ;  % normalize again inside range and scale to map
    else
        data = round((data - dmin)/(dmax - dmin)*( length(map)-1 ) ) +1 ;                           % normalize between 0 and 1 and scale to map
    end
    rgbImage = ind2rgb(data,map);
end

function ax = plot_masked(ux,uy,quantity,mask,map,picture_title,zsymmetry,massimo,gamma)
    if nargin > 6
        if zsymmetry == "symmetric"
            range = [-1,1] * max(abs(quantity(:)));
            rgbImage = getRGB(quantity,map,range);            
        elseif zsymmetry == "unitary"
            range = [0,1] ;
            rgbImage = getRGB(quantity,map,range); 
        end
        if nargin > 7
            if zsymmetry == "symmetric"
                range = [-1 1]*massimo;
                rgbImage = getRGB(quantity,map,range); 
            elseif zsymmetry == "unitary"
                range = [0 1]*massimo;
                rgbImage = getRGB(quantity,map,range); 
            end
        end
    else
        rgbImage = getRGB(quantity,map);
    end

    minimo = min(mask(:));

    minimo = 0; % force minimo to 0 so higlight the presence of a baseline
    massimo= max(mask(:));
    mask = (mask - minimo) / (massimo - minimo); % normalize mask
    % binary_mask = repmat(mask < threshold,1,1,3);
    
    % grayImage =  repmat(mean(rgbImage,3),1,1,3);
    % grayImage = grayImage * 0.8;
    % rgbImage(binary_mask) = grayImage(binary_mask);
    
    hsvImage = rgb2hsv(rgbImage);
    satImage = hsvImage(:,:,2);
    satImage(mask>2/3)= ((1+2/3*3/2)-3/2*mask(mask>2/3)).^gamma;
    valImage = hsvImage(:,:,3);
    valImage(mask<2/3)= (3/2*mask(mask<2/3)).^gamma;
    hsvImage(:,:,2) = satImage;
    hsvImage(:,:,3) = valImage;
    % hsvImage(:,:,3) = mask.^gamma;
    rgbImage_new = hsv2rgb(hsvImage);

    imshow(rgbImage_new);
    
    ax = gca;
    % set(ax, 'XTickLabel', string(linspace(-0.2,0.2,7)), 'XTick', linspace(1, length(image_to_show(1,:)),7),...
    %         'YTickLabel', string(linspace(-0.2,0.2,7)), 'YTick', linspace(1, length(image_to_show(:,1)),7));
    set(ax,'YDir','normal','tickdir','both', 'visible','on') 
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
                end
            end
        end
    end
    % colorbar
    xlabel("ux");
    ylabel('uy');
    axis('square')
    
    %% colobar2D
    p = ax.Position';
    rel_size = 0.4;
    ax_map = axes('Position',[p(1)+p(3)*(1-rel_size*3/4), p(2)+(p(3)+p(4))/2, p(3)*rel_size, p(4)*rel_size ]);
    x_map = 1:length(map);
    y_map = linspace(0,1,256/2*3);

    [X_MAP,Y_MAP] = meshgrid(x_map,y_map);
    hsv_map = rgb2hsv(ind2rgb(X_MAP,map));      % converti la lista di indici in rgb e quello che risulta in hsv
    sat_map = hsv_map(:,:,2);
    sat_map(Y_MAP>2/3) = ((1+2/3*3/2)-3/2*Y_MAP(Y_MAP>2/3)).^gamma;
    val_map = hsv_map(:,:,3);
    val_map(Y_MAP<2/3)= (3/2*Y_MAP(Y_MAP<2/3)).^gamma;
    hsv_map(:,:,2) = sat_map;
    hsv_map(:,:,3) = val_map;
    % hsv_map(:,:,3) = Y_MAP;
    imshow(hsv2rgb(hsv_map))
    ax_map.Position = [p(1)+p(3)*(1-rel_size*3/4), p(2)+p(3)*1.1, p(3)*rel_size, p(4)*rel_size ];
    ax_map.XLabel.String = "S3/S0";
    ax_map.XAxisLocation = 'top';
    ax_map.XTick = linspace(1, length(x_map),3);
    ax_map.XTickLabel = [{"LHC"}, {"Linear"}, {"RHC"}];

    ax_map.YDir = 'normal';
    ax_map.YLabel.String = "S0";
    ax_map.YLabel.Rotation = 0;
    ax_map.YAxisLocation = 'right';
    ax_map.YTick = linspace(1, length(y_map),5);
    ax_map.YTickLabel = compose("%1.2e",linspace(minimo, massimo,5));


    set(ax_map,'visible','on','tickdir','both',...
               'position',[p(1)+p(3)*(1-rel_size*3/4), p(2)+p(4)-(p(4))/4, ax_map.Position(3), ax_map.Position(4) ]) 
end

function ax = plot_surf(ux,uy,quantity,map,picture_title,zsymmetry,massimo)
%     s=surface(ux,uy,quantity);
%     s.EdgeColor='none';
    imagesc(ux,uy,quantity);
    ax = gca;
    set(ax,'YDir','normal','tickdir','both') 
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