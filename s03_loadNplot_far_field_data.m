clear,
close all
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
folder="SIM04_complex_outcouplers/far_field_data/";
top_charge=0;
n_g=5;
add_detail =["_exponential_design_TM"];%_Hdipole_moreTurns,"_TM_larger_domain","_TM_larger_domain_oring"]
    
for top_charge=0
    
%     details = strcat(add_detail,'_charge',num2str(top_charge),'_negative');
    details = strcat(add_detail,'_charge',num2str(top_charge),'_tilt45','_negative');
%     load(strcat(folder,"far_field_data_design_gd3_onSiO2_2_N35_positive"));
%     load(strcat(folder,"far_field_data","_positive_charge",string(top_charge),"_N",string(n_g)));
    load(strcat(folder,"far_field_data",details));
%     load(strcat(folder,"far_field_data","_charge",string(top_charge)));
    % load(strcat(folder,"far_field_dataRing_positive"));

    
    % convert to matlab reference frame
    Ex=transpose(Ex);
    Ey=transpose(Ey);
%     [ux,uy]=meshgrid(ux,uy);
%     ux=ux';
%     uy=uy';
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
    
    S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
    S0 = (abs(Ex).^2+abs(Ey).^2);
    chi = 0.5*asin( real(S3)./S0);
    
    E = sqrt(real(Ex).^2+real(Ey).^2)+1i*sqrt(imag(Ex).^2+imag(Ey).^2);

    % color maps
    col_vec=linspace(0,1,256);
    map_wave=[col_vec' col_vec' ones(256,1);
              ones(255,1) col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'];
    map_wave_dark=[zeros(255,1) zeros(255,1) col_vec(end-1:-1:1)';
                   col_vec' zeros(256,1) zeros(256,1)];
    map_yell_dark=[col_vec(end-1:-1:1)' col_vec(end-1:-1:1)' zeros(255,1) ;
                   col_vec' zeros(256,1) zeros(256,1)];
    map_oran_dark=[col_vec(end-1:-1:1)' col_vec(end-1:-1:1)'.^2 zeros(255,1) ;
                   col_vec' col_vec([1:end/2,end/2:-1:1])' zeros(256,1)];
    map_intensity=[ones(256,1) ones(256,1) col_vec(end:-1:1)';
                   ones(255,1) col_vec(end-1:-1:1)' zeros(255,1)];
    %%
    fig=figure('units','normalized','outerposition',[0 0 1 1]);
    subplot(2,3,1)
    plot_surf(ux,uy,abs(ER).^2,'hot',"Right circular polarization intensity");
    subplot(2,3,2)
    plot_surf(ux,uy,abs(EL).^2,'hot',"Left circular polarization intensity");
    subplot(2,3,4)
    plot_surf(ux,uy,angle(ER),'hot',"Right circular polarization phase");
    subplot(2,3,5)
    plot_surf(ux,uy,angle(EL),'hot',"Left circular polarization phase");
    subplot(2,3,3)
    plot_surf(ux,uy,abs(Ex).^2+abs(Ey).^2,'hot','Total Intensity');
%     plot_surf(ux,uy,real(S3./max(max(S0))),'jet',"S3/max(S0) Stokes parameter",1);
%     plot_surf(ux,uy,real(S3),map_yell_dark,"S3 Stokes parameter",'symmetric');
    subplot(2,3,6)
%     plot_surf(ux,uy,angle(E),'hot',"E=\surd(Ex^2+Ey^2) phase");
    plot_surf(ux,uy,real(tan(chi)),map_yell_dark,"eccentricity tan\chi",'symmetric',1);

    sgtitle({strcat('{\fontsize{8} ','design',strrep(details,'_','\_'),'}');...
        ['Topological charge ',num2str(top_charge)]},'fontsize',18,'fontweight','bold');
    
%     saveas(fig,strcat(folder,"far_field_PLOT",details),'png')    
end



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