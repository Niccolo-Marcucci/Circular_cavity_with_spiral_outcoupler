clear,
% close all
folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
folder="SIM04_complex_outcouplers/far_field_data/";

top_charge=0;
n_g=5;
add_detail =["_metasurface_design_TM_filled"];%_TM_finerMesh_moreTurns_Hdipole,"_TM_larger_domain","_TM_larger_domain_oring"]
    
animated_quantity = 'E'; % options are E, RHC, LHC
for top_charge=0
    
    details = strcat(add_detail,'_charge',num2str(top_charge),'_negative');
    details = strcat(add_detail,'_charge',num2str(top_charge),'_tilt22','_negative');

    load(strcat(folder,"far_field_data",details));
    
    % convert to matlab reference frame (from Lumerical 2019)
    Ex=transpose(Ex);
    Ey=transpose(Ey);
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
%     Extmp=Ex;
%     Ex=Ex+Ey;
%     Ey=Ey-Extmp;
%     Ex=Ex+Ey*exp(1i*pi/2);
%     Ey=Ey-Extmp*exp(1i*pi/2);
%     Ex=Ey;
%     Ey=-Extmp;


    ER = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(-1i*pi/2);
    EL = sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(+1i*pi/2);

    ERx=    (Ex-1i*Ey)/2;
    ERy=+1i*(Ex-1i*Ey)/2;
    ELx=    (Ex+1i*Ey)/2;
    ELy=-1i*(Ex+1i*Ey)/2;
    
    %% color maps
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
    N = 50;
    dePhi = 2*pi/N;
    fig=figure;
    sp=20;
    title_name={strcat('{\fontsize{8} ','design',...
                       strrep(details,'_','\_'),'}');...
                strcat('\fontsize{18}\bfTopological charge ',...
                num2str(top_charge),", ",animated_quantity," intensity")};
    vecUX=ux(sp:sp:end);%,sp:sp:end);
    vecUY=uy(sp:sp:end);%,sp:sp:end);
    
    switch animated_quantity
        case "LHC"
            Ex=ELx;
            Ey=ELy;
        case "RHC"
            Ex=ERx;
            Ey=ERy;            
        case "E"
        otherwise
            error("cannot animate requested quantity")
    end
    E = sqrt(real(Ex).^2+real(Ey).^2)+1i*sqrt(imag(Ex).^2+imag(Ey).^2);
    massimo = max(max(max(real(E).^2)),max(max(imag(E).^2)));    
    for kk = 1:N*10
        E = sqrt(real(Ex).^2+real(Ey).^2)+1i*sqrt(imag(Ex).^2+imag(Ey).^2);
        figure(fig)
        hold off
        plot_surf(ux,uy,real(E).^2,'hot',title_name,'unitary',massimo)
        hold on
        vecEx=real(Ex(sp:sp:end,sp:sp:end));
        vecEy=real(Ey(sp:sp:end,sp:sp:end));
        quiver(vecUX,vecUY,vecEx,vecEy);
        Ex=Ex*exp(-1i*dePhi);
        Ey=Ey*exp(-1i*dePhi);
        nicePlot
        xlim([-0.2 0.2])
        ylim([-0.2 0.2])
        pause(0.05)
%         % Capture the plot as an image 
%         frame = getframe(fig); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         gifname=strcat(folder,"far_field_coloredAnimation",details,"_",...
%                                 animated_quantity,'.gif');
%         if kk == 1 
%           imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',0.05); 
%         else 
%           imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',0.05); 
%         end 
    end
%    saveas(figure(3),strcat(folder,"far_field_quiverPLOT",details),'png')
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