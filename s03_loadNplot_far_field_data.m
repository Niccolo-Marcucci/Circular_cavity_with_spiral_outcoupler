clear,
close all
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge_and_ngrating/far_field_data/";
% folder="SIM03_circular_cavity_spiral_outcoupler/far_field_data/";
% folder="SIM02_no_cavity_spiral_outcoupler/sweep_charge/far_field_data/";
folder="SIM02_no_cavity_spiral_outcoupler/far_field_data/";
top_charge=0;
n_g=5;
for add_detail =["_TM_finerMesh"]%,"_TM_larger_domain","_TM_larger_domain_oring"]
    
    
details = strcat(add_detail,'_charge',num2str(top_charge),'_negative');
%     load(strcat(folder,"far_field_data_design_gd3_onSiO2_2_N35_positive"));
%     load(strcat(folder,"far_field_data","_positive_charge",string(top_charge),"_N",string(n_g)));
    load(strcat(folder,"far_field_data",details));
%     load(strcat(folder,"far_field_data","_charge",string(top_charge)));
    % load(strcat(folder,"far_field_dataRing_positive"));
    %
    % useful_ux = abs(ux)<0.2;
    % useful_uy = abs(uy)<0.2;
    
    [ux,uy]=meshgrid(ux,uy);
    ux=ux';
    uy=uy';
    % ux = ux(useful_ux,useful_uy);
    % uy = uy(useful_ux,useful_uy);
    % E_phi = E_phi(useful_ux,useful_uy);
    % E_theta = E_theta(useful_ux,useful_uy);
    
    % since cos(theta) = uz
%     theta = real( acos( sqrt(1 - ux.^2 - uy.^2)));
%     cos_phi = ux./sin(theta);
%     sin_phi = uy./sin(theta);
%     
%     % compute Ex Ey from Etheta and Ephi
%     Ex = E_theta.*cos_phi- E_phi.*sin_phi;
%     Ey = E_theta.*sin_phi+ E_phi.*cos_phi;
    
    % add Ez if it not negligible
%     Ex = Ex + Ez./cos_phi;
%     Ey = Ey + Ez./sin_phi;

    EL = +sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(-1i*pi/2);
    ER = -sqrt(2)/2*Ex + sqrt(2)/2*Ey*exp(-1i*pi/2);
    
    S3 = 1i*(Ex.*conj(Ey)-Ey.*conj(Ex));
    S0 = (abs(Ex).^2+abs(Ey).^2);
    chi = 0.5*asin( real(S3)./S0);
    
    E = sqrt(Ex.^2+Ey.^2);
%     figure
%     plot_surf(ux,uy,sin_phi,'hot',"something about e_phi");
%     plot_surf(ux,uy,theta,'hot',"something about e_theta");
    
    figure
%     subplot(1,2,1)
    plot_surf(ux,uy,abs(Ey).^2   +abs(Ex).^2,'hot','Intensity from Ex and Ey');
%     subplot(1,2,2)
    % plot_surf(ux,uy,abs(E_phi).^2+abs(E_theta).^2,'hot','Intensity from E_\phi and E_\theta');
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
%     plot_surf(ux,uy,real(S3./max(max(S0))),'jet',"S3/max(S0) Stokes parameter",1);
    plot_surf(ux,uy,real(S3),'jet',"S3 Stokes parameter",1);
    subplot(2,3,6)
%     plot_surf(ux,uy,angle(E),'hot',"E=\surd(Ex^2+Ey^2) phase");
    plot_surf(ux,uy,abs(tan(chi)),'hot',"eccentricity |tan\chi|");
    
    sgtitle({strcat('{\fontsize{8} ','design',strrep(details,'_','\_'),'}');...
        ['Topological charge ',num2str(top_charge)]},'fontsize',18,'fontweight','bold');
%     saveas(figure(2),strcat(folder,"far_field_PLOT",details),'png')
    
    N = 100;
    phi = 2*pi/N;
    fig=figure   ; 
    sp=20;
    E2=abs(sqrt(Ex.^2+Ey.^2)).^2;
    h1=contour(ux,uy,E2);
    title({strcat('{\fontsize{8} ','design',strrep(details,'_','\_'),'}');...
        ['Topological charge ',num2str(top_charge)]},'fontsize',18,'fontweight','bold');
    axis('square')
    hold on
    usedUX=ux(sp:sp:end,sp:sp:end);
    usedUY=uy(sp:sp:end,sp:sp:end);
    usedEx=real(Ex(sp:sp:end,sp:sp:end));
    usedEy=real(Ey(sp:sp:end,sp:sp:end));
    h2=quiver(usedUX,usedUY,usedEx,usedEy);
    nicePlot
    for kk = 1:N
        Ex=Ex*exp(-1i*phi);
        Ey=Ey*exp(-1i*phi);
        usedEx=real(Ex(sp:sp:end,sp:sp:end));
        usedEy=real(Ey(sp:sp:end,sp:sp:end));
        h2.UData = usedEx;
        h2.VData = usedEy;
        xlim([-0.2 0.2])
        ylim([-0.2 0.2])
        pause(0.05)
%         % Capture the plot as an image 
%         frame = getframe(fig); 
%         im = frame2im(frame); 
%         [imind,cm] = rgb2ind(im,256); 
%         % Write to the GIF File 
%         gifname=strcat(folder,"far_field_quiverAnimation",details,'.gif');
%         if kk == 1 
%           imwrite(imind,cm,gifname,'gif', 'Loopcount',inf,'DelayTime',0.05); 
%         else 
%           imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',0.05); 
%         end 
    end
%    saveas(figure(3),strcat(folder,"far_field_quiverPLOT",details),'png')
end



function plot_surf(ux,uy,quantity,map,picture_title,symmetric)
%     s=surf(ux,uy,Quantity);
%     s.EdgeColor='none';
%     view(2)
    imagesc(ux(:,1),uy(1,:),quantity);
    ax=gca;
    if nargin > 3
        colormap(ax,map)
    if nargin > 4
        title(picture_title)
    if nargin > 5 && symmetric
        c = max(abs([min(min(quantity)),max(max(quantity))]));
        caxis([-c c])
    end;end;end
    colorbar
    xlabel("ux");
    ylabel('uy');
    axis('square')
end