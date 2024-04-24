clear,
% close all
clear all
folder="SIM05_metasurface_outcoupler/scatterTests_PMMA_topped_negative/far_field_data/";%";%

power_ratio = [];
total_power = [];
max_chi = [];
chi_at_max = [];
sc_width_v  = [25, 50, 75, 100, 125, 150];
sc_length_v = [250, 275, 300, 325];
for dphi = -60%-60:120:60
    % for sigma = 1%-1:2:1
    %     for charge = 0%-1:1
    k = 0;
    sigma = 1;
    charge = 0;
    for sc_width  = sc_width_v
        k = k+1;
        l = 0;
        for sc_length = sc_length_v
            l=l+1;
            details = ['_TM_AlOTiO2_N10negative_filled_scShapeI_Dphi',num2str(dphi),'_N12_sigma',num2str(sigma),'_charge', num2str(charge), '_scWidth', num2str(sc_width), '_scLength', num2str(sc_length)];
            name =  string(details);
    load(strcat(folder,"far_field_data",name))

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

    R = sqrt(Ux.^2 +Uy.^2);
    THETA = atan2(Uy,Ux);
    
    N = 200;
    r_grid = linspace(0,1,N+1); % max(max(R)),N+1);
    theta_grid = linspace(0,2*pi,N+1);
    beta = r_grid(1:end-1);
    ER_r = zeros(1,N);
    EL_r = zeros(1,N);
    chi_r = zeros(1,N);
    S3_r = zeros(1,N);
    counter = zeros(1,N);
    for j = 1:length(R(:))
        r = R(j);
        idx = find((r <= r_grid(2:end)) .* (r > r_grid(1:end-1)));  % inserts R(j) into an regularly spaced grid
        if ~isempty(idx)
            ER_r(idx) = ER_r(idx) + abs(ER(j)).^2;
            EL_r(idx) = EL_r(idx) + abs(EL(j)).^2;
            chi_r(idx) = chi_r(idx) + abs(chi(j));
            S3_r(idx) = S3_r(idx) + abs(S3(j));
            counter(idx) = counter(idx) + 1;
        end
        
    end
    ER_r(counter>0) = ER_r(counter>0) ./ counter(counter>0);%.*(2*pi*beta); %if one would want to integrate in the azimutal direction
    EL_r(counter>0) = EL_r(counter>0) ./ counter(counter>0);%.*(2*pi*beta); %instead of averaging only, it should moltyply by (2*pi*beta)
    S3_r(counter>0) = S3_r(counter>0) ./ counter(counter>0);%.*(2*pi*beta); %(counter>0) is required to avoid 0/0 division  
    chi_r(counter>0)= chi_r(counter>0)./ counter(counter>0);
    E_r = ER_r+EL_r;

    % plot(beta, S3_r, beta, E_r)
    % [~, idx_peak] = max(ER_r<0.13);
    % chi_at_max (k,l) = chi_r(idx_peak);
    max_chi(k,l) = max(chi_r(E_r>E_max*0.1));% & beta<0.13));
    [~, idx_peak] = max(abs(S3_r));
    max_S3(k,l) = S3_r(idx_peak);
    [~, idx_peak] = max(E_r);
    S3_at_max(k,l) = S3_r(idx_peak);
    power_ratio(k,l) = sum(ER_r)/sum(EL_r);
    total_power(k,l) = sum(E_r); % if ER and ER have been averaged and not 
                                 % integrated, this power is not an actual 
                                 % power, but still gives an estimate of
                                 % the total power useful to us
                                 
    disp([k,l])
end;end;end
%%
y = [sc_width_v(1) sc_width_v(end)]; %- (sc_width_v(2) -sc_width_v(1))/2;
x = [sc_length_v(1) sc_length_v(end)];%- (sc_length_v(2) -sc_length_v(1))/2;
figure
subplot(2,2,1)
imagesc(x,y,S3_at_max);
% s.EdgeColor='none';
title('S3 at max')
xlabel('scatter width - nm')
ylabel('scatter length - nm')
colorbar
axs(1) = gca;

subplot(2,2,2)
imagesc(x,y,max_S3);
% s.EdgeColor='none';
title('max S3' )
xlabel('scatter width - nm')
ylabel('scatter length - nm')
colorbar
axs(2) = gca;

subplot(2,2,3)
imagesc(x,y,power_ratio);
title('ER/EL')
xlabel('scatter width - nm')
ylabel('scatter length - nm')
colorbar
axs(3) = gca;

subplot(2,2,4)
imagesc(x,y,total_power);
title('total power')
xlabel('scatter width - nm')
ylabel('scatter length - nm')
colorbar
axs(4) = gca;

for ax = axs
    ax.YAxis.TickValues = sc_width_v;
    ax.YAxis.TickLabels = string(sc_width_v);
    ax.YAxis.TickDirection = 'out';
    ax.XAxis.TickValues = sc_length_v;
    ax.XAxis.TickLabels = string(sc_length_v);
    ax.XAxis.TickDirection = 'out';
end

function rgbImage = getRGB(data,map)
    dmin = min(data(:));
    dmax = max(data(:));
    data = round((data - dmin)/(dmax - dmin)*( length(map)-1 ) ) +1 ;
    rgbImage = ind2rgb(data,map);
end

function plot_masked(ux,uy,quantity,mask,map,picture_title,zsymmetry,massimo,threshold)
    rgbImage = getRGB(quantity,map);
    mask = (mask - min(mask(:))) / (max(mask(:)) - min(mask(:)));
    binary_mask = repmat(mask < threshold,1,1,3);
    
    grayImage =  repmat(mean(rgbImage,3),1,1,3);
    grayImage = grayImage * 0.8;
    rgbImage(binary_mask) = grayImage(binary_mask);
    
    imshow(rgbImage(abs(ux)<0.2,abs(uy)<0.2,:));
    
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