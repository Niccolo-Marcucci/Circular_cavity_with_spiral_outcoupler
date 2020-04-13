clear, close all
folder="SIM01_circular_cavity_without_outcoupler/spectra/";

for i=1:12   
    switch i
        case 1
            details = '_design_gd3_onSiO2_N50_positive';
        case 2
            details = '_design_gd3_N40_negative';
        case 3
            details = '_design_gd3_N40_positive';
        case 4
            details = '_design_gd3_onSiO2_N50_positive';
        case 5
            details = '_design_gd4_onSiO2_N40_positive';
        case 6
            details = '_design_Bad_N45_negative';
        case 7
            details = '_design_Bad9_N45_negative';
        case 8
            details = '_design8_N45_negative';
        case 9
            details = '_designDesc_N30_negative';
        case 10
            details = '_designDesc_N35_negative';
        case 11
            details = '_design3_N30_positive';
        case 12
            details = '_design3_N30_negative';
    end

    load(strcat(folder,"spectrum",details,".mat"))
    
    f=figure('outerposition',[0 0 1280 720]);
    subplot(1,2,1)
    plot(lambda_out*1e9,spectrum_out)
    xlim([530 610])
    xlabel('wavelength [nm]')
    ylabel('power density [AU] - log scale')
    title('log scale')
    set(gca,'yscale','log')
    nicePlot
    
    subplot(1,2,2)
    plot(lambda_in*1e9,spectrum_in)
    xlabel('wavelength [nm]')
    ylabel('power density [AU]')
    title('linear scale')
    set(gca,'xminorgrid','on');
    [pks,idxs] = findpeaks(spectrum_out);
    [~,pk_ix] = max(pks);
    idx = idxs(pk_ix);
    peak = lambda_in(idx)*1e9;
    xlim(round([peak-10 peak+10]))
	text(lambda_in(idx)*1e9+2,spectrum_in(idx)/5*4,...
        ['\lambda = ',num2str(peak),'nm',...
        newline,'\Delta\lambda = ',num2str(peak/quality),'nm',...
        newline,'Q = ',num2str(round(quality))],'linewidth',2,...
        "fontsize",18,'backgroundcolor','w','edgecolor','r');
    nicePlot
    sgtitle(['Spacer = ',num2str(D*1e9),'nm;  DBR period =',...
             num2str(period_DBR*1e9),'nm;  N_{rings}=',num2str(N)],...
             "fontsize",18,'fontweight','bold');
    saveas(f,strcat(folder,"spectrum",details),'fig')
    close(f)
end