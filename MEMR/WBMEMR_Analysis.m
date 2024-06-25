%% WBMEMR_Analysis

% DPOAE swept Analysis
% Author: Samantha Hauser modified code from Hari Bharadwaj
% Created: August 2023
% Last Updated: August 23, 2023
% Purpose:
% Helpful info:

%% Import data
cwd = pwd;
if exist(datapath)
    cd(datapath)
else
    sprintf('No such datapath - Chin: %s, Cond: %s', subj, condition)
    return
end

datafile = dir(fullfile(cd,['WBMEMR*.mat']));
if length(datafile) < 1
    fprintf('No file for subject %s...Quitting!\n', subj);
    return
end

numOfFiles = size(datafile,1);
for fileNum = 1:numOfFiles
    cd(datapath)
    load(datafile(fileNum).name);
    file = datafile(fileNum).name;
    
    if strcmp(file(1), 'W')
        stim = data.stim;
        stim.resp = data.resp.AllBuffs; 
    end
   
    resp = stim.resp(:, 1:stim.Averages, :, :);
    
    cd(cwd);
    
    % %Call artifact rejection function
    % [stim] = artifact_rejection(stim);
    % %Compare GREATtrials to stim.rms-->NEW STIM.RMS
    % for L = 1:stim.nLevels
    %     for T = 1:stim.nTrials
    %         if stim.reject(L,T)
    %             stim.resp(L,T,:,:) = NaN;
    %         end
    %     end
    % end
    
    %% Analysis loop
    res = MEMRbyLevel_wCalib(stim);
    
    %% Additional calculations
    % smoothmem = true;
    % if smoothmem
    %     for k = 1:stim.nLevels
    %         res.MEMs(k, :) = sgolayfilt(res.MEM(k, :), 2, 35);
    %     end
    % end
    
    power = mean(abs(res.MEMs(:, res.ind)), 2);
    deltapow = power - mean(power(1:2));
    
    [val, fittedparams] = memfit(res.elicitor', deltapow);
    growthfit = memgrowth(res.elicitor, fittedparams);
    res.threshold = memgrowthinv(0.1, fittedparams)
    %% Plot
    figure(35);
    figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
    figure_prop_val = {'auto', 'inches', [1 1 12 4]}; % xcor, ycor, xwid, yheight
    set(gcf,figure_prop_name,figure_prop_val);
    if strcmp(datapath(1:3), '/Vo') | strcmp(datapath(1), 'D')
       sgtitle(sprintf('%s | MEMR - WB | %s | %s', subj, condition, file(end-23:end-4)))
    else
        suptitle(sprintf('%s | MEMR - WB | %s | %s', subj, condition, file(end-23:end-4)))
    end
    % Colorblind friendly continuous hue/sat changes
    cols = [103,0,31;
        178,24,43;
        214,96,77;
        244,165,130;
        253,219,199;
        247, 247, 247;
        209,229,240;
        146,197,222;
        67,147,195;
        33,102,172;
        5,48,97];
    cols = cols(end:-1:1, :)./255;
    
    % Plot
    subplot(1,3,1:2)
    semilogx(res.freq / 1e3, res.MEMs, 'linew', 2);
    xlim([0.2, 8]);
    xticks([0.25, 0.5, 1, 2, 4, 8])
    set(gca,'ColorOrder', cols)
    legend(num2str(res.elicitor'), 'location', 'eastoutside');
    xlabel('Frequency (kHz)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('\Delta Absorbed Power (dB)', 'FontSize', 14, 'FontWeight', 'bold');
    
    subplot(1,3,3)
    plot(res.elicitor, deltapow, 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold, [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow+.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 14)
    
    drawnow;
    hold off;
    
    %% Data to save
    data.res = res;
    data.stim = stim; 
    

    %% Export:
    
    cd(datapath);
    cd ..
    cd('Processed')
    if stim.fc == 7000
        fname = [subj,'_MEMR_HP_',condition,  file(end-24:end-4)];
    else
        fname = [subj,'_MEMR_WB_',condition,  file(end-24:end-4)];
    end
    
    print(gcf,[fname,'_figure'],'-dpng','-r300');
    save(fname,'data')
    cd(cwd);
    
    
end
