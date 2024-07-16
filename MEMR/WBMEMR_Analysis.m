%% WBMEMR_Analysis

% DPOAE swept Analysis
% Author: Samantha Hauser modified code from Hari Bharadwaj
% Created: August 2023
% Last Updated: August 23, 2023
% Purpose:
% Helpful info:

%% Import data
% Function to recursively search for files matching the pattern
cwd = '/Users/lizzyjensen/Desktop/Code/EBP/MEMR'

% Specify the main directory where you want to start searching
mainDirectory = '/Volumes/drive/MEMR Data from Sam/';

% Specify the pattern of filenames you are looking for
filenamePattern = 'WBMEMR_Full*.mat';  % Example: Looking for all .txt files

% Use dir to list all directories within the main directory
dirInfo = dir(mainDirectory);

% Initialize a cell array to store full paths of files found
fileList = {};

% Loop through each item in the directory
for i = 1 : length(dirInfo)
    % Skip '.' and '..' directories
    if dirInfo(i).isdir && ~strcmp(dirInfo(i).name, '.') && ~strcmp(dirInfo(i).name, '..')
        % Construct the full path to the current subdirectory
        subdirectory = fullfile(mainDirectory, dirInfo(i).name);
        
        % Use dir to list files in the current subdirectory that match the pattern
        files = dir(fullfile(subdirectory, filenamePattern));
        
        % Append the full paths of matching files to the fileList
        for j = 1 : length(files)
            fileList = [fileList; fullfile(subdirectory, files(j).name)];
        end
    end
end

% Display the list of files found
fprintf('Files found matching "%s" in directory "%s" and its subdirectories:\n', filenamePattern, mainDirectory);
for k = 1 : length(fileList)
    fprintf('%s\n', fileList{k});
end

% Example: Load .mat files identified in fileList

% Loop through each file in fileList
for i = 1:length(fileList)
    % Get the current file path
    matFilePath = fileList{i};
    
    % Load variables from .mat file into the workspace
    loadedData = load(matFilePath);
    
    % Access variables from the loaded .mat file
    % Example: Display the variables loaded
    disp(['Variables loaded from file: ', matFilePath]);
    disp(loadedData.data);  % Display all variables loaded from the file
    
    
    % Access specific variables if known
    % Example: Access a specific variable named 'data'
    if isfield(loadedData.data, 'stim')
        stim = loadedData.data.stim;  % Replace 'data' with your variable name
        disp(['Variable "stim" loaded from file: ', matFilePath]);
        disp(stim);  % Display the contents of the variable 'data'
        
        stim.resp = loadedData.data.resp.AllBuffs;
        subj = loadedData.data.info.subj.ID;
    else
        disp('No variable named "stim" found in file.');
    end
    
    % Further processing as required
    
    
    %% Analysis loop
    res = MEMRbyLevel_wCalib(stim);
    
    %% Additional calculations
    % smoothmem = true;
    % if smoothmem
    %     for k = 1:stim.nLevels
    %         res.MEMs(k, :) = sgolayfilt(res.MEM(k, :), 2, 35);
    %     end
    % end
    
    clinicfreq = [226, 678, 1000];
    for i = 1:length(clinicfreq)
        for j = 1:size(res.MEMs, 1)
            clinres(i,j) = interp1(res.freq, res.MEMs(j,:), clinicfreq(i), 'linear');
        end
    end
    
    power_WB = mean(abs(res.MEMs(:, res.ind)), 2);
    deltapow_WB = power_WB - mean(power_WB(1:2));
    deltapow_clin = clinres - mean(clinres(:,1:2),2);
    
    [val, fittedparams] = memfit(res.elicitor', deltapow_WB);
    growthfit = memgrowth(res.elicitor, fittedparams);
    res.threshold_WB = memgrowthinv(0.1, fittedparams)
    
    for i = 1:length(clinicfreq)
        [val, fittedparams] = memfit(res.elicitor', deltapow_clin(i,:)');
        growthfit = memgrowth(res.elicitor, fittedparams);
        res.threshold_clin(i) = memgrowthinv(0.1, fittedparams)
    end
    
    res.deltapow_clin = deltapow_clin
    res.deltapow_WB = deltapow_WB
    
    %% Plot
    figure(35);
    figure_prop_name = {'PaperPositionMode', 'units', 'Position'};
    figure_prop_val = {'auto', 'inches', [1 1 12 4]}; % xcor, ycor, xwid, yheight
    set(gcf,figure_prop_name,figure_prop_val);
%     if strcmp(datapath(1:3), '/Vo') | strcmp(datapath(1), 'D')
%         sgtitle(sprintf('%s | MEMR - WB | %s | %s', subj, condition, file(end-23:end-4)))
%     else
%         suptitle(sprintf('%s | MEMR - WB | %s | %s', subj, condition, file(end-23:end-4)))
%     end
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
    plot(res.elicitor, deltapow_WB, 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold_WB, [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow_WB+.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 14)
    
    drawnow;
    hold off;
    
    %% Data to save
    data.res = res;
    data.stim = stim;
    
    
    %% Export:

    cd('Results')
    
    fname = [subj,'_MEMR_Result'];
    
    print(gcf,[fname,'_figure'],'-dpng','-r300');
    save(fname,'data')
    cd(cwd);
    
end