%%Clinical Threshold Results
%Evaluates the WBMEMR at WB, 226, 678, and 1000 Hz probe tone
%Made by LJ
%Created 7/3/24

%% Import Data
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

    %% Pull Variables from WBMEMR_Analysis
    res = MEMRbyLevel_wCalib(stim);
    
    clinicfreq = [226, 678, 1000];
    for i = 1:length(clinicfreq)
        for j = 1:size(res.MEMs, 1)
            clinres(i,j) = interp1(res.freq, abs(res.MEMs(j,:)), clinicfreq(i), 'linear');
        end
    end
    
    %individual max
    %find biggest peak for each elicitor level in res.MEMs
    for i = 1:size(res.MEMs, 1)
        [peak, index] = max(abs(res.MEMs(i, :))); %how to limit : to be from 100-4000 Hz?
        peaks(i) = peak
        peakfreq(i) = res.freq(index)
    end
    
    deltapow_peak = peaks - mean(peaks(1:2));
    
 
    [val, fittedparams] = memfit(res.elicitor', deltapow_peak');
    growthfit = memgrowth(res.elicitor, fittedparams);
    res.threshold_peak = memgrowthinv(0.1, fittedparams)
    
    res.peakfreq = peakfreq
 
    
    % do all the stuff below again (deltapow and lower) to create
    % -- threshold, growth curve 
    % also make sure you save all the freq that were the peaks per each
    % elicitor
    
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
    figure; 
    subplot(2,3,1) %1 row, 4 columns, first position
    plot(res.elicitor, deltapow_clin(1, :), 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold_clin(1,1), [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    title('226 Hz Probe Tone')
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow_clin(1, :) +.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 7)
    
    subplot(2,3,2) %1 row, 4 columns, first position
    plot(res.elicitor, deltapow_clin(2, :), 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold_clin(1,2), [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    title('678 Hz Probe Tone')
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow_clin(2, :)+.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 7)
    
    subplot(2,3,3) %1 row, 4 columns, first position
    plot(res.elicitor, deltapow_clin(3, :), 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold_clin(1,3), [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    title('1000 Hz Probe Tone')
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow_clin(3, :)+.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 7)
    
    subplot(2,3,4)
    plot(res.elicitor, deltapow_WB, 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold_WB, [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    title('WB Probe Tone')
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow_WB+.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 7)
    
    subplot(2,3,5)
    plot(res.elicitor, deltapow_peak, 'ok-', 'linew', 2);
    hold on;
    plot(res.threshold_peak, [0.1], 'or', 'MarkerSize', 7)
    plot(res.elicitor, growthfit, 'k--', 'linew', 2);
    title('Peak of Each Level')
    xticks([45:10:105])
    xlabel('Elicitor Level (dB FPL)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel({'\Delta Absorbed Power (dB)', '0.5-2 kHz Average'}, 'FontSize', 14, 'FontWeight', 'bold');
    ymax = max(deltapow_peak+.05);
    ylim([0,ymax])
    xlim([40,110])
    set(gca, 'FontSize', 7)
    
    hold off;
    
    %% %% Data to save
    data.res = res;
    data.stim = stim;
    
    
    %% Export:
    cd(cwd)
    cd('Results')
    
    fname = [subj,'_MEMR_Result'];
    
    print(gcf,[fname,'_figure'],'-dpng','-r300');
    save(fname,'data')
    cd(cwd);
    
end
    
