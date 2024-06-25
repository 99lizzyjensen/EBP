%% MEMR by level helper function
function [res] = MEMRbyLevel(stim)

endsamps = ceil(stim.clickwin*stim.Fs*1e-3)+1; 

freq = linspace(200, 8000, 1024); 
MEMband = [500, 2000];
ind = (freq >= MEMband(1)) & (freq <= MEMband(2));

resp_freq = zeros(stim.nLevels, numel(freq)); 
bline_freq = zeros(stim.nLevels, numel(freq)); 
for k = 1:stim.nLevels
    fprintf(1, 'Analyzing level # %d / %d ...\n', k, stim.nLevels);
    temp = reshape(squeeze(stim.resp(k, :, 2:end, 1:endsamps)),...
        (stim.nreps-1)*stim.Averages, endsamps);
    tempf = pmtm(temp', 4, freq, stim.Fs)';
    resp_freq(k, :) = median(tempf, 1);
    
    blevs = k; % Which levels to use as baseline (consider 1:k)
    temp2 = squeeze(stim.resp(blevs, :, 1, 1:endsamps));
    
    if(numel(blevs) > 1)
        temp2 = reshape(temp2, size(temp2, 2)*numel(blevs), endsamps);
    end
    
    temp2f = pmtm(temp2', 4, freq, stim.Fs)';
    bline_freq(k, :) = median(temp2f, 1);
end

elicitor = 105 - stim.noiseatt;

MEM = pow2db(resp_freq ./ bline_freq);

res.freq = freq; 
res.MEM = MEM; 
res.elicitor = elicitor; 
res.ind = ind; 
res.nTrials = stim.Averages; 
end 