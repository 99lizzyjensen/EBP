%% MEMR by level helper function
function [res] = MEMRbyLevel_wCalib(stim)

load('Thevenin_ER10X.mat');
% Ps = calib.Ps_1; 
% Zs = calib.Zs_1;
% att = calib.Attenuation; 
% f = calib.freq; 
% vo = calib.vo; 

fcalib = f;
endsamps = ceil(stim.clickwin*stim.Fs*1e-3)+1; 

freq = linspace(200, 8000, 1024); 
MEMband = [500, 2000];
ind = (freq >= MEMband(1)) & (freq <= MEMband(2));

%---------
% Additional info to match calibEar.m
delay = (round(endsamps / 3) + 97) / stim.Fs;
mat2volts = 5.0;
Vo = rfft(vo)*mat2volts*db2mag(-1 * att);
factor = stim.mat2Pa / sqrt(2);
%---------

resp_freq = zeros(stim.nLevels, numel(freq)); 
bline_freq = zeros(stim.nLevels, numel(freq)); 
for k = 1:stim.nLevels
    fprintf(1, 'Analyzing level # %d / %d ...\n', k, stim.nLevels);
    temp = reshape(squeeze(stim.resp(k, :, 2:end, 1:endsamps)),...
        (stim.nreps-1)*stim.Averages, endsamps);
    
    
    %----------------
    % Additional computation to get the Conductances
    energy = sum(temp.^2, 2);
    good = energy < (median(energy) + 0.5*mad(energy, 1));
    avg = mean(temp(good, :), 1) * -1 * factor;
    EarRespH =  rfft(avg').*exp(1j*2*pi*fcalib*delay) ./ Vo;
    Zec_mem = ldimp(Zs, Ps, EarRespH);
    Gmem(k, :) = interp1(fcalib, real(1./Zec_mem), freq);
    %----------------

    tempf = pmtm(temp', 4, freq, stim.Fs)';
    resp_freq(k, :) = median(tempf, 1);
    
    blevs = k; % Which levels to use as baseline (consider 1:k)
    temp2 = squeeze(stim.resp(blevs, :, 1, 1:endsamps));
    
    if(numel(blevs) > 1)
        temp2 = reshape(temp2, size(temp2, 2)*numel(blevs), endsamps);
    end
    
    %----------------
    % Additional computation to get the Conductances
    energy = sum(temp2.^2, 2);
    good = energy < (median(energy) + 2*mad(energy, 1));
    avg = mean(temp2(good, :), 1) * -1 * factor;
    EarRespH =  rfft(avg').*exp(1j*2*pi*fcalib*delay) ./ Vo;
    Zec_bline = ldimp(Zs, Ps, EarRespH);
    Gbline(k, :) = interp1(fcalib, real(1./Zec_bline), freq);
    %---------------
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


% Adjust for small calibration errors that sometimes lead to negative (but
% close to 0 errors to avoid numerical problems)
if any(Gbline(:) <= 0) || any(Gmem(:) <= 0)
    const = min([min(min(Gbline)), min(min(Gmem))]) - 1e-5;
    Gmem = Gmem - const;
    Gbline = Gbline -  const;
end


% Adjust measured ear-canal dB change with conductance change
MEMp = pow2db(resp_freq ./ bline_freq);
MEMg = softclip(pow2db(Gmem ./ Gbline));
MEM = MEMp + MEMg;

nfilt = 65;
for k = 1:stim.nLevels
    MEMs(k, :) = sgolayfilt(MEM(k, :), 2, nfilt);
    MEMps(k, :) = sgolayfilt(MEMp(k, :), 2, nfilt);
    MEMgs(k, :) = sgolayfilt(MEMg(k, :), 2, nfilt);
end

res.MEMs = MEMs; 
res.MEM = MEM; 
end 