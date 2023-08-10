function [rel_bp,n_chan]=calc_band_power(pxx,freq_bands)

fi= 1:0.5:((size(pxx,2)+1)/2);

% check pxx and freq dimensions)
if size(pxx,2) ~= length(fi)
    error('Number of columns in pxx must match the length of freq')
end

% number of channels
n_chan = size(pxx,1);

% number of frequency bands
n_bands = size(freq_bands,1);

% initialise array for storing band power
bp = zeros(n_chan,n_bands);

% compute band power from PSD
for b=1:n_bands
    bp(:,b) = bandpower(pxx',fi,freq_bands(b,:),'psd');
end

% log transform
bp = log10(1+bp);

% relative band power
rel_bp = bp./repmat(sum(bp,2),1,n_bands);

end