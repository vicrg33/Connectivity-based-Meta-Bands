function AEC = calculo_AEC(signal1)
% Computes the Amplitude Envelope Correlation (AEC)

% Input: M/EEG signal (nSamples*nChannel)
% Output: AEC values

%% Creamos el vector de salida con las correlaciones
signal2 = signal1;
AEC=NaN(1,size(signal1,2));

% Estimate the envelope of each channel
hilb1 = hilbert(signal1(:,:));
hilb2 = hilbert(signal2(:,:));
envelope1 = squeeze(abs(hilb1));
envelope2 = squeeze(abs(hilb2));
env1=log(envelope1.^2);
env2=log(envelope2.^2);

% Calculate the correlation between each pair of channels
AEC = corr(env1,env2);
end
