function [signal_ort] = ortogonalization_optimizada(signal1)
% Orthogonalize each channel regarding the other channels. Based on O'Neill
% et al. (2015)

signal2 = signal1;
signal_ort=NaN(size(signal1,1),size(signal1,2),size(signal1,2));

for channel1=1:size(signal1,2)
    signal_ort(:,channel1,channel1)=signal1(:,channel1);
    for channel2=1:size(signal1,2)
        if channel1~=channel2
            signal1_canal=signal1(:,channel1);
            signal2_canal=signal2(:,channel2);
            beta=signal2_canal.'*pinv(signal1_canal.');
            signal_ort(:,channel2,channel1)=signal2_canal-beta*signal1_canal;
        end
    end
end
end