% Note: "temporal_activation_sequence" is "frequency_activation_sequence"
% but with the original notation by Nuñez et al. 2021 (Neuroimage).

%% Attraction strength

for ii = 1:size(temporal_activation_sequence,1)
    if remove_50_hz
        for jj = 1:freq_bins
            if ~ismember(jj,freqs_artifact)
                AS_curve(jj,ii) = instantaneous_correlation_tensor(temporal_activation_sequence(ii,jj),jj,ii);
            else % The frequencies in the removed (due to the artifact) band, 
                % will be set to the average, but they have no importance (AS do not exist in these frequencies)
                AS_curve(jj,ii) = NaN;
            end
        end
    else
        for jj = 1:freq_bins
            AS_curve(jj,ii) = instantaneous_correlation_tensor(temporal_activation_sequence(ii,jj),jj,ii);            
        end
    end
end
%AS_curve = (nanmean(AS_curve,2) - nanmean(nanmean(AS_curve,2))) / nanstd(nanmean(AS_curve,2));    
AS_curve = nanmean(AS_curve,2);


%% Degree of dominance

total = sum(instantaneous_correlation_tensor(:));
for ii = 1:size(temporal_activation_sequence,1)
    if remove_50_hz
        for jj = 1:freq_bins
            if ~ismember(jj,freqs_artifact)
                max_corr = instantaneous_correlation_tensor(temporal_activation_sequence(ii,jj),jj,ii);
                min_corr = instantaneous_correlation_tensor(:,jj,ii);
                min_corr(temporal_activation_sequence(ii,jj)) = [];
                min_corr = mean(min_corr);
                DoD_curve(jj,ii) = max_corr - min_corr;
            else % The frequencies in the removed (due to the artifact) band, 
                % will be set to the average, but they have no importance (AS do not exist in these frequencies)
                DoD_curve(jj,ii) = NaN;
            end
        end            
    else
        for jj = 1:freq_bins
            corrs = instantaneous_correlation_tensor(:,jj,ii);
            DoD_curve(jj,ii) = sum(corrs) / total;
        end
    end
end
DoD_curve = nanmean(DoD_curve,2);