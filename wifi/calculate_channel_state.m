
function [new_AP_SINR, new_MS_IN] = calculate_channel_state(MS, AP, APofMS, MS_frequency, AP_frequency, frequency_reuse_factor, NumFDMA,currentMS,shift,side)
    %calculate AP_SINR and MS_IN
    %calculate AP_SINR
    new_AP_SINR = zeros(1, length(MS));
    transmitting = (currentMS == 1);
    for i = 1: frequency_reuse_factor
        for j = 1: NumFDMA
            %consider MS with the same frequency and is currently
            %transmitting
            now = (MS_frequency == (i - 1) * NumFDMA + j) & transmitting;
            new_AP_SINR(now) = calculate_AP_SINR(MS(now), AP(AP_frequency == i), AP(APofMS(now)),shift,side);
        end
    end
    %calculate MS_IN
    new_MS_IN = zeros(1, length(MS));
    for i = 1: frequency_reuse_factor
        for j = 1: NumFDMA
            %consider MS with the same frequency and is currently
            %transmitting
            now = (MS_frequency == (i - 1) * NumFDMA + j) & transmitting;
            new_MS_IN(now) = calculate_MS_IN(MS(now), shift, side); %transmitting MS 
            now2 = (MS_frequency == (i - 1) * NumFDMA + j) & (~transmitting);
            new_MS_IN(now2) = calculate_notransmitting_MS_IN(MS(now), MS(now2), shift, side); %not transmitting MS
        end
    end
end