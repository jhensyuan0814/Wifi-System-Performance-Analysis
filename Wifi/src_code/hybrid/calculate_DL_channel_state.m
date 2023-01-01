function [capacity, SINR] = calculate_DL_channel_state(MS, AP, APofMS, shift, side, MS_frequency, AP_frequency, bandwidth, frequency_reuse_factor, NumFDMA, downlink_buffer)
    %find out the AP every MS belongs to (APlist) and the Shannon capacity for each MS
    capacity = zeros(1, length(MS));
    SINR = zeros(1, length(MS));
    %for i = 1: frequency_reuse_factor
     %   for j = 1: NumFDMA
            %consider MS with the same frequency and is currently
            %transmitting
      %      now = (MS_frequency == (i - 1) * NumFDMA + j) & (downlink_buffer > 0);
       %     [capacity(now), SINR(now)] = calculate_DL_SINR(MS(now), AP(AP_frequency == i), AP(APofMS(now)), shift, side, bandwidth, frequency_reuse_factor, NumFDMA);
        %end
    %end
    for i = 1: length(AP)
        now = (APofMS == i);
        samef = sum(AP_frequency == AP_frequency(i));
        [capacity(now), SINR(now)] = calculate_DL_SINR(MS(now), AP(i), AP(APofMS(now)), shift, side, bandwidth, frequency_reuse_factor, NumFDMA);
        capacity(now) = capacity(now) / sum(now) / samef;
    end
end