function [capacity,SINR] = calculate_LTE_DL_SINR(MS, AP, APofMS, shift, side, bandwidth, frequency_reuse_factor, NumFDMA)
    %find out the AP every MS belongs to (APlist) and the Shannon capacity for each MS
    numberAP = length(AP);
    numberMS = length(MS);
    SINR = zeros(1, numberMS);
    noise = 1.38 * 10 ^ (-23) * 300 * 10 * 10 ^ 6; 
    AP_collection = zeros(9, numberAP);
    AP_collection(1, :) = AP;
    for i = 2: 9
        AP_collection(i, :) = AP + shift(i - 1);
    end
    AP_collection = reshape(AP_collection, 1, 9 * numberAP);
    [~, ID] = ismember(APofMS, AP);
    tot_downlink = zeros(1, numberMS);
    for i = 1: numberMS
        %wrap-around style to calculate SINR
        a = abs(real(AP_collection - MS(i))) < side / 2;
        b = abs(imag(AP_collection - MS(i))) < side / 2;
        t_AP = AP_collection(a & b);
        t_sinr = calculate_received_power(MS(i)- t_AP, 0);
        tot_downlink(i) = sum(t_sinr,'all');
    end
    for k = 1: numberMS
        %wrap-around style to calculate SINR
        a = abs(real(AP_collection - MS(i))) < side / 2;
        b = abs(imag(AP_collection - MS(i))) < side / 2;
        t_AP = AP_collection(a & b);
        a = calculate_received_power(MS(k) - t_AP(ID(k)), 0);
        SINR(k) = a /(tot_downlink(k)- a + noise);
    end
    capacity = (bandwidth / frequency_reuse_factor / NumFDMA) * log2 (1 + SINR);
end