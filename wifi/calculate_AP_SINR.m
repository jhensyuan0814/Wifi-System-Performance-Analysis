function SINR = calculate_AP_SINR(MS, AP, APofMS,shift,side)  %all inputs are positions
    %AP_SINR(i): MS(i)'s SINR at the correponding AP
    numberAP = length(AP);
    numberMS = length(MS);
    SINR = zeros(1, numberMS);
    noise = 1.38 * 10 ^ (-23) * 300 * 10 * 10 ^ 6; 
    MS_collection = zeros(9, numberMS);
    MS_collection(1, :) = MS;
    for i = 2: 9
        MS_collection(i, :) = MS + shift(i - 1);
    end
    MS_collection = reshape(MS_collection, 1, 9 * numberMS);
    [~, ID] = ismember(APofMS, AP); 
    for i = 1: numberAP
        %wrap-around style to calculate SINR
        a = abs(real(MS_collection - AP(i))) < side / 2;
        b = abs(imag(MS_collection - AP(i))) < side / 2;
        t_MS = MS_collection(a & b);
        t_sinr = calculate_received_power(AP(i)- t_MS,1);
        tot_uplink = sum(t_sinr,'all');
        for k = 1:length(MS)
            if(ID(k)==i)
                SINR(k) = t_sinr(k)/(tot_uplink-t_sinr(k)+noise);
            end
        end
    end
end