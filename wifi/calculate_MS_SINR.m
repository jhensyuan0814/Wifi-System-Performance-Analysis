function SINR = calculate_MS_SINR(MS, AP, APofMS,shift,side)  %all inputs are positions
    %AP_SINR(i): MS(i)'s SINR at the correponding AP
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
    %display(APofMS);
    %display(AP);
    %display(ID(0));
    for i = 1: numberMS
        %wrap-around style to calculate SINR
        a = abs(real(AP_collection - MS(i))) < side / 2;
        b = abs(imag(AP_collection - MS(i))) < side / 2;
        t_AP = AP_collection(a & b);
        t_sinr = calculate_received_power(MS(i)- t_AP, 0);
        tot_uplink = sum(t_sinr,'all');
        SINR(i) = t_sinr(ID(i))/(tot_uplink-t_sinr(ID(i))+noise);
    end
end