function IN = calculate_AP_IN(AP, shift, side)
    %MS_IN(i): MS(i)'s IN (inteference and noise) at the MS
    numberMS = length(AP);
    IN = zeros(1, numberMS);
    AP_collection = zeros(9, numberMS);
    AP_collection(1, :) = AP;
    for i = 2: 9
        AP_collection(i, :) = AP + shift(i - 1);
    end
    AP_collection = reshape(AP_collection, 1, 9 * numberMS);
    for i = 1: numberMS
        %wrap-around style to calculate SINR
        a = abs(real(AP_collection - AP(i))) < side / 2;
        b = abs(imag(AP_collection - AP(i))) < side / 2;
        now = AP_collection(a & b);
        now(i) = []; %delete the current AP from now
        g = sum(calculate_received_power(now, 1)); %power from AP
        noise = 1.38 * 10 ^ (-23) * 300 * 10 * 10 ^ 6;  %Noise = k * Temperature * Bandwidth
        IN(i) = noise + g;  %in Watt, not in dB
    end
end