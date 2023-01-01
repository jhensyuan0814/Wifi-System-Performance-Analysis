function IN = calculate_notransmitting_AP_IN(tAP, ntAP, shift, side)
    %AP_IN(i): AP(i)'s IN (inteference and noise) at the MS
    numberAP = length(ntAP);
    IN = zeros(1, numberAP);
    AP_collection = zeros(9, length(tAP));
    AP_collection(1, :) = tAP;
    for i = 2: 9
        AP_collection(i, :) = tAP + shift(i - 1);
    end
    AP_collection = reshape(AP_collection, 1, 9 * length(tAP));
    for i = 1: numberAP
        %wrap-around style to calculate SINR
        a = abs(real(AP_collection - ntAP(i))) < side / 2;
        b = abs(imag(AP_collection - ntAP(i))) < side / 2;
        now = AP_collection(a & b);
        g = sum(calculate_received_power(now, 1)); %power from AP
        noise = 1.38 * 10 ^ (-23) * 300 * 10 * 10 ^ 6;  %Noise = k * Temperature * Bandwidth
        IN(i) = noise + g;  %in Watt, not in dB
    end
end