function IN = calculate_MS_IN(MS, shift, side)
    %MS_IN(i): MS(i)'s IN (inteference and noise) at the MS
    numberMS = length(MS);
    IN = zeros(1, numberMS);
    MS_collection = zeros(9, numberMS);
    MS_collection(1, :) = MS;
    for i = 2: 9
        MS_collection(i, :) = MS + shift(i - 1);
    end
    MS_collection = reshape(MS_collection, 1, 9 * numberMS);
    for i = 1: numberMS
        %wrap-around style to calculate SINR
        a = abs(real(MS_collection - MS(i))) < side / 2;
        b = abs(imag(MS_collection - MS(i))) < side / 2;
        now = MS_collection(a & b);
        now(i) = []; %delete the current MS from now
        g = sum(calculate_received_power(now, 0)); %power from MS
        noise = 1.38 * 10 ^ (-23) * 300 * 10 * 10 ^ 6;  %Noise = k * Temperature * Bandwidth
        IN(i) = noise + g;  %in Watt, not in dB
    end
end