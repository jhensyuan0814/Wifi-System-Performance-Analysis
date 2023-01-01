function g = calculate_received_power(d, state)
    %calculate received power when the distance is d (an array)
    ht = 51.5;  %the height of the base station
    hi = 1.5;  %the height of the mobile station
    if(state == 0)
        Ptx = 3;  %the power of the AP, in dB
    elseif(state ==1)
        Ptx = -7; %the power of the MS, in dB
    end
    Gtx = 14;  %the power of the transimitter gain, in dB
    Grx = 14;  %the power of the receiver gain, in dB
    numerator = (hi * ht) ^ 2;
    g = numerator * abs(d) .^ (-4);  %two-ray ground model
    g = g * db2pow( Ptx + Gtx + Grx);  %in Watt
end