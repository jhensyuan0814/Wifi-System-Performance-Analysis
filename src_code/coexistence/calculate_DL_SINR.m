function [capacity,SINR] = calculate_DL_SINR(MS,AP_collection,APNum,bandwidth,side, AP_frequency, MS_frequency,frequency_reuse_factor, NumFDMA, APofMS, active_APs, active_freq)
    %find out the AP every MS belongs to (APlist) and the Shannon capacity for each MS
    MSperAP = zeros(1,APNum);%counter recording the current number of MS under that AP
    forbidlist = [];
    capacity = [];% save the shannon capacity for each MS
    uplink_capacity = [];
    SINR = [];
    for i = 1: length(MS)
        %wrap-around style to calculate SINR
        lowerbound = MS(i) - side*0.5*(1+1j);
        upperbound = MS(i) + side*0.5*(1+1j);
        lx = real(lowerbound); ly = imag(lowerbound);
        ux = real(upperbound); uy = imag(upperbound);
        t_BS = []; t_idx = [];
        for j = 1:length(AP_collection)
            tx = real(AP_collection(j)); ty = imag(AP_collection(j));
            if(tx <= ux && tx >= lx && ty >= ly && ty<= uy)
                t_BS = [t_BS AP_collection(j)];
                if(mod(j,APNum) == 0)
                    t_idx = [t_idx APNum];
                else
                    t_idx = [t_idx mod(j,APNum)];
                end
            end
        end
        t_BS_calculate = [];
        for j = 1:length(active_APs)
            if active_freq(j) == MS_frequency(i)
                if sum(t_BS_calculate == t_BS(t_idx==active_APs(j)))==0
                    t_BS_calculate = [t_BS_calculate t_BS(t_idx==active_APs(j))];
                end
            end 
        end
        % only calculate sinr with the same FRF
        t_SINR = calculate_SINR( MS(i), t_BS_calculate, 0);
        if length(t_SINR)==0
            sinr = 0;
        else
            sinr = max(t_SINR);
        end
        SINR = [SINR sinr];
        capacity = [capacity bandwidth/(frequency_reuse_factor*NumFDMA)*log2(1 + sinr)];
    end
     SINR = pow2db(SINR);
end