function [new_MS_frequency, APlist, MSofAP,APofMS2] = find_AP(MS,AP_collection,APNum,side,MSlimit, AP_frequency,  NumFDMA)
    %find out the AP every MS belongs to (APlist) and the Shannon capacity for each MS
    MSperAP = zeros(1,APNum);%counter recording the current number of MS under that AP
    APofMS2 = zeros(1,APNum);
    forbidlist = [];
    APlist = []; % save the AP for each MS
    MSofAP = zeros(APNum*2,MSlimit);
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
                    t_idx = [t_idx 16];
                else
                    t_idx = [t_idx mod(j,APNum)];
                end
            end
        end
        t_SINR = calculate_SINR( MS(i), t_BS,0);
        %forbidlist is to limit the number of MS under each AP
        for m = 1:length(forbidlist)
            [~,idx2] = find(t_idx == forbidlist(m));
            t_SINR(idx2) = -1;
        end
        [~, selectedAP] = max(t_SINR);
        APlist = [APlist t_idx(selectedAP)];
        MSperAP(t_idx(selectedAP)) = MSperAP(t_idx(selectedAP))+1;
        [~,forbidlist] = find(MSperAP>=MSlimit);
    end
     new_MS_frequency = zeros(1, length(MS));
     % Round Robin assign frequency in every MS under one AP
     prev_allocation = zeros(1,length(APlist));
     for i = 1: length(MS)
         new_MS_frequency(i) = (AP_frequency(APlist(i)) - 1) * NumFDMA + (1+prev_allocation(APlist(i)));
         idx = 2*APlist(i)-1+prev_allocation(APlist(i));
         MSofAP(idx,find(MSofAP(idx,:)==0, 1, 'first')) = i;
         APofMS2(i) = idx;
         prev_allocation(APlist(i)) = mod(prev_allocation(APlist(i)) + 1, NumFDMA);
     end
end