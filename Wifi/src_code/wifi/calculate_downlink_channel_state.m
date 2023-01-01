function [new_MS_SINR, new_AP_IN] = calculate_downlink_channel_state(MS, AP, APofMS, frequency_reuse_factor, NumFDMA,currentAP,shift,side,APtotfreq,MSofAP)
    %calculate MS_SINR and AP_IN
    %APtotfreq:16*2 map, AP to freq(1~6)
    %MSofAP:32*MSlimit, map channel to MS
    %currentAP:32*1
    %MS_SINR:160*1
    %AP_IN:32*1
    %calculate MS_SINR
    new_MS_SINR = zeros(1, length(MS));
    currentAP = reshape(currentAP,NumFDMA,length(AP))';
    transmitting = (currentAP == 1); %16*2
    for i = 1: frequency_reuse_factor
        for j = 1: NumFDMA
            %consider MS with the same frequency and is currently
            %transmitting
            now = (APtotfreq == (i - 1) * NumFDMA + j) & transmitting; %16*2
            now2 = any(now, 2)'; %1*16
            %now2 = now(:,1)|now(:,2);
            %now2 = now2'; %1*16
            now = reshape(now',1,[]); %1*32
            for k = 1:length(now)
                if(now(k))
                    idx = MSofAP(k,:);%ms under ap(k)
                    %display(idx(idx~=0));
                    new_MS_SINR(idx(idx~=0)) = calculate_MS_SINR(MS(idx(idx~=0)), AP(now2), AP(APofMS(idx(idx~=0))),shift,side);%mistake????
                end
            end
        end
    end
    %calculate AP_IN
    new_AP_IN = zeros(1,length(AP)*2);%1*32
    for i = 1: frequency_reuse_factor
        for j = 1: NumFDMA
            %consider MS with the same frequency and is currently
            %transmitting
            now = (APtotfreq == (i - 1) * NumFDMA + j) & transmitting;
            noww = now(:,1)|now(:,2);
            noww = noww';
            now = reshape(now',1,[]); %1*32

            new_AP_IN(now) = calculate_AP_IN(AP(noww), shift, side); %transmitting AP 
            now2 = (APtotfreq == (i - 1) * NumFDMA + j) & (~transmitting);
            noww2 = now2(:,1)|now2(:,2);
            noww2 = noww2';
            now2 = reshape(now2',1,[]); %1*32
            new_AP_IN(now2) = calculate_notransmitting_AP_IN(AP(noww), AP(noww2), shift, side); %not transmitting AP
        end
    end
    currentAP = reshape(currentAP',1,[]);
   % if(length(find(currentAP))~=0)
      %  display(currentAP);
      %  display( new_MS_SINR);
      %  pause(100);
    %end
end