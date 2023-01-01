function [a, b, c, d, e, f] = downlink_update_timer(MS_collision, windowN, backoff_timer, DIFS_timer, resultMS, backoffreset, downlink_buffer, remaining_time,downlinkI, APofMS2,Ith, DIFSth, initial, max_time, CWminN, CWmaxN, slottime)
%this function updates the contention window size(windowN), the backoff timers, and which APs are transmitting(currentAP)
     %input interference instead of SINR
   if (initial==1) %initialization
      windowN = windowN * CWminN;
      backoff_timer = randi([0,2^CWminN-1],1,numel(downlinkI))*slottime;
      DIFS_timer = ones(1,numel(downlinkI))*DIFSth;
      backoffreset = zeros(1,numel(downlinkI));
      resultMS = zeros(1,numel(downlinkI));
   else
       for i = 1:length(downlinkI) %downlinkI can be switched to downlink interference
            if(resultMS(i)==1) %AP is transmitting
                remaining_time(i) =  remaining_time(i) -1;
                if(remaining_time(i)==0)
                    resultMS(i) = 0;
                end
            else
               is_collision =  sum(MS_collision(APofMS2==i));
               has_pkt = sum(downlink_buffer(APofMS2==i));
               if(DIFS_timer(i)>0 && has_pkt && downlinkI(i)<=Ith) %still DIFS+has pkt to send+channel idle
                    if(backoffreset(i)==1)
                       if(is_collision>=1) %last round:collision
                            if(windowN(i) ~= CWmaxN)
                                 windowN(i) = windowN(i) + 1;
                            end
                            backoff_timer(i) = randi([0,2^ windowN(i)-1])*slottime;
                            MS_collision(APofMS2==i)=0;
                            %fprintf("collision %d \n",i);
                       else
                            windowN(i) = CWminN; 
                            backoff_timer(i) = randi([0,2^CWminN-1])*slottime;
                            MS_collision(APofMS2==i)=0;
                            %fprintf("success %d \n",i);
                       end
                       backoffreset(i) = 0;
                       remaining_time(i)=0;
                    end
                   DIFS_timer(i) = DIFS_timer(i) -1;
               elseif(DIFS_timer(i)>0 && has_pkt && downlinkI(i)>Ith)%still DIFS+has pkt to send+channel busy
                   DIFS_timer(i) = DIFSth;
               elseif(DIFS_timer(i)==0)%contention
                   if(backoff_timer(i) == 0)
                           resultMS(i) = 1;
                           remaining_time(i) = max_time;
                           DIFS_timer(i) = DIFSth;
                           backoffreset(i) = 1;
                   else 
                        if(downlinkI(i) <= Ith) %channel is idle
                           backoff_timer(i) = backoff_timer(i) -1;
                        else
                            DIFS_timer(i) = DIFSth;
                        end
                       %channel is busy->pause the backoff_timer
                   end 
               end
           end
       end
   end
   a = windowN;
   b = backoff_timer;
   c = DIFS_timer;
   d = resultMS;
   e = backoffreset;
   f = remaining_time;
end