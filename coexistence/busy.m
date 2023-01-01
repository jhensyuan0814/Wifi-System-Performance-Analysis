function [new_buffer, new_packet_loss, new_uplink_state, new_MS_collision, remaining_data, new_currentMS, MS_pkt_loss] = busy(currentMS, buffer, capacity, AP_SINR, sinr_thres,  MS_collision, payload_size, uplink_state, remaining_data, remaining_time)
    %AP is transmitting data
    %need to update packet loss and buffer of MS
    %payload_size: the size of a packet
    %uplink_state == 0: not transmitting any packet now
    %uplink_state == 1: collision have not occurred
    %uplink_state == 2: collision occurs, meaning that the trasmitted
    %packets should be calculated as packet loss
    %remaining_data: remaining data for the current packet
    number = length(buffer);  %number of MSs served by the AP
    new_uplink_state = uplink_state;
    new_currentMS = currentMS;
    new_packet_loss = 0; 
    remaining_capacity = capacity;  %initialize capacity
    new_buffer = buffer;
    new_MS_collision = MS_collision;
    MS_pkt_loss = zeros(1,number);
    %remaining_time is not enough
    for i = 1: number
        %unit of buffer_size: data number, not packet number
        %ready to transmit new packet (new_current_loss(i) == 2)
        if new_uplink_state(i) == 0 && remaining_time(i) * remaining_capacity(i) > payload_size && buffer(i) > 0
            new_uplink_state(i) = 1;
            new_buffer(i) = new_buffer(i) - 1;
            remaining_data(i) = payload_size;
        end
        %not enough time to transmit another packet or there are not
        %packet => stop transmitting
        if new_uplink_state(i) == 0
            new_currentMS(i) = 0;
        end
        if new_uplink_state(i) > 0
            %AP: collision occurs
            if AP_SINR(i) < sinr_thres
                new_uplink_state(i) = 2;
                new_MS_collision(i) = 1;
            end
            %MS thinks that collision occurs
            %if MS_IN(i) > Ith && new_MS_collision(i) == 0
            %   new_MS_collision(i) = 1;
            %end
                
            remaining_data(i) = max(0, remaining_data(i) - floor(remaining_capacity(i)));
            %finish transmitting current data
            if remaining_data(i) == 0
                if (new_uplink_state(i) == 2)  %packet loss
                    new_packet_loss = new_packet_loss + 1;
                    MS_pkt_loss(i) = 1;
                    
                end
                new_uplink_state(i) = 0;
            end
        end
    end

    
end