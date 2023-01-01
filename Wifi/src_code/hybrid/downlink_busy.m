function [new_buffer, new_packet_loss, new_downlink_collision, new_remaining_data] = downlink_busy(buffer, capacity, downlink_SINR, sinr_thres, downlink_collision, payload_size, remaining_data)
    %AP is transmitting data
    %need to update packet loss and buffer of MS
    %payload_size: the size of a packet
    %uplink_state == 0: not transmitting any packet now
    %uplink_state == 1: collision have not occurred
    %uplink_state == 2: collision occurs, meaning that the trasmitted
    %packets should be calculated as packet loss
    %remaining_data: remaining data for the current packet
    number = length(buffer);  %number of MSs served by the AP
    new_packet_loss = 0; 
    remaining_capacity = capacity;  %initialize capacity
    new_buffer = buffer;
    new_downlink_collision = downlink_collision;
    new_remaining_data = remaining_data;
    %remaining_time is not enough
    for i = 1: number
        %unit of buffer_size: data number, not packet number
        %ready to transmit new packet (new_current_loss(i) == 2)
        if new_remaining_data(i) == 0 && buffer(i) > 0
            new_buffer(i) = new_buffer(i) - 1;
            new_remaining_data(i) = payload_size;
            new_downlink_collision(i) = 0;
        end
        
        if new_remaining_data(i) > 0
            %MS: collision occurs
            if downlink_SINR(i) < sinr_thres
                new_downlink_collision(i) = 1;
            end
           
            new_remaining_data(i) = max(0, new_remaining_data(i) - floor(remaining_capacity(i)));
            %finish transmitting current data
            if new_remaining_data(i) == 0
                if (new_downlink_collision(i) == 1)  %packet loss
                    new_packet_loss = new_packet_loss + 1;
                end
                new_downlink_collision(i) = 0;
            end
        end
    end
end