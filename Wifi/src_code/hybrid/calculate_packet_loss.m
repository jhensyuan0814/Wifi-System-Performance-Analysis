function new_packet_loss = calculate_packet_loss(currentMS, remaining_time, uplink_state)
    %update remaining time and calculate packet loss when timeout
    %in reality, this should not happen
    number = length(uplink_state);
    new_packet_loss = 0;
    for i = 1: number
        %transmitting, but timeout => packet loss
        if (currentMS(i) ~= 0) && (remaining_time(i) == 0) && (uplink_state(i) ~= 0)
            new_packet_loss = new_packet_loss + 1;
            uplink_state(i) = 0;
        end
    end
end