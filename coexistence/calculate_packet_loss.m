function [new_packet_loss, MS_pkt_loss] = calculate_packet_loss(currentMS, remaining_time, uplink_state)
    %update remaining time and calculate packet loss when timeout
    %in reality, this should not happen
    number = length(uplink_state);
    MS_pkt_loss = zeros(1,number);
    new_packet_loss = 0;
    for i = 1: number
        %transmitting, but timeout => packet loss
        if (currentMS(i) ~= 0) && (remaining_time(i) == 0) && (uplink_state(i) ~= 0)
            new_packet_loss = new_packet_loss + 1;
            MS_pkt_loss(i) = MS_pkt_loss(i) + 1;
            uplink_state(i) = 0;
        end
    end
end