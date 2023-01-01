function [new_buffer, new_packet_loss] = idle(buffer, buffer_size, payload_size,packet,remaining_packet)
    %put coming packets into buffer
    %need to update packet loss and buffer of MS
    %payload_size: the size of a packet
    %remaining_data: remaining data for the current packet
    new_packet_loss = 0; 
    new_buffer = buffer;
    number = length(buffer);  %number of MSs served by the AP
    %remaining_time is not enough
    for i = 1: number
        %put new packets into buffer
        %unit of buffer_size: data number, not packet number
        if packet(i) > 0  %buffer_size is for each MS
            previous = packet(i);
            remaining_buffer_size = buffer_size - new_buffer(i) * payload_size - remaining_packet(i);
            packet(i) = max(0, packet(i) - fix(remaining_buffer_size/ payload_size));
            new_buffer(i) = new_buffer(i) + (previous - packet(i));
        end
        if packet(i) > 0
            new_packet_loss = new_packet_loss + packet(i);
        end
    end
end