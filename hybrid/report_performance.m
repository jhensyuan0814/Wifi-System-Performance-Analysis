function report_performance(total_packet, buffer, buffer_size, payload_size, packet_loss, avg_curMS, mode)
    number = length(buffer);
    disp(['total packet: ', num2str(total_packet)]);
    disp(['total transmitted packet: ', num2str(total_packet - sum(buffer))]);
    disp(['total collision packet: ', num2str(packet_loss)]);
    disp(['packet loss rate: ', num2str(packet_loss / (total_packet - sum(buffer)) * 100), '%']);
    if mode == 1
        disp(['total buffer packet: ', num2str(sum(buffer))]);
        disp(['buffer utilization rate: ', num2str(sum(buffer) * payload_size/ buffer_size / number * 100), '%']);
        disp(['number of transmitting MS: ', num2str(avg_curMS)]);
    end
end