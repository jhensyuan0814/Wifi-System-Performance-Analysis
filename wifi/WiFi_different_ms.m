clear all;
clear plt;

%set some parameters


mscoll = [6:2:16];
packet_loss = zeros(4, length(mscoll));
throughput = zeros(4, length(mscoll));
buffutil = zeros(4, length(mscoll));
for lm = 1: length(mscoll)
    %parameters for busy and idle
    bandwidth = 10 ^ 7;
    side = 1000; %length of square boundary 
    APNum = 16; %total number of APs
    MSperCell = mscoll(lm);  %average number of MSs per AP
    MSNum = MSperCell * APNum;  %total number of MSs
    downlink_buffer_size = 10 ^ 5;  %
    uplink_buffer_size = 10^5;%292*10 ^ 5;  %
    boxNum = 9; %9 box wrap around
    MSlimit = 10; %the maximum number of MSs under one AP
    sinr_thres = 10;%0.01;%-10000; %the SINR threshold for collision
    slot_time = 2;
    DIFSth =  5; %128 / 28
    PIFSth = 3;
    SIFSth = 1; %28
    Ith = 0.01;
    payload_size = 500;%292; %8184/28, the size of a packet, may convert to uniform distribution
    max_time = 10; %the maximum time an AP can transmit
    CWminN = 4; %contention window parameter- 2^CWminN
    CWmaxN = 7; 
    %(1,5):packet loss = 4.5%  transmitting = 7
    %(4,7):packet loss = 8.4%  transmitting = 10
    %initial construction for our scenario
    %AP = 2*side*(complex(rand(1, APNum), rand(1,APNum))-0.5); %random version
    rowNum = sqrt(APNum);
    gap = side/(2 * rowNum);
    [X, Y] = meshgrid([-rowNum + 1: 2: rowNum - 1] * gap, [-rowNum + 1: 2: rowNum - 1] * gap);
    AP = X + 1j * Y;
    AP = reshape(AP, [1, APNum]); 
    noiseamp = 0.1;
    AP = AP + noiseamp * side * (complex(rand(1, APNum) - 0.5, rand(1, APNum) - 0.5));  %uniform version with skews
    MS = side * (complex(rand(1, MSNum)-0.5, rand(1, MSNum) - 0.5));
    labels = reshape(split(num2str(1: 1: APNum)), 1, APNum);
    %figure(1); plot(MS,'x'); title('Overview'); xlabel('x-axis(m)'); ylabel('y-axis(m)'); hold on;
    %plot(AP,'ro'); hold on;
    %text(real(AP), imag(AP), labels, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left'); hold off;

    %replicate the information to 9 boxes
    [X2, Y2] = meshgrid([-side 0 side], [-side 0 side]);
    shift =  reshape(X2+1j*Y2,[1,9]);
    shift = shift(shift ~= (0+0j));
    AP_collection = AP; MS_collection = MS;
    for i = 1: boxNum - 1
        AP_collection = [AP_collection AP + shift(i)];
        MS_collection = [MS_collection MS + shift(i)];
    end
    clearvars X Y X2 Y2 MS_collection;
    downlink_lambda = 0.005;
    uplink_lambda = 0.005;%1000;%1000
    uplink_windowN = ones(1,MSNum);
    uplink_backoff_timer = zeros(1,MSNum);
    downlink_packet_loss = 0;
    uplink_packet_loss = 0;
    poll_UL_packet_loss = 0;
    poll_DL_packet_loss = 0;
    poll_total_UL_packet = 0;
    poll_total_DL_packet = 0;
    downlink_total_packet = 0;
    uplink_total_packet = 0;
    uplink_curMS = 0;
    downlink_curAP = 0;
    downlink_buffer = zeros(1, MSNum);  %downlink_buffer(i): the buffer of MS(i) in downlink
    uplink_buffer = zeros(1, MSNum);  %downlink_buffer(i): the buffer of MS(i) in downlink
    downlink_packet = zeros(1, MSNum);  %packet(i): the number of packets of MS(i)
    uplink_packet = zeros(1, MSNum);  %uplink_packet(i): the number of packets of MS(i) in uplink
    currentAP = zeros(1, APNum);  %currentAP: if AP(i) is transmitting, currentAP(i) = 1
    currentMS = zeros(1, MSNum);  %currentMS: if MS(i) is transmitting, currentMS(i) = 1  

    %some attributes
    %uplink_state == 0: not transmitting any packet now
    %uplink_state == 1: collision have not occurred
    %uplink_state == 2: collision occurs, meaning that the trasmitted
    %packets should be calculated as packet loss
    uplink_state = zeros(1, MSNum);
    downlink_state= zeros(1, MSNum);
    %frequency
    frequency_reuse_factor = 3;
    NumFDMA = 2;
    %AP_frequency = randi([1, frequency_reuse_factor], 1, APNum); %1, 2, or 3
    AP_frequency = 1:APNum;
    AP_frequency = mod(AP_frequency,3)+1;
    APtotfreq =[AP_frequency'*2-1,AP_frequency'*2];
    MS_frequency = zeros(1, MSNum); %1, 2, ..., or 6
    %movable MS
    walk_t = zeros(1, MSNum); %walk_t(i): MS(i) move for how long
    velocity = zeros(1, MSNum); %velocity(i): MS(i)'s conjugate velocity
    scale_v = 0.001; %scale the velocity
    max_velocity = 6;
    max_walktime = 6;
    %SINR and collision
    AP_SINR = zeros(1, MSNum);  %each MS's SINR at the AP
    MS_IN = zeros(1, MSNum);  %each MS's IN at the MS
    AP_IN = zeros(1, APNum*NumFDMA);
    MS_collision = zeros(1, MSNum); %MS_collision(i): MS(i) think the packet has collision
    remaining_time = zeros(1, MSNum); %remaining_time(i): remaining time of a current transmitting MS(i)
    uplink_remaining_packet = zeros(1, MSNum); %remaining_packet(i): remaining packet size MS(i) is transmitting
    downlink_remaining_packet = zeros(1, MSNum); %remaining_packet(i): remaining # of packet MS(i) is transmitting [0.1)
    remaining_poll_UL_packet = zeros(1, MSNum); %remaining_poll_UL_packet(i): # of remaining poll packet  MS(i) is transmitting [0,1)
    remaining_poll_DL_packet = zeros(1, MSNum); %remaining_poll_UL_packet(i): # of remaining poll packet  MS(i) is transmitting [0,1)
    counted_dl_loss = zeros(1,MSNum);
    counted_UL_poll_loss = zeros(1,MSNum);
    counted_DL_poll_loss = zeros(1,MSNum);
    %update timer
    backoffreset = zeros(1,2*APNum);
    backoff_timer = zeros(1,2*APNum);
    DIFS_timer = zeros(1,2*APNum);
    downlink_DIFS_timer = zeros(1,2*APNum);
    downlink_backoffreset = zeros(1,2*APNum);
    downlink_backoff_timer = zeros(1,2*APNum);
    downlink_windowN = ones(1,2*APNum);
    downlink_remaining_time = zeros(1,2*APNum);
    remain_timeDL = zeros(1,MSNum);
    APofMS2 = [];
    AP_collision = zeros(1,MSNum);
    [uplink_windowN, backoff_timer,  DIFS_timer,  currentMS,  backoffreset] = update_timer(MS_collision, uplink_windowN, backoff_timer,  DIFS_timer,  currentMS,  backoffreset, downlink_buffer, remaining_time,MS_IN, Ith, DIFSth,1, max_time,CWminN, CWmaxN, slot_time);%initialize
    [downlink_windowN, downlink_backoff_timer,  downlink_DIFS_timer,  currentAP,  downlink_backoffreset, downlink_remaining_time] =  downlink_update_timer(AP_collision, downlink_windowN, downlink_backoff_timer,  downlink_DIFS_timer,  currentAP,  downlink_backoffreset, downlink_buffer, downlink_remaining_time,AP_IN, APofMS2,Ith, DIFSth,1, max_time,CWminN, CWmaxN, slot_time);

    simulation_time = 0;iteration_time =  25*10^-6;
    DL_time = 1; %downlink throughtput is greater than uplink
    UL_time = 1;
    PCF_UL_time = 1;
    PCF_DL_time = 1;
    dl_prob_thres = 1; % dl_prob_thres% probability to transmit DL PCF data
    ul_prob_thres = 1; % ul_prob_thres% probability to transmit UL PCF data



    while simulation_time < 5000
        simulation_time = simulation_time + 1;
        %move MS
        for i = 1: MSNum
            if (walk_t(i) == 0)
                walk_t(i) = randi([1, max_walktime]);
                velocity(i) = scale_v * randi([1, max_velocity])  * exp (1j*rand(1) * 2 * pi);
            end
        end
        walk_t = walk_t - 1;
        MS = MS + velocity;
        MS = in_boundary(MS, side); %deal with out of boundary

        %packet arrival
        downlink_packet = poissrnd(downlink_lambda, MSNum,1);  %packet arrival
        uplink_packet = poissrnd(uplink_lambda, MSNum,1);  %packet arrival
        for i = 1: MSNum %if buffer is full, forbid generating new packets
            if (uplink_packet(i) + uplink_buffer(i)) * payload_size + uplink_remaining_packet(i) > uplink_buffer_size  
                uplink_packet(i) = 0;
            end
        end
        downlink_total_packet = downlink_total_packet + sum(downlink_packet, 'all');  %calculate total number of packets
        uplink_total_packet = uplink_total_packet + sum(uplink_packet, 'all');  %calculate total number of packets
        [downlink_buffer, temp1] = idle(downlink_buffer, downlink_buffer_size, payload_size, downlink_packet, uplink_remaining_packet);
        downlink_packet_loss = downlink_packet_loss + temp1; 
        [uplink_buffer, temp1] = idle(uplink_buffer, uplink_buffer_size, payload_size,uplink_packet,uplink_remaining_packet);
        uplink_packet_loss = uplink_packet_loss + temp1;    
        %calculate SINR and deal with handoff
        if(mod(simulation_time,100)==1)
            [MS_frequency, APofMS,MSofAP,APofMS2] = find_AP(MS,AP_collection,APNum,side,MSlimit, AP_frequency,  NumFDMA);
        end

        %downlink
        if mod(simulation_time, 4) == 0
            %First, calculate channel state
            [MS_SINR,AP_IN]=calculate_downlink_channel_state(MS, AP, APofMS, frequency_reuse_factor, NumFDMA,currentAP,shift,side,APtotfreq,MSofAP);
            %Second, push all the new packets into buffer
            downlink_capacity = bandwidth/(NumFDMA*frequency_reuse_factor)*log2(1+MS_SINR);
            downlink_capacity = iteration_time*downlink_capacity*DL_time;
            %Third, transmit packets in the buffer
            [~,tidx] = find(currentAP == 1);        
            currentDL = ismember(APofMS2, tidx);
            for i = 1: length(downlink_remaining_time)
                [~,tidx2] = find(APofMS2 == i);
                remain_timeDL(tidx2) = downlink_remaining_time(i);
                downlink_capacity(tidx2) = downlink_capacity(tidx2) / length(tidx2);
            end
            [downlink_buffer(currentDL), temp2, downlink_state(currentDL), AP_collision(currentDL), downlink_remaining_packet(currentDL),currentDL(currentDL)] = busy(currentDL(currentDL), downlink_buffer(currentDL), downlink_capacity(currentDL), MS_SINR(currentDL), sinr_thres,  AP_collision(currentDL), payload_size, downlink_state(currentDL), downlink_remaining_packet(currentDL), remain_timeDL(currentDL));
            for i = 1: length(currentAP)
                currentAP(i) = sum(currentDL(APofMS2 == i)) > 0;
            end
            temp3 = calculate_packet_loss(currentDL, remain_timeDL, downlink_state);
            downlink_packet_loss = downlink_packet_loss + temp2 + temp3;
            %Fourth, update timer
            [downlink_windowN, downlink_backoff_timer,  downlink_DIFS_timer,  currentAP,  downlink_backoffreset, downlink_remaining_time] =  downlink_update_timer(AP_collision, downlink_windowN, downlink_backoff_timer,  downlink_DIFS_timer,  currentAP,  downlink_backoffreset, downlink_buffer, downlink_remaining_time,AP_IN, APofMS2,Ith, DIFSth,0, max_time,CWminN, CWmaxN, slot_time); 
        %uplink
        elseif mod(simulation_time, 4) == 1
            %First, calculate channel state
            [AP_SINR,MS_IN]=calculate_channel_state(MS, AP, APofMS, MS_frequency, AP_frequency, frequency_reuse_factor, NumFDMA,currentMS,shift,side);
            %Second, push all the new packets into buffer
            uplink_capacity = bandwidth/(NumFDMA*frequency_reuse_factor)*log2(1+AP_SINR);
            uplink_capacity = iteration_time*uplink_capacity*UL_time;
            %Third, transmit packets in the buffer
            [uplink_buffer(currentMS == 1), temp2, uplink_state(currentMS == 1),MS_collision(currentMS == 1), uplink_remaining_packet(currentMS == 1),currentMS(currentMS == 1)] = busy(currentMS(currentMS == 1), uplink_buffer(currentMS == 1), uplink_capacity(currentMS == 1), AP_SINR(currentMS == 1), sinr_thres,  MS_collision(currentMS == 1), payload_size, uplink_state(currentMS == 1), uplink_remaining_packet(currentMS == 1), remaining_time(currentMS == 1));
            temp3 = calculate_packet_loss(currentMS, remaining_time, uplink_state);
            uplink_packet_loss = uplink_packet_loss + temp2 + temp3;
            %Fourth, update timer
            [uplink_windowN, backoff_timer,  DIFS_timer,  currentMS,  backoffreset, remaining_time] =  update_timer(MS_collision, uplink_windowN, backoff_timer,  DIFS_timer,  currentMS,  backoffreset, uplink_buffer, remaining_time,MS_IN, Ith, DIFSth,0, max_time,CWminN, CWmaxN, slot_time);
        %PCF uplink
        elseif mod(simulation_time, 4) == 2
            time = 1;
            while time <= PCF_UL_time
                start_MS = 1;
                index = 1 : length(MS);
                while start_MS + NumFDMA-1 <= MSperCell
                    % polled_MS: stored the MS to be transmitted, dim1: # of
                    % frequencies, dim2: index of the MS of that frequency
                    polled_MS = zeros(NumFDMA*frequency_reuse_factor,length(AP));
                    Num_polled_MS = ones(NumFDMA*frequency_reuse_factor,1);
                    for i = 1 : length(AP)
                        MS_under_AP = index(APofMS == i);
                        for j = start_MS : start_MS + NumFDMA-1
                            if j<= length(MS_under_AP)
                                freq = MS_frequency(MS_under_AP(j));
                                polled_MS(freq,Num_polled_MS(freq)) = MS_under_AP(j);
                                Num_polled_MS(freq) = Num_polled_MS(freq)+1;
                            end
                        end
                    end
                    start_MS = start_MS + NumFDMA;
                    for i = 1: NumFDMA*frequency_reuse_factor
                        now = false(1,length(MS));
                        for j = 1: sum(polled_MS(i,:)~=0)
                            if randi([1,100])<=ul_prob_thres
                                now(polled_MS(i,j)) = true;
                            end
                            %now(polled_MS(i,j)) = true;
                        end
                        poll_SINR = calculate_AP_SINR(MS(now), AP(AP_frequency == floor((i-1)/NumFDMA)+1), AP(APofMS(now)),shift,side);
                        poll_capacity = bandwidth/(NumFDMA*frequency_reuse_factor)*log2(1 + poll_SINR);
                        poll_capacity = iteration_time*poll_capacity;
                        % transmitted_data: unit is # of packet
                        for j =1 : length(poll_SINR)
                            % check remaining polled packet
                            transmitted_size = poll_capacity(j)/payload_size; % # of transmitted payload unit time
                            if remaining_poll_UL_packet(polled_MS(i,j))>0
                                if poll_SINR(j)<sinr_thres
                                    if counted_UL_poll_loss(polled_MS(i,j))==0
                                        poll_UL_packet_loss = poll_UL_packet_loss + 1;
                                        counted_UL_poll_loss(polled_MS(i,j))=1;
                                    end                                
                                end
                                if transmitted_size < remaining_poll_DL_packet(polled_MS(i,j))
                                    remaining_poll_UL_packet(polled_MS(i,j)) = remaining_poll_UL_packet(polled_MS(i,j)) - transmitted_size;
                                    continue
                                else
                                    transmitted_size = transmitted_size - remaining_poll_UL_packet(polled_MS(i,j));
                                    counted_UL_poll_loss(i)=0;
                                end
                            end    
                            transmitted_data = min(transmitted_size,uplink_buffer(polled_MS(i,j)));
                            remaining_poll_UL_packet(polled_MS(i,j)) = ceil(transmitted_data) - transmitted_data;
                            uplink_buffer(polled_MS(i,j)) = uplink_buffer(polled_MS(i,j)) - ceil(transmitted_data);
                            uplink_total_packet = uplink_total_packet - ceil(transmitted_data);
                            poll_total_UL_packet = poll_total_UL_packet + ceil(transmitted_data);
                            if poll_SINR(j)<sinr_thres
                                poll_UL_packet_loss = poll_UL_packet_loss + ceil(transmitted_data);
                                counted_UL_poll_loss(polled_MS(i,j))=1;
                            end
                        end
                    end  
                end
                time = time + 1;
            end
        %PCF downlink
        else
            time = 1;
            while time <= PCF_DL_time
                start_MS = 1;
                index = 1 : length(MS);
                while start_MS + NumFDMA-1 <= MSperCell
                    % polled_MS: stored the MS to be transmitted, dim1: # of
                    % frequencies, dim2: index of the MS of that frequency
                    polled_MS = zeros(NumFDMA*frequency_reuse_factor,length(AP));
                    Num_polled_MS = ones(NumFDMA*frequency_reuse_factor,1);
                    for i = 1 : length(AP)
                        MS_under_AP = index(APofMS == i);
                        for j = start_MS : start_MS + NumFDMA-1
                            if j<= length(MS_under_AP)
                                freq = MS_frequency(MS_under_AP(j));
                                polled_MS(freq,Num_polled_MS(freq)) = MS_under_AP(j);
                                Num_polled_MS(freq) = Num_polled_MS(freq)+1;
                            end
                        end
                    end
                    start_MS = start_MS + NumFDMA;
                    for i = 1: NumFDMA*frequency_reuse_factor
                        now = false(1,length(MS));
                        for j = 1: sum(polled_MS(i,:)~=0)
                            if randi([1,100])<=dl_prob_thres
                                now(polled_MS(i,j)) = true;
                            end
                        end
                        if sum(now)==0
                            continue
                        end
                        poll_SINR = calculate_MS_SINR(MS(now), AP(APofMS(now)), AP(APofMS(now)),shift,side);
                        poll_capacity = bandwidth/(NumFDMA*frequency_reuse_factor)*log2(1 + poll_SINR);
                        poll_capacity = iteration_time*poll_capacity;
                        % transmitted_data: unit is # of packet
                        for j =1 : length(poll_SINR)
                            % check remaining polled packet
                            transmitted_size = poll_capacity(j)/payload_size; % # of transmitted payload unit time
                            if remaining_poll_DL_packet(polled_MS(i,j))>0
                                if poll_SINR(j)<sinr_thres
                                    if counted_DL_poll_loss(polled_MS(i,j))==0
                                        poll_DL_packet_loss = poll_DL_packet_loss + 1;

                                        counted_DL_poll_loss(polled_MS(i,j))=1;
                                    end                                
                                end
                                if transmitted_size < remaining_poll_DL_packet(polled_MS(i,j))
                                    remaining_poll_DL_packet(polled_MS(i,j)) = remaining_poll_DL_packet(polled_MS(i,j)) - transmitted_size;
                                    continue
                                else
                                    transmitted_size = transmitted_size - remaining_poll_DL_packet(polled_MS(i,j));
                                    counted_DL_poll_loss(i)=0;
                                end
                            end    
                            transmitted_data = min(transmitted_size,downlink_buffer(polled_MS(i,j)));
                            remaining_poll_DL_packet(polled_MS(i,j)) = ceil(transmitted_data) - transmitted_data;
                            downlink_buffer(polled_MS(i,j)) = downlink_buffer(polled_MS(i,j)) - ceil(transmitted_data);
                            downlink_total_packet = downlink_total_packet - ceil(transmitted_data);
                            poll_total_DL_packet = poll_total_DL_packet + ceil(transmitted_data);
                            if poll_SINR(j)<sinr_thres
                                poll_DL_packet_loss = poll_DL_packet_loss + ceil(transmitted_data);
                                counted_DL_poll_loss(polled_MS(i,j))=1;
                            end
                        end
                    end  
                end
                time = time + 1;
            end
        end

        downlink_curAP = downlink_curAP+sum(currentAP);
        uplink_curMS = uplink_curMS+length(find(currentMS));
    end
    packet_loss(1, lm) = downlink_packet_loss / (downlink_total_packet - sum(downlink_buffer)) * 100;
    packet_loss(2, lm) = uplink_packet_loss / (uplink_total_packet - sum(uplink_buffer)) * 100;
    packet_loss(3, lm) = poll_DL_packet_loss / poll_total_DL_packet * 100;
    packet_loss(4, lm) = poll_UL_packet_loss / poll_total_UL_packet * 100;
    throughput(1, lm) = (downlink_total_packet-sum(downlink_buffer))* payload_size / iteration_time / simulation_time;
    throughput(2, lm) = (uplink_total_packet-sum(uplink_buffer))* payload_size / iteration_time / simulation_time;
    throughput(3, lm) = poll_total_DL_packet* payload_size / iteration_time / simulation_time;
    throughput(4, lm) = poll_total_UL_packet* payload_size / iteration_time / simulation_time;
    buffutil(1, lm) = sum(downlink_buffer) * payload_size/ downlink_buffer_size / length(downlink_buffer) * 100;
    buffutil(2, lm) = sum(uplink_buffer) * payload_size/ uplink_buffer_size / length(uplink_buffer) * 100;

end

%clearvars temp1 temp2 temp3
disp("downlink")
report_performance(downlink_total_packet, downlink_buffer, downlink_buffer_size, payload_size, downlink_packet_loss, downlink_curAP/simulation_time, 1);
disp("uplink")
report_performance(uplink_total_packet, uplink_buffer, uplink_buffer_size, payload_size, uplink_packet_loss, uplink_curMS/simulation_time,1);
disp("poll uplink")
report_performance(poll_total_UL_packet, 0, 0, payload_size, poll_UL_packet_loss, 0, 0);
disp("poll downlink")
report_performance(poll_total_DL_packet, 0, 0, payload_size, poll_DL_packet_loss, 0, 0);

%plot histogram
figure(2);
plot(mscoll, packet_loss(1, :), '--*r');
hold on;
plot(mscoll, packet_loss(2, :), '--*g');
hold on;
plot(mscoll, packet_loss(3, :), '--*b');
hold on;
plot(mscoll, packet_loss(4, :), '--*c');
hold on;
title('packet loss rate to number of MS per cell');
legend('downlink DCF', 'uplink DCF', 'downlink PCF', 'uplink PCF');
xlabel('number of MS per cell');
ylabel('packet loss rate (%)');
hold off;
%saveas(gcf,'ms_loss.jpg');
figure(3);
plot(mscoll, throughput(1, :), '--*r');
hold on;
plot(mscoll, throughput(2, :), '--*g');
hold on;
plot(mscoll, throughput(3, :), '--*b');
hold on;
plot(mscoll, throughput(4, :), '--*c');
hold on;
title('throughput to number of MS per cell');
legend('downlink DCF', 'uplink DCF', 'downlink PCF', 'uplink PCF');
xlabel('number of MS per cell');
ylabel('throughput (bit/s)');
hold off;
%saveas(gcf,'ms_throughput.jpg');
figure(4);
plot(mscoll, buffutil(1, :), '--*r');
hold on;
plot(mscoll, buffutil(2, :), '--*g');
hold on;
title('buffer utilization rate to number of MS per cell');
legend('downlink', 'uplink');
xlabel('number of MS per cell');
ylabel('buffer utilization rate (%)');
hold off;
%saveas(gcf,'ms_buffer.jpg');