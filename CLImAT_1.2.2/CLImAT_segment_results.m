function results = CLImAT_segment_results(State_pre,normal_state)
%this function is used to generate abnormal segments from a sequence of assigned
%states
results = [];
start_pos = 0;

% for now_pos=1:length(State_pre)
%     if start_pos == 0
%         start_pos = now_pos;
%     elseif State_pre(start_pos) ~= State_pre(now_pos)
%         %and a new segment starts
%         results = [results; [start_pos now_pos-1 State_pre(start_pos)]];
%         start_pos = now_pos; %current position is the start
%     end
% end
% if start_pos <= now_pos
%     results = [results; [start_pos now_pos State_pre(start_pos)]];
% end

for now_pos = 1:length(State_pre)
    if start_pos == 0 %not in a abnormal segment yet
        if State_pre(now_pos) ~= normal_state % a new abnormal segment starts now
            start_pos = now_pos;
        end
    elseif State_pre(start_pos) ~= State_pre(now_pos) % %in a abnormal segment right now
        %and a new segment starts
        results = [results; [start_pos now_pos-1 State_pre(start_pos)]];
        if State_pre(now_pos) ~= normal_state %a new abnormal segment
            start_pos = now_pos; %current position is the start
        else %the new segment has normal state
            start_pos = 0;
        end
    end

    if (now_pos == length(State_pre)) && (start_pos > 0)
        %The end of this chromosome is a abnormal segment
        results = [results; [start_pos now_pos State_pre(start_pos)]];
    end
end