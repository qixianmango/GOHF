function local_time = convert_to_local_time(beijing_time)
    time_difference = minutes(16);
    local_time = beijing_time - time_difference;
end