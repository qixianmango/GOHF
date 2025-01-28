function [total_energy, unit_energy, avg_efficiency, results_table] = evaluate_performance(install_heights, heights, widths, circle_num, min_radius, spacing, tower_loc)

    mirror_positions = arrange_concentric_circles(tower_loc, circle_num, min_radius, spacing);
    mirror_positions = [mirror_positions, install_heights];

    dates = {'2025-1-21', '2025-2-21', '2025-3-21', '2025-4-21', '2025-5-21', '2025-6-21', '2025-7-21', '2025-8-21', '2025-9-21', '2025-10-21', '2025-11-21', '2025-12-21'};
    times = {'9:00:00', '10:30:00', '12:00:00', '13:30:00', '15:00:00'};
    datetimeMatrix = repmat(datetime('now'), 12, 5);
    for i = 1:length(dates)
        for j = 1:length(times)
            datetimeString = [dates{i} ' ' times{j}];
            datetimeMatrix(i, j) = datetime(datetimeString, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        end
    end

    num_mirrors = size(mirror_positions, 1);
    num_entries = num_mirrors * length(dates) * length(times);
    results_table = table('Size', [num_entries, 5], 'VariableTypes', {'datetime', 'double', 'double', 'double', 'double'}, 'VariableNames', {'DateTime', 'MirrorID', 'OpticalEfficiency', 'TruncationEfficiency', 'OutputPower'});
    idx = 1;

    total_energy = 0;
    total_optical_efficiency = 0;
    for j = 1:length(dates)
        for i = 1:length(times)
            results = compute_mirror_performance(datetimeMatrix(j, i), mirror_positions, heights, widths, tower_loc);
            for m = 1:num_mirrors
                results_table.DateTime(idx) = datetimeMatrix(j, i);
                results_table.MirrorID(idx) = m;
                results_table.OpticalEfficiency(idx) = results(m).optical_efficiency;
                results_table.TruncationEfficiency(idx) = results(m).truncation_efficiency;
                results_table.OutputPower(idx) = results(m).output_power;
                total_energy = total_energy + results(m).output_power;
                total_optical_efficiency = total_optical_efficiency + results(m).optical_efficiency;
                idx = idx + 1;
            end
        end
    end

    total_energy = total_energy / (length(dates) * length(times));
    avg_efficiency = total_optical_efficiency / (num_mirrors * length(dates) * length(times));
    total_area = sum(widths .* heights);
    unit_energy = total_energy / total_area;
end