clear;
clc;

tower_loc = [0, 0, 80];
min_radius = 100;
install_height_factor = 0;
height_factor = 0;
width_factor = 0;
a1 = 2.875;
a0 = 5.75;
init_circle_num = 75;
spacing = 5.75 + 5;

[~, ~, radii, ~] = arrange_concentric_circles(tower_loc, init_circle_num, min_radius, spacing);
[valid, install_heights, heights, widths, spacing, radii] = adjust_mirror_params(install_height_factor, height_factor, width_factor, radii, tower_loc, a0, a1);
assert(valid)

[total_energy, unit_energy, avg_efficiency] = evaluate_performance(install_heights, heights, widths, init_circle_num, min_radius, spacing, tower_loc);

flag = true;
epoch = 0;
while flag
    epoch = epoch + 1;
    if epoch ~= 1
        step = 2.5 + 15 * rand();
        step_sizes = [-step, 0, step];

        neighbor_points = [];
        for dx = step_sizes
            for dy = step_sizes
                new_x = tower_loc(1) + dx;
                new_y = tower_loc(2) + dy;
                neighbor_points = [neighbor_points; [new_x, new_y]];
            end
        end

        unit_energy_up_layer = zeros(1, length(step_sizes)^2);
        energy_up_layer = zeros(1, length(step_sizes)^2);
        eta_up_layer = zeros(1, length(step_sizes)^2);
        for i = 1:length(step_sizes)^2
            [~, install_heights, heights, widths, spacing] = adjust_mirror_params(install_height_factor, height_factor, width_factor, radii, [neighbor_points(i, :), 80], a0, a1);
            [~, data_len] = arrange_concentric_circles([neighbor_points(i, :), 80], init_circle_num, min_radius, spacing);

            [energy_up_layer(i), unit_energy_up_layer(i), eta_up_layer(i)] = evaluate_performance(install_heights, heights, widths, init_circle_num, min_radius, spacing, [neighbor_points(i, :), 80]);
                    end

        [~, index] = max(eta_up_layer);
        tower_loc = [neighbor_points(index, :), 80];
            end

    step = 0.00025;
    range = 0.001;
    for inner_loop = 1:2
        best_install_height_factor = install_height_factor;
        best_unit_energy = -inf;
        best_radii = radii;
        for new_install_height_factor = install_height_factor-range:step:install_height_factor+range
            [valid, install_heights, heights, widths, spacing, tmp_radii] = adjust_mirror_params(new_install_height_factor, height_factor, width_factor, radii, tower_loc, a0, a1);
            if ~valid
                continue;
            end
            [res_energy, res_unit_energy, ~] = evaluate_performance(install_heights, heights, widths, init_circle_num, min_radius, spacing, tower_loc);
            if res_energy > 60 && res_unit_energy > best_unit_energy
                best_unit_energy = res_unit_energy;
                best_install_height_factor = new_install_height_factor;
                best_radii = tmp_radii;
            end
            end

        install_height_factor = best_install_height_factor;
        radii = best_radii;
        
    end

    if install_height_factor == best_install_height_factor && height_factor == best_height_factor && width_factor == best_width_factor
        flag = false;
    end

    [valid, install_heights, heights, widths, spacing, radii] = adjust_mirror_params(install_height_factor, height_factor, width_factor, radii, tower_loc, a0, a1);
    [total_energy, unit_energy, avg_efficiency] = evaluate_performance(install_heights, heights, widths, init_circle_num, min_radius, spacing, tower_loc);
    
end