function [shadow_eff, trunc_eff, cosine_eff, atm_trans] = evaluate_efficiency(mirror_pos, mirror_index, sunlight_dir, mirror_normal, tower_loc, mirror_corners)

    num_points = 20;
    nearby_distance = 20;
    num_epochs = 1;

    target_mirror = mirror_pos(mirror_index, 1:3);
    distances = sqrt(sum((target_mirror - mirror_pos).^2, 2));
    nearby_indices = distances <= nearby_distance & distances > 0;
    nearby_mirror_corners = mirror_corners(:, :, nearby_indices);
    points = sample_points(mirror_corners(:, :, mirror_index), num_points);
    counter_not_blocked = 0;
    counter_reach_collector = 0;

    for e = 1:num_epochs
        for i = 1:size(points, 1)
            blocked = false;
            new_direction = perturb_direction(sunlight_dir, 4.65e-3);
            reflection_direction = reflect_ray(-sunlight_dir, mirror_normal);
            blocked = detect_cylinder_overlap(points(i, :), new_direction, [tower_loc(1), tower_loc(2), 0], 84, 3.5);
            for j = 1:size(nearby_mirror_corners, 3)
                if blocked
                    break;
                end
                blocked = blocked | detect_rectangle_overlap(points(i, :), new_direction, nearby_mirror_corners(:, :, j));
                blocked = blocked | detect_rectangle_overlap(points(i, :), reflection_direction, nearby_mirror_corners(:, :, j));
            end
            if ~blocked
                counter_not_blocked = counter_not_blocked + 1;
                if detect_cylinder_overlap(points(i, :), reflection_direction, [tower_loc(1), tower_loc(2), 0], 84, 3.5)
                    counter_reach_collector = counter_reach_collector + 1;
                end
            end
        end
    end

    shadow_eff = counter_not_blocked / (num_epochs * size(points, 1));

    if counter_not_blocked == 0
        trunc_eff = 0;
    else
        trunc_eff = counter_reach_collector / counter_not_blocked;
    end

    cosine_eff = cos(0.5 * acos(dot(real(sunlight_dir), real(reflection_direction))));
    assert(isreal(cosine_eff));
    d_HR = norm(tower_loc - target_mirror);
    atm_trans = 0.99321 - 0.0001176 * d_HR + 1.97e-8 * d_HR^2;
end