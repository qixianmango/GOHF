function new_direction = perturb_direction(original_direction, max_offset_angle)

    sigma = 1/3 * atan(max_offset_angle);
    original_direction = original_direction / norm(original_direction);
    if original_direction(1) ~= 0 || original_direction(2) ~= 0
        temp_vector = [0; 0; 1];
    else
        temp_vector = [1; 0; 0];
    end

    basis1 = cross(original_direction, temp_vector);
    basis1 = basis1 / norm(basis1);

    basis2 = cross(original_direction, basis1);
    basis2 = basis2 / norm(basis2);

    rand_point = sigma * randn(2, 1);

    point = (basis1 * rand_point(1, :) + basis2 * rand_point(2, :))';
    new_direction = original_direction + point';
    new_direction = new_direction / norm(new_direction);
end