function [corner1, corner2, corner3, corner4] = define_mirror_corners(target_mirror, mirror_height, mirror_width, mirror_normal)

    up_vector = [0, 0, 1];
    side_vector1 = cross(mirror_normal, up_vector);
    side_vector1 = side_vector1 / norm(side_vector1);
    side_vector2 = cross(mirror_normal, side_vector1);
    side_vector2 = side_vector2 / norm(side_vector2);

    corner1 = target_mirror - (mirror_width / 2) * side_vector1 - (mirror_height / 2) * side_vector2;
    corner2 = target_mirror + (mirror_width / 2) * side_vector1 - (mirror_height / 2) * side_vector2;
    corner3 = target_mirror + (mirror_width / 2) * side_vector1 + (mirror_height / 2) * side_vector2;
    corner4 = target_mirror - (mirror_width / 2) * side_vector1 + (mirror_height / 2) * side_vector2;
end