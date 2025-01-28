function [isIntersected] = detect_cylinder_overlap(ray_start, ray_dir, cylinder_center, cylinder_height, cylinder_radius)
    isIntersected = false;
    intersection_point = [NaN, NaN, NaN];

    t1 = (cylinder_center(3) - ray_start(3)) / ray_dir(3);
    t2 = (cylinder_center(3) + cylinder_height - ray_start(3)) / ray_dir(3);

    p1 = ray_start + t1 .* ray_dir;
    p2 = ray_start + t2 .* ray_dir;

    if t1 > 0 && (norm(p1(1:2) - cylinder_center(1:2)) <= cylinder_radius)
        isIntersected = false;
        intersection_point = [];
        return;
    elseif t2 > 0 && (norm(p2(1:2) - cylinder_center(1:2)) <= cylinder_radius)
        isIntersected = true;
        intersection_point = p2;
        return;
    end

    d = ray_start(1:2) - cylinder_center(1:2);
    A = ray_dir(1)^2 + ray_dir(2)^2;
    B = 2 * (ray_dir(1)*d(1) + ray_dir(2)*d(2));
    C = d(1)^2 + d(2)^2 - cylinder_radius^2;

    discriminant = B^2 - 4*A*C;

    if discriminant >= 0
        t3 = (-B + sqrt(discriminant)) / (2*A);
        t4 = (-B - sqrt(discriminant)) / (2*A);

        p3 = ray_start + t3 .* ray_dir;
        p4 = ray_start + t4 .* ray_dir;

        if t3 > 0 && (cylinder_center(3) <= p3(3) && p3(3) <= cylinder_center(3) + cylinder_height)
            isIntersected = true;
            intersection_point = p3;
        elseif t4 > 0 && (cylinder_center(3) <= p4(3) && p4(3) <= cylinder_center(3) + cylinder_height)
            isIntersected = true;
            intersection_point = p4;
        end
    end
end