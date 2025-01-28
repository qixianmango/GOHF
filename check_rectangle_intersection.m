function [isIntersect, intersectionPoint] = detect_rectangle_overlap(ray_start, ray_dir, rect_corners)

    edge1 = rect_corners(2, :) - rect_corners(1, :);
    edge2 = rect_corners(3, :) - rect_corners(1, :);
    normal = cross(edge1, edge2);
    denom = dot(normal, ray_dir);
    t = dot(rect_corners(1,:) - ray_start, normal) / denom;

    if t < 0
        isIntersect = false;
        intersectionPoint = [];
        return;
    end

    intersectionPoint = ray_start + t * ray_dir;

    AB = rect_corners(2,:) - rect_corners(1,:);
    AM = intersectionPoint - rect_corners(1,:);
    BC = rect_corners(3,:) - rect_corners(2,:);
    BM = intersectionPoint - rect_corners(2,:);

    if (dot(AB, AM) >= 0) && (dot(AB, AM) <= dot(AB, AB)) && ...
       (dot(BC, BM) >= 0) && (dot(BC, BM) <= dot(BC, BC))
        isIntersect = true;
    else
        isIntersect = false;
        intersectionPoint = [];
    end
end