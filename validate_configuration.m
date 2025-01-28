function [valid, flag] = validate_configuration(install_height, height, width, radii, point_num)

    if width < height || height > 2 * install_height
        valid = false;
        flag = 1;
        return;
    end
    if width < 2 || width > 8 || height < 2 || height > 8 || install_height < 2 || install_height > 6
        valid = false;
        flag = 2;
        return;
    end
    if ~all(radii >= 100)
        valid = false;
        flag = 3;
        return;
    end
    for i = 2:length(radii)
        if radii(i) - radii(i-1) < width + 5
            valid = false;
            flag = 4;
            return;
        end
    end
    for i = 1:length(radii)
        if 2 * sin(pi / point_num(i)) * radii(i) < width + 5
            valid = false;
            flag = 100 + i;
            return;
        end
    end
    valid = true;
    flag = 0;
end