function [x_cropped, y_cropped] = delimit_window(x, y, window_start, window_end)
    indexes = find(x >= window_start & x <= window_end);
    x_cropped = x(indexes);
    y_cropped = y(indexes);
end
