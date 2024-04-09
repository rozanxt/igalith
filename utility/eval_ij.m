% Evaluate a raster tile map using nearest-neighbor interpolation
% and the 'ij' axes convention.
%
% The coordinate system origin is at the upper left corner.
% The i-axis is vertical and is numbered from top to bottom.
% The j-axis is horizontal and is numbered from left to right.
function e = eval_ij(data, i, j, tile_size)
    [m, n] = size(data);
    [h, w] = deal(tile_size(1), tile_size(2));
    a = min(max(round(h * i + 0.5), 1), m);
    b = min(max(round(w * j + 0.5), 1), n);
    e = data(sub2ind([m, n], a, b));
end
