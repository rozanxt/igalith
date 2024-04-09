% Evaluate a raster tile map using nearest-neighbor interpolation
% and the 'uv' axes convention.
%
% The coordinate system origin is at the lower left corner.
% The u-axis is horizontal and is numbered from left to right.
% The v-axis is vertical and is numbered from bottom to top.
function e = eval_uv(data, u, v, tile_size)
    [m, n] = size(data);
    [w, h] = deal(tile_size(1), tile_size(2));
    a = min(max(round(m - h * v + 0.5), 1), m);
    b = min(max(round(w * u + 0.5), 1), n);
    e = data(sub2ind([m, n], a, b));
end
