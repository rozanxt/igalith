% Parametrization of a single quadrilateral patch
function nrb = nrb_quad_tile(vertices)
    pts = zeros(4, 2, 2);
    pts(:, 1, 1) = [vertices(1, :), 0, 1].';
    pts(:, 2, 1) = [vertices(2, :), 0, 1].';
    pts(:, 1, 2) = [vertices(3, :), 0, 1].';
    pts(:, 2, 2) = [vertices(4, :), 0, 1].';
    knt = {[0, 0, 1, 1], [0, 0, 1, 1]};
    nrb = nrbmak(pts, knt);
end
