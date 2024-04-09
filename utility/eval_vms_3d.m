% Evaluate the three-dimensional von Mises stress
function vms = eval_vms_3d(sts)
    tmp = 0.5 * ((sts(1, 1, :, :) - sts(2, 2, :, :)).^2 + (sts(2, 2, :, :) - sts(3, 3, :, :)).^2 + (sts(3, 3, :, :) - sts(1, 1, :, :)).^2) + ...
          3.0 * (sts(1, 2, :, :) .* sts(2, 1, :, :) + sts(2, 3, :, :) .* sts(3, 2, :, :) + sts(3, 1, :, :) .* sts(1, 3, :, :));
    vms = squeeze(sqrt(max(tmp, 0)));
end
