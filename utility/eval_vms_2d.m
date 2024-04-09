% Evaluate the two-dimensional von Mises stress
function vms = eval_vms_2d(sts)
    tmp = sts(1, 1, :, :).^2 + sts(2, 2, :, :).^2 - ...
          sts(1, 1, :, :) .* sts(2, 2, :, :) + ...
          3 * sts(1, 2, :, :) .* sts(2, 1, :, :);
    vms = squeeze(sqrt(max(tmp, 0)));
end
