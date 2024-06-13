function norm = compute_norm(dU)
sqr_norm = dU(:,1).^2 + dU(:,2).^2;
norm = sqrt(sqr_norm);
end
