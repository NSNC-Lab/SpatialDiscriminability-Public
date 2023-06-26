% for comparing figure 3f and supplementary figure s8 to see if decrease in
% performance is significantly different between arch and non-arch subjects

% experimental group (figure 3f)
figure(1);
subplot(1,2,1); a = findobj(gca,'type','line'); a = a(3:end);
for i = 1:length(a), exp_clean(i,:) = a(i).YData; end
subplot(1,2,2); a = findobj(gca,'type','line'); a = a(3:end);
for i = 1:length(a), exp_masked(i,:) = a(i).YData; end

% control group (supplementary figure s8)
figure(2);
subplot(1,2,1); a = findobj(gca,'type','line'); a = a(3:end);
for i = 1:length(a), ctrl_clean(i,:) = a(i).YData; end
subplot(1,2,2); a = findobj(gca,'type','line'); a = a(3:end);
for i = 1:length(a), ctrl_masked(i,:) = a(i).YData; end

[h,p] = ttest2(diff(exp_clean,1,2),diff(ctrl_clean,1,2))
d = computeCohen_d(diff(exp_clean,1,2),diff(ctrl_clean,1,2))

[h,p] = ttest2(diff(exp_masked,1,2),diff(ctrl_masked,1,2))
d = computeCohen_d(diff(exp_masked,1,2),diff(ctrl_masked,1,2))
