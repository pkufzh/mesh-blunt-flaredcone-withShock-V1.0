
%% Program initization
clear all
clf
close all
clc

%% Data introduction
G = load("grid_xy_plot.dat");
nx_tot = G(1, 1);
ny_tot = G(1, 2);
xx = G((2 : end), 1);
yy = G((2 : end), 2);

A = load("xa.dat");
B = load("xb.dat");
C = load("xc.dat");

nx_part = 4;
nx_buff = 0;
n_eff = length(A) - nx_buff;

id = A((1 : n_eff), 2);

xa = A((1 : n_eff), 3);
ya = A((1 : n_eff), 4);
xb = B((1 : n_eff), 3);
yb = B((1 : n_eff), 4);
xc = C((1 : n_eff), 3);
yc = C((1 : n_eff), 4);

figure(77);
clf
hold on

% plot(xa(id == 1), ya(id == 1), '-', 'LineWidth', 2.0, 'Color', [0.00, 0.00, 0.00], 'HandleVisibility', 'off');
% plot(xa(id == 2), ya(id == 2), '-', 'LineWidth', 2.0, 'Color', [0.00, 0.45, 0.74], 'HandleVisibility', 'off');
% plot(xa(id == 3), ya(id == 3), '-', 'LineWidth', 2.0, 'Color', [0.85, 0.33, 0.10], 'HandleVisibility', 'off');

% for i = 1 : n_eff - 1
%     if (id(i) == 1)
%         plot(xa(), '-', 'LineWidth', 2.0, 'Color', [0.00, 0.00, 0.00], 'HandleVisibility', 'off');
%     elseif (id(i) == 2)
%         plot([xa(i), xa(i + 1)], [ya(i), ya(i + 1)], '-', 'LineWidth', 2.0, 'Color', [0.00, 0.45, 0.74], 'HandleVisibility', 'off');
%     elseif (id(i) == 3)
%         plot([xa(i), xa(i + 1)], [ya(i), ya(i + 1)], '-', 'LineWidth', 2.0, 'Color', [0.85, 0.33, 0.10], 'HandleVisibility', 'off');
%     end
% 
% end

plot(xa, ya, 'm-', 'LineWidth', 2.0);
plot(xb, yb, 'b-', 'LineWidth', 2.0);
plot(xc, yc, 'r-', 'LineWidth', 2.0);

for k = 1 : nx_part
    id_part = find(id == k);
    plot(xa(id_part(1)), ya(id_part(1)), 'k.', 'MarkerSize', 16, 'HandleVisibility', 'off');
    plot([xa(id_part(1)), xc(id_part(1))], [ya(id_part(1)), yc(id_part(1))], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([xc(id_part(1)), xb(id_part(1))], [yc(id_part(1)), yb(id_part(1))], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end
% for k = 1 : 2
%     id_part = find(id == k);
%     plot(xa(id_part(end)), ya(id_part(end)), 'k.', 'MarkerSize', 16, 'HandleVisibility', 'off');
%     plot([xa(id_part(end)), xc(id_part(end))], [ya(id_part(end)), yc(id_part(end))], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
%     plot([xc(id_part(end)), xb(id_part(end))], [yc(id_part(end)), yb(id_part(end))], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% end

for i = 1 : n_eff
    plot([xa(i), xc(i)], [ya(i), yc(i)], '-', 'Linewidth', 0.5, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');
    plot([xc(i), xb(i)], [yc(i), yb(i)], '-', 'Linewidth', 0.5, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');
end
hold off
axis equal
legend('wall surface', 'farfield', 'shock wave');


figure(88);

hold on
for j = 1 : ny_tot

    ks = (j - 1) * nx_tot + 1;
    kt = (j - 1) * nx_tot + nx_tot;
    plot(xx(ks : kt), yy(ks : kt), '-', 'Linewidth', 0.5, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');

%     for i = 1 : nx_tot - 1
%         k = i + (j - 1) * nx_tot;
%         plot([xx(k), xx(k + 1)], [yy(k), yy(k + 1)], '-', 'Linewidth', 0.5, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');
%     end

end

plot(xa, ya, 'm-', 'LineWidth', 2.0);
plot(xb, yb, 'b-', 'LineWidth', 2.0);
plot(xc, yc, 'r-', 'LineWidth', 2.0);

for k = 1 : nx_part
    id_part = find(id == k);
    plot(xa(id_part(1)), ya(id_part(1)), 'k.', 'MarkerSize', 16, 'HandleVisibility', 'off');
    plot([xa(id_part(1)), xc(id_part(1))], [ya(id_part(1)), yc(id_part(1))], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot([xc(id_part(1)), xb(id_part(1))], [yc(id_part(1)), yb(id_part(1))], 'k-', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

for i = 1 : n_eff
    plot([xa(i), xc(i)], [ya(i), yc(i)], '-', 'Linewidth', 0.5, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');
    plot([xc(i), xb(i)], [yc(i), yb(i)], '-', 'Linewidth', 0.5, 'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');
end
hold off
axis equal
legend('wall surface', 'farfield', 'shock wave');

hold off
axis equal

