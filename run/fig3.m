% Figure 3 - manifolds

color2 = [0 0 .5];
color1 = [.3 .5 .8];
gray2 = [.3 .3 .3];
gray1 = [.8 .8 .8];

%% Oval manifold

s = linspace(0, 2*pi, 1001);
Tx = 2*cos(s) + 3;
Ty = 2*cos(s+pi/4) + 3;

A = [1 0; 0 1]';
x1 = [1 4]';
Ax = A*x1;

subplot(2, 3, [1 4]);
hold on;
plot(Tx, Ty, 'linewidth', 2, 'color', [.5 0 .5]);
plot(0, 0);
for i = 1:100:length(s)
    plot([Ax(1) Tx(i)], [Ax(2) Ty(i)], 'Color', gray2);
end
plot(Ax(1), Ax(2), 'o', 'Color', color2, 'MarkerFaceColor', color2);
axis equal; axis square; axis off;
set(gca, 'xtick', [], 'ytick', []);

%% Squiggle manifold

s = linspace(0, 2*pi, 1001);
s2 = pi/2 + 4*pi ./ (1 + exp(-linspace(-10, 10, length(s))));
Tx = 3*cos(s) - sin(s2) + 3;
Ty = 3*cos(s+pi/3) - sin(s2+pi/4) + 3;

A = [1 0; 0 1]';
x1 = [4.8 4.5]';
Ax1 = A*x1;
x2 = [2 3]';
Ax2 = A*x2;

subplot(2, 3, [2 5]);
hold on;
plot(Tx, Ty, 'linewidth', 2, 'color', [.5 0 .5]);
plot(0, 0);
for i = 1:100:length(s)
    plot([Ax1(1) Tx(i)], [Ax1(2) Ty(i)], 'Color', gray1);
    plot([Ax2(1) Tx(i)], [Ax2(2) Ty(i)], 'Color', gray2);
end
plot(Ax1(1), Ax1(2), 'o', 'Color', color1, 'MarkerFaceColor', color1);
plot(Ax2(1), Ax2(2), 'o', 'Color', color2, 'MarkerFaceColor', color2);
axis equal; axis square; axis off;
set(gca, 'xtick', [], 'ytick', []);

dx1 = Ax1(1) - Tx;
dy1 = Ax1(2) - Ty;
pdf1 = exp(-1/2*(dx1.^2 + dy1.^2));
pdf1 = pdf1 / sum(pdf1(:));

subplot(2, 3, 3);
plot(s, pdf1, 'color', color1);
set(gca, 'xtick', [], 'ytick', 0);
ylabel('p(s|x)');
xlim([-inf inf]);

dx2 = Ax2(1) - Tx;
dy2 = Ax2(2) - Ty;
pdf2 = exp(-1/2*(dx2.^2 + dy2.^2));
pdf2 = pdf2 / sum(pdf2(:));

subplot(2, 3, 6);
plot(s, pdf2, 'color', color2);
set(gca, 'xtick', [], 'ytick', 0);
ylabel('p(s|x)');
xlabel('s');
xlim([-inf inf]);