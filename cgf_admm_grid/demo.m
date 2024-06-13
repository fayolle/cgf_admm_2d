clear; close all; clc;

I  = im2double(imread('../data/mushroom.png'));
%I  = im2double(imread('../data/keyhole.png'));
%I  = im2double(imread('../data/bagel.png'));

% cgf admm
r = 1.0;
niter = 100;
nupdate = 1; % n: Solve a Poisson problem every n iterations 
I_cgfadmm = cgf_admm(I, niter, r, nupdate);
[gx, gy] = gradient(I_cgfadmm);
gr_cgfadmm = sqrt(gx.*gx+gy.*gy);

figure,
imshow([I_cgfadmm./max(I_cgfadmm(:)),gr_cgfadmm]),
title('CGF-ADMM and its gradient');

I_cgfadmm_ud = flipud(I_cgfadmm);
figure,
contourf(I_cgfadmm_ud, 10), title('CGF-ADMM');
colormap jet;

figure,
[C, h] = contourf(I_cgfadmm_ud, 25);
title('CGF-ADMM');
colormap jet;
set(h, 'LineColor', 'none');

% interesting options
%contourf(I_cgfadmm,10); %10 contour levels
%[M,c]=contourf(I_cgfadmm); c.LineWidth=3;
