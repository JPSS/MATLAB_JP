function [ gamma ] = estimate_gamma_bands(dd, da, aa, ref)
%ESTIMATE_GAMMA Summary of this function goes here
%   Detailed explanation goes here

%% Select areas of d-only and d+a band
display('First, select donor-only band. Then D+A band.')
[I pos]  = integrate_areas({dd, da, ref}, 2, 1); % colum 1 = dd, column 2 = da, column 3 = ref

gamma_app  = (I(2,2)-I(1,2) ) / (I(1,1)-I(2,1)) % approximated gamma factor without reference

ratio = I(2,3) / I(1,3)
gamma = (I(2,2)/I(2,3)-I(1,2)/I(1,3) ) / (I(1,1)/I(1,3)-I(2,1)/I(2,3)) % approximated gamma factor, normalizing with gamma-factor
gamma = (I(2,2)-ratio*I(1,2) ) / (ratio*I(1,1)-I(2,1)) % approximated gamma factor, normalizing with gamma-factor


%% plot
close all
subplot(1,2,1)
h = bar(I);
set(h(1),'FaceColor','g')
set(h(2),'FaceColor','b')
set(h(3),'FaceColor','k')
ylabel('Integrated Intensity [a.u.]')
legend({'D->D', 'D->A', 'Reference'})
set(gca, 'XTickLabel', {'D-only', 'D+A'})

subplot(1,2,2)
h = bar(I');
set(h(1),'FaceColor','k')
set(h(2),'FaceColor','k')
%ylabel('Integrated Intensity [a.u.]')
legend({'D-only', 'D+A'})
set(gca, 'XTickLabel', {'D->D', 'D->A', 'Reference'})
pause(2)
close all
%%



end

