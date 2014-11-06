%%

s1 =   [10 11.78 3.34 0.20;
        20 22.88 4.54 0.26;
        30 32.46 5.58 0.32]; % old sequences
    
s2 =   [30 36.56 5.15 0.30;
        40 50.20 7.29 0.42;
        50 60.35 5.74 0.33]; % orthogonal sequences
    
v2data_new = [30 24.98 5.60 0.40];
%%

structure = cell(0,4);
folder = '/Users/Administrator/Documents/Jonas/FRET_STAGE/';

structure = [structure; get_structure([folder 'Results_8bp_v2.txt'] , '8bp Spacer',8)];
structure = [structure; get_structure([folder 'Results_10bp_v2.txt'] , '10bp Spacer',10)];
structure = [structure; get_structure([folder 'Results_12bp_v2.txt'],'12bp Spacer', 12)];
structure = [structure; get_structure([folder 'Results_14bp_v2.txt'] , '14bp Spacer',14)];
structure = [structure; get_structure([folder 'Results_16bp_v2.txt'] , '16bp Spacer',16)];
structure = [structure; get_structure([folder 'Results_18bp_v2.txt'] , '18bp Spacer',18)];
structure = [structure; get_structure( [folder 'Results_20bp_v2.txt'], '20bp Spacer',20)];
structure = [structure; get_structure( [folder 'Results_22bp_v2.txt'], '22bp Spacer',22)];
structure = [structure; get_structure( [folder 'Results_24bp_v2.txt'], '24bp Spacer',24)];
structure = [structure; get_structure( [folder 'Results_30bp_v2.txt'], '30bp Spacer',30)];

v2data = zeros(size(structure,1), 4);
for i=1:size(structure,1)
    v2data(i,:) = [structure{i,3} mean(structure{i,1})  std(structure{i,1}) std(structure{i,1})/sqrt(length(structure{i,1})) ];
    
end
%%
close all
plot(s1(:,1), s1(:,2), 'b.-', s2(:,1), s2(:,2), 'r.-')


set(gca, 'XLim', [0 70])
set(gca, 'YLim', [0 70])

xlabel('Contour length [bp]')
ylabel('Angle [deg]')

%axis square
%%
xplot = 0:1:70;

close all
fig_dim =[20 25];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

subplot(2, 1, 1)
errorbar(s1(:,1), s1(:,2), s1(:,3), 'b.-'), hold on
errorbar(s2(:,1), s2(:,2), s2(:,3), 'r.-')
errorbar(v2data(:,1), v2data(:,2), v2data(:,3), 'k.-' )
errorbar(v2data_new(:,1), v2data_new(:,2), v2data_new(:,3), 'c.-')
plot(xplot, 2*asin(xplot/2/47.5)*180/pi, 'g')
legend({'ssHinge - old seq', 'ssHinge - orthogonal seq', 'dsHinge - old seq',  'dsHinge - old seq - new measured','Expected'}, 'location', 'southeast')

set(gca, 'XLim', [0 70])
set(gca, 'YLim', [0 70])

xlabel('Contour length [bp]')
ylabel('Angle [deg]')


%
subplot(2, 1, 2)
errorbar(s1(:,1), sin(s1(:,2)*pi/180/2), sin(s1(:,3)*pi/180/2), 'b.-'), hold on
errorbar(s2(:,1), sin(s2(:,2)*pi/180/2), sin(s2(:,3)*pi/180/2), 'r.-')
errorbar(v2data(:,1), sin(v2data(:,2)*pi/180/2), sin(v2data(:,3)*pi/180/2), 'k.-' )
errorbar(v2data_new(:,1), sin(v2data_new(:,2)*pi/180/2), sin(v2data_new(:,3)*pi/180/2), 'c.-')



c_s1 = polyfit(s1(:,1) , sin(s1(:,2)*pi/180/2), 1);
c_s2 = polyfit(s2(:,1) , sin(s2(:,2)*pi/180/2), 1);
c_v2 = polyfit(v2data(:,1) , sin(v2data(:,2)*pi/180/2), 1);
plot(xplot, xplot/2/47.5, 'g', xplot, c_s1(2)+xplot*c_s1(1), 'b--', xplot, c_s2(2)+xplot*c_s2(1), 'r--', xplot, c_v2(2)+xplot*c_v2(1), 'k--')



legend({'ssHinge - old seq', 'ssHinge - orthogonal seq', 'dsHinge - old seq',  'dsHinge - old seq - new measured','Expected'}, 'location', 'southeast')

set(gca, 'XLim', [0 70])
set(gca, 'YLim', sin([0 70]*pi/180/2))

xlabel('Contour length [bp]')
ylabel('sin(alpha/2)')

%

print(cur_fig, '-dtiff', '-r600' , [data_dir filesep 'Jonas' filesep '2013-10-09_angle_distributions_02.tif'])
print(cur_fig, '-depsc2', [data_dir filesep 'Jonas' filesep '2013-10-09_angle_distributions_02.eps'])





