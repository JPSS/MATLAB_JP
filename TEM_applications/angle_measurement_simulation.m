%%
n = 300; 


y0 = 36; % mean of prob distr deg
sigma_y = 4; % real width of prob distribution deg
%A = sqrt(2*pi)*sigma_y*56; % amplitude factor for prob distribution 
sigma_exp = 1.5; % experimental error

y = normrnd(y0, sigma_y, n, 1); % prob distribution values
dy = normrnd(0, sigma_exp, n,1); % experimental error
data = y + dy; % measured angles

sqrt(sigma_y^2+sigma_exp^2)

% FIT ONE GAUSS
options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',10000,'TolFun',1e-9,'MaxIter',10000, 'TolX', 1e-9); %'Algorithm','levenberg-marquardt',
%%
dx = 2;
xhist = 0:dx:80;
xplot = 0:0.1:80;
n = hist(data, xhist);
n0 = hist(y, xhist);

% fit measured data
p_init = [mean(data) std(data) max(n)];
[p,r,J,COVB,mse] = nlinfit(xhist, n, @gauss1d, p_init, options);
ci = nlparci(p,r,'Jacobian',J);
p(2) = abs(p(2));
p_err =  abs(p-ci(:,1)')  ;

% fit original data
p_init = [mean(y) std(y) max(n0)];
[p0,r,J,COVB,mse] = nlinfit(xhist, n0, @gauss1d, p_init, options);
ci = nlparci(p,r,'Jacobian',J);
p0(2) = abs(p0(2));
p0_err =  abs(p0-ci(:,1)')  ;

% plot results
close all
plot( xhist, n0, 'b.-', xplot, gauss1d(p0, xplot), 'b--', xhist, n, 'r.-', xplot, gauss1d(p, xplot), 'r'), hold on
legend({ 'Original', 'Original', 'Measured-fit','Original-fit'})



%% check if error of measurment is NORMAL distributed
a = [1	99	88	5	78	-86.332	12	78.160
2	110	84	47	65	-53.707	12	79.942
3	98	89	5	76	-86.236	12	75.854
4	109	85	47	60	-51.459	12	75.749
5	98	89	4	79	-87.101	12	78.752
6	107	82	48	65	-53.556	12	81.001
7	99	92	4	70	-85.914	12	70.467
8	108	83	43	57	-53.448	12	71.467
9	99	89	4	72	-86.028	12	72.130
10	108	84	46	61	-52.980	12	76.134
11	98	91	5	74	-86.135	12	74.192
12	109	86	41	52	-51.746	12	65.958
13	99	89	4	68	-85.795	12	67.805
14	107	83	50	66	-52.306	12	82.738
15	99	90	3	64	-87.316	12	64.070
16	108	84	47	62	-52.836	12	78.067
17	99	89	2	73	-88.431	12	73.027
18	108	84	48	66	-53.556	12	81.537
19	98	90	4	77	-87.026	12	76.789
20	106	84	56	73	-52.017	12	92.209
21	100	91	2	89	-89.349	12	88.682
22	108	83	46	62	-53.427	12	77.667
23	97	90	7	81	-85.061	12	81.302
24	109	87	56	69	-51.340	12	89.124
25	98	89	7	83	-85.865	12	83.267
26	109	84	46	63	-53.276	12	78.203
27	99	89	2	75	-89.236	12	75.019
28	109	84	48	67	-53.820	12	82.614
29	99	90	5	80	-86.424	12	79.823
30	108	84	43	58	-54.090	12	72.003
31	98	92	5	81	-86.468	12	80.802
32	105	82	47	62	-52.836	12	77.536
33	99	87	4	85	-87.306	12	85.094
34	107	84	47	66	-54.545	12	80.753
35	99	88	3	85	-87.979	12	84.720
36	109	84	50	67	-53.267	12	83.333
37	98	87	4	77	-87.769	12	76.754
38	110	85	45	63	-54.028	12	77.150
39	98	89	3	77	-88.493	12	76.713
40	108	84	45	63	-54.462	12	77.150
41	100	90	1	88	-89.341	12	87.677
42	109	85	45	61	-53.584	12	75.605
43	100	90	4	73	-87.647	12	73.092
44	108	84	43	61	-54.819	12	74.360
45	98	89	7	83	-85.179	12	82.991
46	108	83	56	74	-52.883	12	92.535
47	99	88	4	73	-87.647	12	73.092
48	110	85	41	56	-53.791	12	69.405
49	98	91	6	72	-84.447	12	71.946
50	109	85	48	65	-53.130	12	80.534
51	99	88	6	86	-86.009	12	85.877
52	109	84	43	57	-52.970	12	70.933
53	98	89	4	83	-87.207	12	82.780
54	107	82	48	65	-53.556	12	81.001
55	99	91	3	79	-87.852	12	79.390
56	108	83	51	67	-53.673	12	84.267
57	99	91	4	87	-87.368	12	87.092
58	108	83	44	62	-54.638	12	76.220
59	98	89	3	67	-87.436	12	67.067
60	110	86	46	61	-52.980	12	76.400
61	98	89	4	65	-86.479	12	65.123
62	110	84	43	63	-55.257	12	76.190
63	100	91	3	75	-87.709	12	74.714
64	109	85	45	61	-53.584	12	76.269
65	100	91	3	69	-88.340	12	69.385
66	109	83	42	58	-54.090	12	71.880
67	99	91	3	76	-87.739	12	75.740
68	110	86	47	62	-53.276	12	78.067
69	98	93	6	77	-85.601	12	77.566
70	108	83	45	60	-53.130	12	75.467
71	99	89	4	80	-87.138	12	80.417
72	110	85	46	64	-53.865	12	78.546
73	99	92	3	69	-87.510	12	68.748
74	109	83	48	64	-53.130	12	79.800
75	97	88	5	81	-86.468	12	81.487
76	109	84	44	59	-53.915	12	73.134
77	100	92	3	74	-87.678	12	74.394
78	110	84	45	62	-54.638	12	76.414
79	97	90	1	83	-90	12	83.336
80	108	83	44	59	-53.286	12	73.401
81	98	85	4	75	-86.135	12	74.792
82	107	84	42	54	-52.125	12	68.878
]; % FS-4sp1-30-ortho (particle 12)

b = abs(a(1:2:end,6)-a(2:2:end,6));
data = b;
%%
dx = 1;
xhist = 0:dx:80;
xplot = 0:0.1:80;
n = hist(b, xhist);

% fit measured data
p_init = [mean(data) std(data) max(n)/sum(n)];
[p,r,J,COVB,mse] = nlinfit(xhist, n/sum(n), @gauss1d, p_init, options);
ci = nlparci(p,r,'Jacobian',J);
p(2) = abs(p(2));
p_err =  abs(p-ci(:,1)')  ;

% plot results
close all
fig_dim =[10 7.5];
cur_fig = figure('Visible','on', 'PaperPositionMode', 'manual','PaperUnits','centimeters','PaperPosition', [0 0 fig_dim(1) fig_dim(2)], 'Position', [0 scrsz(4) fig_dim(1)*40 fig_dim(2)*40]);

 
bar(  xhist, n/sum(n)), hold on
%plot(  xhist, n/sum(n), 'b.-', xplot, gauss1d(p, xplot), 'r', xplot, gauss1d([mean(b) std(b) 1/sqrt(2*pi*std(b))], xplot), 'g'), hold on
plot(   xplot, gauss1d(p, xplot), 'r', xplot, dx*gauss1d([mean(b) std(b) 1/sqrt(2*pi*std(b))], xplot), 'g', 'linewidth', 2), hold on

legend({ 'Data', 'Gaussian fit', 'Normal distribution'})
ylabel('rel Frequency')
xlabel('Angle [deg]')

set(gca, 'XLim', [25 45])
set(gca, 'YLim', [0 1.5*max(n/sum(n))])

print(cur_fig, '-dtiff', '-r600' , [cd filesep 'Estimation_of_Measurement_error.tif'])

