function [ gamma ] = estimate_gamma(dd, da, aa)
%ESTIMATE_GAMMA Summary of this function goes here
%   Detailed explanation goes here
display('Select donor-only trace!')
donly = cell(1, 7);
[l1, y , pos] = integrate_lane(dd, 'd-only (D -> D)');
donly{1, 1} =  +l1;
donly{1, 2} =     +transpose( sum(transpose(  da(y  ,   pos(1):pos(1)+pos(3)   )  )) );
donly{1, 3} =    +transpose( sum(transpose(  aa(y  ,   pos(1):pos(1)+pos(3)   )  )) );
donly{1, 4} = double(y);
donly{1,5} = pos;
lane_name = inputdlg({'Name of lane:', 'Length of spacer:'}, 'Lane properties' , 1, {'d-only', '10'} );
donly{1,6} = lane_name{1};
donly{1,7} = str2double(lane_name{2});

display('Select corresponding fret-trace!')
lanes = cell(1, 7);
[l1, y , pos] = integrate_lane(dd, 'D+A (D -> D)');
lanes{1, 1} =  +l1;
lanes{1, 2} =     +transpose( sum(transpose(  da(y  ,   pos(1):pos(1)+pos(3)   )  )) );
lanes{1, 3} =    +transpose( sum(transpose(  aa(y  ,   pos(1):pos(1)+pos(3)   )  )) );
lanes{1, 4} = double(y);
lanes{1,5} = pos;
lane_name = inputdlg({'Name of lane:', 'Length of spacer:'}, 'Lane properties' , 1, {'10bp', '10'} );
lanes{1,6} = lane_name{1};
lanes{1,7} = str2double(lane_name{2});


%% fit
p = fit_lane(donly{1,4}, donly{1,1}, 15, 1, 1, 0); % find the limits
p_a = fit_lane(lanes{1,4}, lanes{1,2}, 15, 1, 1, 0); % find the limits
%%

sigma = 1;
y_min = p(1) - sigma * p(2);
y_max = p(1) + sigma * p(2);
   
% integrate each lane
y = lanes{1,4};
i_0 =  floor(length(y)*(y_min-y(1))/(y(end)-y(1)+1) ) +1;
i_end =  floor(length(y)*(y_max-y(1))/(y(end)-y(1)+1) ) +1;

I_dd = sum(lanes{1,1}(i_0:i_end));
I_fret = sum(lanes{1,2}(i_0:i_end));


y_min_a = p_a(1) - sigma * p(2); % same sigma different centreer
y_max_a = p_a(1) + sigma * p(2);

y = donly{1,4};
i_0 =  floor(length(y)*(y_min_a-y(1))/(y(end)-y(1)+1) ) +1;
i_end =  floor(length(y)*(y_max_a-y(1))/(y(end)-y(1)+1) ) +1;


I_donly = sum(donly{1,1}(i_0:i_end));
I_donly_fret = sum(donly{1,2}(i_0:i_end));


gamma = (I_fret-I_donly_fret) / (I_donly-I_dd)

hh = plot( donly{1,4}, donly{1,1} , 'g', donly{1,4}, donly{1,2} , 'b', lanes{1,4}, lanes{1,1}, 'g--', lanes{1,4}, lanes{1,2}, 'b--' ), hold on
%legend( 'D-only D->D', 'D-only D->A','D->D','D->A')
h =vline(y_min, 'k')
set(h, 'LineWidth',1);
h =vline(y_max, 'k')
set(h, 'LineWidth',1);

g =vline(y_min_a, 'k--');
set(g, 'LineWidth',1)

g =vline(y_max_a, 'k--');
set(g, 'LineWidth',1)
legend([hh; h; g] , {'D-only D->D', 'D-only D->A','D->D','D->A', 'D-only', 'D+A'})

hold off
%%



end

