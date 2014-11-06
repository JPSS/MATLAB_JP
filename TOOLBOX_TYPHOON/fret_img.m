function [ fret ] = fret_img( dd, da, gamma, area )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%

dd_sub = dd(area(2):area(2)+area(4) , area(1):area(1)+area(3) );

da_sub = da(area(2):area(2)+area(4) , area(1):area(1)+area(3) );

%% correct for bg
dd_sub = dd_sub - min(min(dd_sub));
da_sub = da_sub - min(min(da_sub));


%{
h_mean = fspecial('average', 5);
dd_sub = imfilter(dd_sub, h_mean);
da_sub = imfilter(da_sub, h_mean);
%}



fret = da_sub ./ (da_sub + gamma * dd_sub);

m_dd = 0.35*mean(mean(dd_sub));
m_da = 0.335*mean(mean(da_sub));

for i=1:size(dd_sub,1)
   for j =1:size(dd_sub,2)
      if da_sub(i,j) < m_da 
          fret(i,j) = 0;
      end
      if dd_sub(i,j) < m_dd 
          fret(i,j) = 0;
      end
   end
end



scrsz = get(0,'ScreenSize');
cur_fig = figure('Visible','on','OuterPosition',[ 1 scrsz(4) scrsz(3)*0.5 scrsz(4)/1.5], 'PaperPositionMode', 'auto'); % figure('Visible','off');%left bottom width height

colormap('Gray');
imagesc(fret, [0 1]), colorbar
title({'FRET efficiency', ['gamma = ' num2str(gamma)]})

end

