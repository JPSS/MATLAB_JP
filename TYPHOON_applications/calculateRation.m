function [ R ] = calculateRation( A, B_in, status )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


    
    
    %[cc, shift, lane_da_shift] = xcorr2_bounded(lane_dd, lane_da, 5); % find best overlay of images
    [cc, shift, B] = xcorr2_bounded(A, B_in, 5, 0); % find best overlay of images
    
    %{
    %first algorith
    profileA = sum(A,2);
    profileB = sum(B,2);
    R(1,1) = max(profileA) ./ max(profileB);
    
    % second algorithm
    R(2,1) = sum(sum(A)) ./ sum(sum(B));
    
    %}
    
    %third algorithm
    p = polyfit(B(:), A(:), 1);
    p_raw = polyfit(B_in(:), A(:), 1);
    
    R = p(1);
    
    
    x = [min(B(:)) max(B(:))];
    
    
   % subplot(2,1,1)
   % plot(1:size(A,1), profileA, 'r', 1:size(B,1), profileB, 'g')
   % legend({'A', 'B'})
     
   if status
        %close all
        subplot(2,1,1)
        plot(B_in(:), A(:), 'b.', x, p_raw(1).*x+p_raw(2), 'r', x,  (sum(sum(A)) ./ sum(sum(B))).*x, 'g')
        legend({[num2str(p_raw(1))], num2str(sum(sum(A)) ./ sum(sum(B))) })
        
        subplot(2,1,2)
        plot(B(:), A(:), 'b.', x, p(1).*x+p(2), 'r', x,  (sum(sum(A)) ./ sum(sum(B))).*x, 'g')
        legend({[num2str(p(1))], num2str(sum(sum(A)) ./ sum(sum(B)))})
   end
    
    
    

end

