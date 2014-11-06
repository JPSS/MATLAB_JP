%%
img = imread([pname filesep fnames{1}]);
ref_full = img(4*area(1):4*area(2), 4*area(3):4*area(4));
imagesc(ref_full), colormap gray, colorbar, axis image



%%
dalpha = 1;
box_size = 4*50;

alpha = 0:dalpha:359;
n_rot = length(alpha); % number of rotations

box_size_template = ceil(2*box_size/2/cos(pi/4));%200;
box_size_template = int16(box_size_template + mod(box_size_template, 2) +1 ) ; % make it even

dx = (box_size_template-box_size-1)/2;



% generate library
lib = zeros(box_size+1, box_size+1, length(alpha));
for j=1:n_rot
    tmp = imrotate(ref_full, alpha(j), 'crop');
    lib(:,:,j) = tmp(dx:dx+box_size, dx:dx+box_size);
end

%%


close all

for i=1:n_img
    tmp = imread([pname filesep fnames{i}]);
    %img = imresize(tmp(1:2048,1:2048),[512 512], 'nearest'); %bin image 4x4 for faster image processing    
    img = tmp(1:2048,1:2048);    

    N = 4*512;
    w = box_size/2;
    p= zeros(2*w+1, 2*w+1, size(peaks{t,i},1));
    for j=1:size(peaks{t,i},1)
        x = 4*peaks{t,i}(j,1);
        y = 4*peaks{t,i}(j,2);
        
        p(max(1,w-y+2):min(2*w+1, 2*w+1-y-w+N),max(1,w-x+2):min(2*w+1, 2*w+1-x-w+N),j) = img(max(1, y-w):min(N,y+w) , max(1, x-w):min(N,x+w) );
        
        
        
        cc = zeros(n_rot,4);
        for r=1:n_rot % loop through rotations
            tmp = normxcorr2(lib(:,:,r), p(:,:,j)); % x-correlate
            %xcor_img(:,:,r) = tmp(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);
           % img2(:,:,i) = img2(:,:,i) + tmp.^2;
           
           tmp2 = tmp(w:box_size+w, w:box_size+w);
           [cc_tmp, cc_imax] = max(tmp2(:));
           cc(r,1) = r;
           cc(r,2) = cc_tmp;
           [cc(r,3), cc(r,4)] = ind2sub(size(tmp2),cc_imax);
           %if mod(r,10)==0
           %     imagesc(tmp(w:box_size+w, w:box_size+w), [-0.2 0.2]), colormap gray, colorbar, hold on
           %     plot(x, y, 'ro')
           %     plot(y, x, 'go')

             %   hold off
             %   pause
          % end

        end    
     
        ind = find((cc(:,3)<110 & cc(:,3)>90 & cc(:,4)<110 & cc(:,4)>90));
        
        [ccmax, imax] = max(cc(ind,2));
        ialpha = cc(ind(imax),1);
        
        subplot(1, 3, 1)
        if size(ind)>0
        imagesc(lib(:,:,ialpha)), colormap gray, colorbar, axis image
        end
      %  title(num2str(peaks{t,i}(j,4)))
        
        subplot(1,3,2)
        imagesc(p(:,:,j)), colormap gray, colorbar, axis image
        title('Partilce')
        subplot(1,3,3)
        plot(cc(:,1), 1000*cc(:,2), 'r', cc(:,1), cc(:,3), 'g', cc(:,1), cc(:,4), 'b')        
        pause
        close all
        
        
    end
    
    
    
    

end
 %%
 close all
 for i = 1:size(p,3)
     cc = zeros(n_rot,1);
        for r=1:n_rot % loop through rotations
            tmp = normxcorr2(lib(:,:,r), p(:,:,i)); % x-correlate
            %xcor_img(:,:,r) = tmp(box_size/2+1:end-box_size/2, box_size/2+1:end-box_size/2);
           % img2(:,:,i) = img2(:,:,i) + tmp.^2;
           %if mod(r,10)==0
           %imagesc(tmp(25:75, 25:75), [-0.2 0.2]), colormap gray, colorbar
           % pause
           %end
           tmp2 = tmp(25:75, 25:75);
           cc(r) = max(tmp2(:));
        end    
     
        [ccmax, imax] = max(cc);
        subplot(1, 2, 1)
        imagesc(lib(:,:,imax)), colormap gray, colorbar, axis image
        subplot(1,2,2)
        imagesc(p(:,:,i)), colormap gray, colorbar, axis image
        pause
 end
 
 %%
    %loop through images
    img2 = zeros(img_size+box_size, img_size+box_size, n_img); % stores sum of quadratic xcoef
   % img2 = zeros(img_size, img_size, n_img); % stores sum of quadratic xcoef
    img3 = zeros(img_size, img_size, n_img); % stores maximum of cor-coef of all rotations
    img3_index = zeros(img_size, img_size, n_img); % stores index of maximum

