classdef movie < handle
    %TRACE_SELECTION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frames; %all images to be read
        N_read = 300; % number of images to read in one portion
        counter = 1; % internal counter for reading the movie
        
        sequence; % sequence to be read, e.g. 101010
        first; % first image to be read
        last; % last image to be read
        
        pname; %pathname of file location
        fname; %filename of file
        
        sizeX; % number of pixel in X-dimension
        sizeY; % number of pixel in Y-dimension
        mov_length; % number of frames in thw whole movue
        
        info; %fits info
        h_min; %minimal heigth for peak fidning
        
        input; % 0=fits, 1=tiff-stack
        fnames; % cell with all filenames, only for tiff-stack
        
        N_frame_per_fits; % stores the number of frames in one fits-file
    end
    
    methods
        %constructor
        function obj = movie(pname, fname, first, last, sequence) % fname is the filename of the first fits-file
            obj.sequence = sequence;
            obj.pname = pname;
            obj.fname = cell(1,1);
            obj.fname{1} = fname;
            
            if strcmp(fname(end-2:end), 'tif')
                obj.input = 1; % read tiff data
            else
                obj.input = 0; % read fits data
            end
            
            if obj.input == 1 % tiff-stack
                tmp = dir([pname filesep '*.tif']);
                obj.fnames = {tmp.name};
                obj.info = 'Tiff-Stack info';
                obj.sizeX = size(imread([pname filesep obj.fnames{1}]),2);
                obj.sizeY = size(imread([pname filesep obj.fnames{1}]),1);
                obj.mov_length = length(obj.fnames);
            else % fits
                obj.info = cell(1,1);
                obj.info{1} = fitsinfo([obj.pname filesep obj.fname{1}]);
                
                obj.N_frame_per_fits = obj.info{1}.PrimaryData.Size(3);
                
                obj.sizeX = obj.info{1}.PrimaryData.Size(1); 
                obj.sizeY = obj.info{1}.PrimaryData.Size(2);
                
                tmp = dir([pname filesep fname(1:end-4) '_X*.fits']); %returns additional  change * to wildcard for 1-2 character/integers
                for i=1:length(tmp)
                    obj.fname = [obj.fname tmp(i).name];
                    obj.info = [obj.info fitsinfo([obj.pname filesep tmp(i).name])];
                end
                
                f_tot = 0; % calculate total number of frames
                for i=1:length(obj.fname)
                    f_tot = f_tot + obj.info{i}.PrimaryData.Size(3);           
                end
                obj.mov_length = f_tot;
            end
            
            
            obj.first = first;
            if last == -1
                obj.last = obj.mov_length;
            else
                obj.last = last;
            end
            obj.frames = obj.getFrames(obj.sequence, obj.first, obj.last);
            
        end
        
        %generate a list of images to be read from the movie
         function frames = getFrames(obj, sequence, first, last)
            frames = [];
            tmp = first:last;
            for i=1:length(tmp)
               if(  sequence(  mod(i-1,size(sequence,2))+1 )  )
                  frames = [frames tmp(i) ];       
               end    
            end
         end
         
         %reads one frame
         function [img] = readFrame(obj,framenumber)
             
              if obj.input == 1 % tif
                    img = double(imread([obj.pname filesep obj.fnames{framenumber}]));
              else
                   i_fits = ceil(framenumber/obj.N_frame_per_fits);    % index of fits file
                   framenumber_effektive = mod(framenumber-1, obj.N_frame_per_fits) +1;
                   img = fitsread([obj.pname filesep obj.fname{i_fits}],  'Info', obj.info{i_fits}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [framenumber_effektive framenumber_effektive] });  % fits
              end
         end
         
         % initialize the counter
         function  initRead(obj)
              obj.counter = 1;
         end
          
         
         %reads the next N_read frames
         function [tmp, frames_out, go_on] = readNext(obj)
            read_stop = obj.counter+obj.N_read-1; % last frame to read
            go_on = 1;
            if read_stop >= length(obj.frames)
                read_stop = length(obj.frames);
                go_on = 0;
            end
            
            
            tmp = zeros(obj.sizeX, obj.sizeY, read_stop-obj.counter+1);
            frames_out = obj.frames(obj.counter:read_stop);
            
            display(['Reading frame ' num2str(frames_out(1)) ' to ' num2str(frames_out(end))])
            for i=1:length(frames_out)
                tmp(:,:,i) = obj.readFrame(frames_out(i));
            end
            
            obj.counter = obj.counter + length(frames_out); 
         end
         
         
         
         
         % trace the movie 
         function [ traces, itraces, avg_frame ] = trace_movie(obj, h_min, r_find, r_integrate, min_length )
            traces = cell(0,1);
            itraces = cell(0,1);
            avg_frame = zeros(obj.sizeX, obj.sizeY);
            
            go_on = 1;
            obj.initRead;
            N = 0;
                     
            %display('Tracing movie... please wait')
            while go_on
                [movie, frames, go_on]  = obj.readNext;
                [traces, itraces] = append_traces(movie, traces, itraces, frames, h_min, r_find, r_integrate, min_length);
            
                avg_frame = avg_frame + sum(movie,3);
                N = N + length(frames);
            end
            avg_frame = avg_frame ./ N; 
            %display('Done tracing movie.')
         end
         
         % integrate specific reagions in movie
         function [itraces ] = traces_movie_position(obj, positions, r_integrate )
            itraces = cell(0,1);
            go_on = 1;
            obj.initRead;
            while go_on
                [movie, frames, go_on]  = obj.readNext;
                itraces = append_traces_to_position(movie, itraces, frames, positions, r_integrate);
            end
         end
         
         % determine peak-finding threshholds
         function [h_min, p_out] = get_h_min(obj, r_find, N_img)
            
             if exist('N_img', 'var') % use average frame
                 img = obj.average_image(N_img); % average frist N_img images
             else % use first frame if N_img is not specified
                if obj.input == 1 % tiff
                    img = double(imread([obj.pname filesep obj.fnames{obj.frames(1)}]));
                else % fits
                    img = fitsread([obj.pname filesep obj.fname{1}],  'Info', obj.info{1}, 'PixelRegion',{[1 obj.sizeX], [1 obj.sizeY], [obj.frames(1) obj.frames(1)] }); % read first frame                
                end
             end
            p = find_peaks2d(img, r_find, 0, 0); % finding all possible peaks p has x, y, height, height-bg, I, I-I_bg
            
            
            close all
            figure('units','normalized','outerposition',[0 0 1 1])
            img_mean = mean(img(:));
            img_std = std(img(:));            
            
            p_7std = p(find(p(:,4)>=7*img_std), :); % this is just an estimate, #of peaks found may vary since peak_find algorithm return height as int not double
            p_5std = p(find(p(:,4)>=5*img_std), :);
            p_3std = p(find(p(:,4)>=3*img_std), :);

            % plot
            subplot(1, 2, 1)
            imagesc(img), colorbar, axis image, colormap gray, hold on
            if size(p_3std,1)>0
                h(1) = plot(p_3std(:,1)+1, p_3std(:,2)+1, 'ro');
            end
            if size(p_5std,1)>0
                h(2) = plot(p_5std(:,1)+1, p_5std(:,2)+1, 'go');
            end
            if size(p_7std,1)>0
                h(3) = plot(p_7std(:,1)+1, p_7std(:,2)+1, 'bo');
            end
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
                  
            
            subplot(1, 2, 2)
            xhist = min(p(:,4)):5:max(p(:,4));
            n = hist(p(:,4), xhist);
            semilogy(xhist, sum(n)-cumsum(n)), hold on
            h(1) = vline(3*img_std, 'r');
            h(2) = vline(5*img_std, 'g');
            h(3) = vline(7*img_std, 'b');
            legend(h, {['3\sigma = ' num2str(round(3*img_std)) ], ['5\sigma = ' num2str(round(5*img_std)) ], ['7\sigma = ' num2str(round(7*img_std)) ]})
            set(gca, 'XLim', [0 xhist(end)])
            xlabel('Minimal height'), ylabel('# of peaks found')
            axis square
            
            
            % promp
            options.WindowStyle='normal';
            prompt={'Enter min heigth (default=5*sigma):'};
            def={num2str(round(5*img_std))};
            threshold = inputdlg(prompt, strcat('Enter threshold:'), 1, def, options);
            h_min = str2double(threshold(1));
            close all
            
            obj.h_min = h_min;
                         
            p_out = p(find(p(:,4)>=h_min),:);
         end
         
         
                
         
         
         % generate average image
         function [ avg_frame ] = average_image(obj, N_max )

            if N_max <= 0 % adjust to full movie lenght
                N_max = obj.mov_length; 
            end
            avg_frame = zeros(obj.sizeX, obj.sizeY);
            
            go_on = 1;
            obj.initRead;
            N = 0;
            while go_on
                [movie, frames, go_on]  = obj.readNext;
            
                if N+length(frames) < N_max
                    avg_frame = avg_frame + sum(movie,3);
                    N = N + length(frames);
                else
                    k = min(N_max-N, length(frames));
                    avg_frame = avg_frame + sum(movie(:,:,1:k),3);
                    N = N + k;
                    go_on = 0;
                end
            end
            avg_frame = avg_frame ./ N; 
         end
         
         
         % Fits a PSF to one spot in each frame in fit_frames, at given
        % positions fit_pos, using fit parameter sigma
        
        function [ pos_trace ] = fit_psf_to_movie(obj, fit_frames, fit_pos, sigma)
        
            w_fit = 3*sigma;
            s_x = sigma;
            s_y = s_x;
            pos_trace = [];
            counter = 0;
            i=1;
            
            while counter<20 && i<=size(fit_frames,1)
                cur_frame = obj.readFrame(fit_frames(i));     
                x_0 = fit_pos(i,1)+1;
                y_0 = fit_pos(i,2)+1;
                A = cur_frame(y_0,x_0);
                [X,Y] = meshgrid(max(x_0-w_fit,1):min(x_0+w_fit,512),max(y_0-w_fit,1):min(y_0+w_fit,512));
                Z = cur_frame(max(y_0-w_fit,1):min(y_0+w_fit,512),max(x_0-w_fit,1):min(x_0+w_fit,512));
                xyz_data=zeros(size(X,1)*size(X,2),3);
                xyz_data(:,1) = reshape(X, size(X,1)*size(X,2), 1); %make a vector out of matrix X
                xyz_data(:,2) = reshape(Y, size(Y,1)*size(Y,2), 1); %make a vector out of matrix X
                xyz_data(:,3) = reshape(Z, size(Z,1)*size(Z,2), 1); %make a vector out of matrix X
                bg = mean(xyz_data(:,3));

                param_init = [x_0 y_0 s_x s_y A-bg bg];
                options = optimset('Algorithm','levenberg-marquardt','display','off', 'MaxFunEvals',50000,...
                    'TolFun',1e-9,'MaxIter',1000, 'TolX', 1e-9); 
                [param, chi2, residual, exitflag, output] = lsqcurvefit(@gauss2d_bg, param_init,xyz_data(:,1:2), xyz_data(:,3),[],[], options);
                display([num2str(output.iterations) ', ' num2str(output.funcCount)])
                if exitflag <= 0
                        display(['WARNING: Fitting gaussian failed. Exitflag: ' num2str(exitflag)])
                        counter = counter + 1
                        if counter == 20
                            display('Moving on to next spot')
                        end
                end
 
                if exitflag > 0
                    counter = 0;
                end
                
            pos_trace = [pos_trace ; [param chi2]];
            i=i+1;
            end
        end


    end
    

    
    
    
end
        