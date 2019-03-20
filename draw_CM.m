function Output_Image = draw_CM(Image, P, R)
% draw_CM - draw Center of Mass graphic at
% location P (2D - x & y coordinates in pixels)
% with circle radius R
%
%   NOTE! THIS VERSION WORKS WITH PIXEL DIMENSIONS

Output_Image    = Image;

PX  = P(1);
PY  = P(2);

for x = (PX-R):(PX+R),
    for y = (PY-R):(PY+R),
        
        color = 'n';    % set default pixel color: NO color
        
        %   Check which Quadrant we are in,
        %   upper left - white,
        %   upper right - black
        %   lower left - black
        %   lower right - white
        %   Since we will be plotting on a HOLD ON image, use value 64 for
        %   "white" - brightest colormap color, 0 for "black"
        rx      = x - PX;
        ry      = y - PY;
        
        tempR   = sqrt(rx^2+ry^2);
        
        if tempR < (R+1),
            if (rx < 0 && ry >= 0),         % upper left
                color = 'w';
            elseif (rx >= 0 && ry >= 0),    % upper right
                color = 'b';
            elseif (rx < 0 && ry < 0),      % lower left
                color = 'b';
            elseif (rx >= 0 && ry < 0),     % lower right
                color = 'w';
            end
        end
        
        %         %   fill in circle in OPPOSITE color
        %         if abs (tempR - R) < 1,
        %             if color == 'w',
        %                 color = 'b';
        %             elseif color == 'b',
        %                 color = 'w';
        %             end
        %         end
        
        %   fill in circle in BLACK color
        if abs (tempR - R) < 1,
            
            color = 'b';
        end
        
        % plot the actual pixel in the correct color
        
        if color == 'w',
            % image(x,y,64)
            Output_Image(y, x) = 64;
        elseif color == 'b',
            % image(x,y,0)
            Output_Image(y, x) = 0;
        end
        
    end % for y
end % for x

