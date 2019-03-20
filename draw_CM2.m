function success = draw_CM2(P, R, dx)
%
%   DOES NOT WORK CORRECTLY YET!
%
%
% draw_CM2 - draw Center of Mass graphic at
% location P (2D - x & y coordinates in pixels)
% with circle radius R
%
%   NOTE! THIS VERSION CM2 WORKS WITH LINEAR DIMENSIONS

PX  = P(1);
PY  = P(2);

for x = (PX-R):dx:(PX+R),
    for y = (PY-R):dx:(PY+R),
        
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
        if abs (tempR - R) < dx,
            
            color = 'b';
        end
        
        % plot the actual pixel in the correct color
        
        if color == 'w',
            image(x,y,64)
        elseif color == 'b',
            image(x,y,0)
        end
        
    end % for y
end % for x


success         = 1;

end

