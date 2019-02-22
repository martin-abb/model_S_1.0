            %----------------------------------------------------------------------
            f4=figure;
            image(img_proc_pick);
            title('MIT Princeton Pick Points')
            
            r1=rectangle('Position',[0,0,648,30]);
            set(r1,'FaceColor','w')
            
            if 0, % the first version gives text interpreter error
                data_dir_text   = mk_str(data_dir);
                f4t = text(10, 10, data_dir_text);
            else
                % f4t = text(10, 10, data_dir, 'interpreter','none');
                f4t = text(10, 10, last_dir, 'interpreter','none');
            end
            
            set(f4t,'FontSize',10);
            set(f4t,'FontWeight','Bold');
            set(f4t,'Color','k');
            
            %----------------------------------------------------------------------
  

            figure(f4);
            axis equal
            
            %***    title('MIT Princeton Pick Points with "MIT Score : ABB Score" (RED - disagreement)')
            title(mk_str([ totename '    MIT : ABB [Score, Freq, Amp] (RED - delta > 30%)' ]))
            for pick_point_no = 1:N_pick,
                
                MIT_Score       = round(100*All_Metrics(pick_point_no,1) );
                ABB_Score       = round(100*All_Metrics(pick_point_no,3) );
                
                Pick_no_string  = num2str( All_Metrics_Percent( pick_point_no, 1 ) );
                MIT_string      = num2str( All_Metrics_Percent( pick_point_no, 2 ) );
                ABB_string    = num2str( All_Metrics_Percent( pick_point_no, 3:5 ) );
                % Note! ABB has Score_metric            = [ Score Score_freq Score_amp ];
                
                
                Score_threshold     = 30;
                if abs(MIT_Score - ABB_Score)>= Score_threshold,
                    score_color = 'r';
                else
                    score_color = 'b';
                end
                
                text_x_offset   = 20;
                text_y_offset   = 0;
                
                t92=text(pixel_file_contents(pick_point_no,1) + text_x_offset, ...
                    pixel_file_contents(pick_point_no,2) - text_y_offset, ...
                    [  MIT_string ' : ' ABB_string ] );
                
                set(t92,'FontSize',12);
                set(t92,'FontWeight','Bold');
                set(t92,'Color',score_color);
            end
            
