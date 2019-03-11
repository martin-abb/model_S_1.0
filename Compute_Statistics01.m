%
%   Compute_Statistics01.m
%
%   Compute statistics for ABB Suction Model as compared with MIT model
%   from real experimental data
%
%   2019-03-01      v01 Initial version

%   Use NA_number as the Not Applicable number for proper column spacing in
%   tables

NA_number       = -9999;

all_MIT_scores_failed  = [];
all_ABB_scores_failed  = [];


%   Set up looping over directories

top_dir         = 'C:\Users\USMAKRU\Documents\ABB Local\Logistics_Project\operation_data';
dirmask_top     = 'empty_mtote_*';
dir_content_top = dir( [ top_dir '\' dirmask_top ] );
num_dir_top     = length(dir_content_top);

base_dir_local  = 'C:\Users\USMAKRU\OneDrive - ABB\2018\Logistics\Suction Cup Modeling\Suction_1.0_ML\data';


%for dir_no = num_dir_top:(-1):1,  % traverse directories from newest to oldest

%disp('*** SPECIAL PROCESSING *** Starting at directory 222 and going downwards...')
for dir_no = 222:(-1):4,  % traverse directories from newest to oldest

%for dir_no = 1:num_dir_top,
%***for dir_no = 201:220, %*** use 3 directories for testing, full number is 247
    
    if 0,
        %----------------------------------------------------------------------
        %   For comparison script
        data_dir_local  = 'C:\Users\USMAKRU\Documents\ABB Local\Logistics_Project\operation_data\empty_mtote_2019-01-28-16-09';
        %data_dir_local  = 'C:\Users\USMAKRU\Documents\ABB Local\Logistics_Project\operation_data\empty_mtote_2019-02-13-09-59';
        
        
    else
        data_dir_local  = [ dir_content_top(dir_no).folder '\' dir_content_top(dir_no).name ];
    end
    
    
    %----------------------------------------------------------------------
    %   Load sensor data images for testing
    %
    
    dirmask     = '*.tote.txt';
    data_dir    = data_dir_local;
    data_dir_parts  = regexp(data_dir, '\','split');
    last_dir    = data_dir_parts{end};
    
    dir_content     = dir( [ data_dir '\' dirmask ] );
    num_images      = length(dir_content);
    
    %----------------------------------------------------------------------
    %   Prepare a directory where to save images & log file
    %----------------------------------------------------------------------
    
    base_dir        = pwd;
    
    log_dir         = [ data_dir '\matlab_results' ];
    
    if exist(log_dir)==7,   %   Make sure matlab_results directory exists
        
        cd(log_dir);
        
        summary_name      = [ 'Summary_' last_dir '.txt' ];
        
        %  Import summary table with the following entries:
        %   pick number (i.e. image number)
        %   MIT Score
        %   ABB Scores
        %   Success flag
        %   Error code
        %   Chosen Pick Point #
        
        Summary_table = import_summary(summary_name);
        tot_num_cols    = 7;    % total number of columns in metrics tables
        
        tot_num_summary_cols    = tot_num_cols + 1; % total number of columns in summary tables, ADD the pick point #
        
        %Summary_table   = zeros(num_images,tot_num_summary_cols);
        
        disp([ 'Processing summary file in directory:  ' log_dir ])
        disp(' ')
        
        
        
        %   Process content of Summary File
        
        col_Success         = 6;
        col_Error_Code      = 7;
        col_PickPoint_no    = 8;
        col_MIT             = 2;
        col_ABB             = 3;
        
        %*** GIVES ONE TOO HIGH NUMBER num_of_failures     =  sum( Summary_table(:,col_Success)==1 | Summary_table(:,col_Success)==0 ) - sum(Summary_table(:,col_Success)==1);
        num_of_failures     =  sum( Summary_table(:,col_Error_Code)~= NA_number ) - sum( Summary_table(:,col_Error_Code)== 0 );
        
        disp([ 'Number of failures: ' num2str(num_of_failures) ])
        disp(' ')
        disp('      Pick #         MIT         ABB        Freq         Amp     Success       Error  PickPoint#')
        disp(Summary_table)
        diary off
        
        %   Extract key metrics & statistics
        
        Error_Codes     = Summary_table(:, col_Error_Code);
        
        ind_Failures    = find( (Error_Codes == -11) | (Error_Codes == -20) );  % -11 failed grasp, -20 item lost in transport
        
        MIT_scores      = Summary_table( ind_Failures, col_MIT);
        ABB_scores      = Summary_table( ind_Failures, col_ABB);
        
        %   prepare for processing of next folder
        cd(base_dir)
        
        close all
        
    end %     if exist(log_dir)==2,   %   Make sure matlab_results directory exists
    
    %   Store results
    all_MIT_scores_failed      = [ all_MIT_scores_failed ; MIT_scores ];
    all_ABB_scores_failed      = [ all_ABB_scores_failed ; ABB_scores ];
    
    
    
    
end % for dir_no = 1:num_dir_top,

%   Plot MIT & ABB Score overview for failed picks

f1  = figure;

set(f1, 'DefaultLineLineWidth', 3);

plot(all_MIT_scores_failed, 'bo');
hold on
plot(all_ABB_scores_failed, 'ro');
legend('MIT','ABB')
ylim([0 100])
ylabel('Pick Score [%]')
xlabel('Failed pick #')
title('Pick Score for failed picks')
grid on

%----------------------------------------------------------------------
%   Return to home directory
%----------------------------------------------------------------------

cd(base_dir)

