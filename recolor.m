function recolor(sys)
% This function will set the screen color of a Simulink model specified by
% the input parameter sys and all the subsystems below sys.
% For example, this function can be combined with the command gcs (get
% current system), to recolor all the blocks below a specified level:
% recolor(gcs)
%
% The off-white color is a personal preference, easier on the eyes than the
% stark black and white default color scheme. KAB

% Set the screen color of the current system
set_param(sys,'ScreenColor','[1.000000, 0.964600, 0.724400]'); %off-white
% set_param(sys,'ScreenColor','white');

% Look for subsystems
blks=find_system(sys);

% Set the screen color of the subsystems
for i=2:length(blks)
    temp=get_param(blks(i),'BlockType');
    if strcmp(temp,'SubSystem')
        %off-white
        set_param(blks{i},'ScreenColor','[1.000000, 0.964600, 0.724400]');
        %        set_param(blks{i},'ScreenColor','white');   %white
    end     %if
end
