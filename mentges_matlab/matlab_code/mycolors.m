%% Prints color RGB values or returns RGB of specified color
% - mycolors(), prints all RGB to command window and opens overview plot.
% - mycolors('blue') returns RGB of blue in vector.
% - use uicolorpicker() to pick new colors.

function [varargout] = mycolors(varargin)
    
colornames = {'darkblue', 'lightdarkblue', 'greyblue', 'paleblue',...
    'verypaleblue', ...
    'blue', 'lightblue', 'darkorange', 'orange', 'lightorange', 'paleorange', ...
    'orangeyellow', 'yellow', 'lightyellow', 'darkpurple', 'purple','lightpurple', ...
    'darkred', 'mediumred', 'red', 'lightred', 'redorange', ...
    'bluegreen', 'darkgreen', 'green', 'green2', 'green3', 'lightgreen', ...
    'palegreen', 'verypalegreen',...
    'pink', 'palepink'};
    
RGB = [0.03, 0.18, 0.41; ... % darkblue
    0.036, 0.26, 0.6; ... % lightdarkblue
    0.31, 0.52, 0.68; ... % greyblue [80, 133, 175]/255
    0.55, 0.75, 0.89; ... % paleblue
    [188 217 242]/255; ... % verypaleblue
    0, 0.44, 0.74; ... % blue
    0.30, 0.75, 0.93; ... % lightblue
    0.55, 0.1, 0.05; ... % darkorange
    0.85, 0.4, 0.15; ... % orange
    0.94, 0.49, 0.2; ... % lightorange [239, 125, 50]/255
    [249, 181, 122]/255; ... % paleorange
    0.92, 0.64, 0.12; ... % orangeyellow
    0.9, 0.7, 0.1; ... % yellow
    1,   0.97, 0.5;... % light yellow
    0.4, 0.17, 0.43;... % darkpurple
    0.55, 0.35, 0.65;... % purple
    0.8, 0.6, 0.9;... % lightpurple
    0.504 0.064 0.144;... % darkred
    0.7, 0.057, 0.09;... % mediumred
    0.95, 0.05, 0.05;... % red
    0.9, 0.6, 0.5;... % lightred
    255/255 80/255 0;... % redorange
    0.25, 0.85, 0.75;... % bluegreen
    0.135 0.45 0.36;... % darkgreen
    0.4, 0.7, 0.2; ...  % green
    [0.63, 0.80, 0.49]*0.9; ... % green2
    0.63, 0.80, 0.49; ... % green3
    0.7, 0.85, 0.2;... % lightgreen
    0.56 ,0.80, 0.5764;... % palegreen
    0.73 ,0.87, 0.74;... % verypalegreen
    228/255, 59/255, 127/255;... % pink
    [240, 199, 235]/255]; % palepink

if ~isempty(varargin) && any(ismember(varargin{:}, colornames))
    ind = ismember(colornames, varargin{:});
    varargout{1} = RGB(ind,:);
elseif isempty(varargin)
    fprintf('Color\t\t');
    fprintf('R    G    B\t\n\n');
    for k = 1:length(colornames)
        fprintf('%s     \t', colornames{k})
        fprintf('%1.2f ', RGB(k,1)) 
        fprintf('%1.2f ', RGB(k,2)) 
        fprintf('%1.2f\n', RGB(k,3))  
    end         
    figure('color', 'white', 'position', [-563,166,268,701])
    imagesc([1:length(RGB)]')
    colormap(RGB)
    text(ones(1, length(RGB))*0.6, 1:length(RGB), colornames)
elseif strcmp(varargin{1}, 'names')
    varargout{1} = colornames;
else
    fprintf('\n\nColor is not valid. Choose from:\n');
    fprintf('  %s\n', colornames{:});
    fprintf('\nOr call mycolors() to see all RGB.\n');
end

end