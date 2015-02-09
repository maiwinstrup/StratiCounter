function matchmaker_evaluate(varargin)
%
% MATCHMAKER_EVALUATE version 2.00, 16 March 2006, 
% Copyright (C) Sune Olander Rasmussen
% MATCHMAKER_EVALUATE is designed to be called from MATCHMAKER.

switch nargin
    case {0, 1, 2}% Only for development purposes
        disp('MATCHMAKER_EVALUATE needs input arguments, and is designed to be called from MATCHMAKER.');
    case {3, 4, 5, 6, 7}
%        try
eval([char(varargin{1}) '(varargin{[2 4:nargin]})']);
%        catch
%            disp(lasterr)
%            dummy = input('ERROR : Log the error and press enter to continue');
%            disp(' ');
%        end;
    otherwise
        disp('Wrong number of arguments when calling matchmaker_evaluate.m');
end;

%---

function evalopen(matchmakerfighandle, mp, core, masterno, current_mp);
for i = 1:length(mp)
    Nmpi134(i) = length(mp{i}(find(mp{i}(:,2)==1 | mp{i}(:,2)==3 | mp{i}(:,2)==4),1));
end;
if min(Nmpi134)<2
    errordlg('There must be at least two first order fix points for each core for the evaluate window to work. Evaluate window will not open.', 'Warning calling MATCHMAKER_EVALUATE', 'modal');
    return
end;

scrsize = get(0, 'screensize');
hgt = scrsize(4);
handles.fig = figure('position', [1 60 scrsize(3) scrsize(4)-100], 'name', 'Matchmaker Evaluation Tool', 'CloseRequestFcn', 'matchmaker_evaluate(''exit_Callback'',gcbo,[],guidata(gcbo))', 'nextplot', 'add', 'color', 0.9*[1 1 1], 'pointer', 'crosshair', 'toolbar', 'figure', 'Numbertitle', 'off', 'KeyPressFcn', 'matchmaker_evaluate(''keypressed_Callback'',gcbo,[],guidata(gcbo))', 'integerhandle', 'off');
handles.matchmakerfighandle = matchmakerfighandle; 
handles.mp = mp;
handles.core = core;
handles.masterno = masterno;
handles.N = length(mp);
handles.current_mp = current_mp([1 end]);;

dummyax = axes('position', [0 0 1 1], 'xlim', [0 1], 'ylim', [0 1], 'visible', 'off', 'nextplot', 'add', 'hittest', 'off');

font1 = 9;
font2 = 14;
eh = 0.029*710/hgt;

handles.title = text(0.005, 0.012, 'MATCHMAKER Evaluation Tool', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Fontsize', font2, 'fontweight', 'bold', 'parent', dummyax);
handles.text1 = text(0.52, 0.04, 'Comparison match point interval', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Fontsize', font1, 'parent', dummyax);
handles.text1 = text(0.69, 0.04, 'Match point ID', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Fontsize', font1, 'parent', dummyax);
handles.return = uicontrol('units', 'normalized', 'position', [0.81 0.012 0.08 eh], 'string', 'Return', 'style', 'pushbutton', 'callback', ['matchmaker_evaluate(''return_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
handles.exit = uicontrol('units', 'normalized', 'position', [0.91 0.012 0.08 eh], 'string', 'Exit', 'style', 'pushbutton', 'callback', ['matchmaker_evaluate(''exit_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
handles.lowmp = uicontrol('units', 'normalized', 'position', [0.41 0.012 0.08 eh], 'string', num2str(current_mp(1)), 'style', 'edit', 'callback', ['matchmaker_evaluate(''xlim_Callback'',gcbo,[],guidata(gcbo), -1)'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
handles.highmp = uicontrol('units', 'normalized', 'position', [0.55 0.012 0.08 eh], 'string', num2str(current_mp(end)), 'style', 'edit', 'callback', ['matchmaker_evaluate(''xlim_Callback'',gcbo,[],guidata(gcbo), 1)'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
handles.inc = uicontrol('units', 'normalized', 'position', [0.50 0.012 0.04 eh], 'string', '1', 'style', 'edit', 'callback', ['matchmaker_evaluate(''inc_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
handles.info = uicontrol('units', 'normalized', 'position', [0.65 0.012 0.08 eh], 'string', [], 'style', 'edit', 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
handles.diff = uicontrol('units', 'normalized', 'position', [0.31 0.012 0.08 eh], 'string', 'Depth diff', 'style', 'togglebutton', 'callback', ['matchmaker_evaluate(''diff_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center', 'value', 1);

handles.ax(1) = axes('position', [0.05 0.1 0.44 0.84], 'nextplot', 'add', 'box', 'on', 'fontsize', font1);
handles.axtitle(1) = title(['Depth vs. ' handles.core{handles.masterno} ' depth'], 'Fontsize', font1);
handles.ax(2) = axes('position', [0.55 0.1 0.44 0.84], 'nextplot', 'add', 'box', 'on', 'fontsize', font1);
handles.axtitle(2) = title(['(slope / ' handles.core{handles.masterno} ' slope) vs. ' handles.core{handles.masterno} 'depth'], 'Fontsize', font1);
xlim_Callback(handles.fig, handles, 1);
diff_Callback(handles.fig, handles)

set(handles.fig, 'handlevisibility', 'callback');
guidata(handles.fig, handles);

matchmaker('evaluate_Callback',handles.fig,[],handles, 'opening_evaluate'); % THIS LINE SHOULD BE UNCOMMENTED !!!!

%---

function evalreuse(evaluatefigurehandle, mp, core, masterno, current_mp);
figure(evaluatefigurehandle);
handles = guidata(evaluatefigurehandle);
handles.mp = mp;
handles.core = core;
handles.masterno = masterno;
handles.N = length(mp);
handles.current_mp = current_mp([1 end]);;

xlim_Callback(handles.fig, handles, 1);
guidata(handles.fig, handles);

%---

function keypressed_Callback(hObject, handles) % Translate keypress to appropriate button actions.
key = double(get(handles.fig, 'currentcharacter'));
if ~isempty(key)
    switch key
    case 28   %<-
        set(handles.lowmp, 'string', num2str(str2num(get(handles.lowmp, 'string'))-str2num(get(handles.inc, 'string'))));
        set(handles.highmp, 'string', num2str(str2num(get(handles.highmp, 'string'))-str2num(get(handles.inc, 'string'))));
        xlim_Callback(handles.fig, handles, 1);
        guidata(handles.fig, handles);
    case 29  %->
        set(handles.lowmp, 'string', num2str(str2num(get(handles.lowmp, 'string'))+str2num(get(handles.inc, 'string'))));
        set(handles.highmp, 'string', num2str(str2num(get(handles.highmp, 'string'))+str2num(get(handles.inc, 'string'))));
        xlim_Callback(handles.fig, handles, 1);
        guidata(handles.fig, handles);
    case {99, 67}   %c, C
        cursor_Callback(hObject, handles);
    case {114, 82}  %r, R
        figure(handles.matchmakerfighandle)
    case {120, 88}  %x, X
        exit_Callback(hObject, handles);
    otherwise % If key not defined, show info window
        disp(['Just for information : MATCHMAKER_EVALUATE undefined key callback : ' num2str(key)]);
        h_help = helpdlg({...
            'Available keyboard commands:';
            '<- = move selected core one frame back and accordianize'
            '-> = move selected core one frame forward and accordianize'
            'C  = change Cursor type'
            'R  = Return to Matchmaker main screen'
            'X  = eXit'});
    end;
end;

%---

function return_Callback(hObject, handles) % Change cursor type from crosshair to fullcrosshair and back
figure(handles.matchmakerfighandle)

%---

function diff_Callback(hObject, handles) % Change cursor type from crosshair to fullcrosshair and back
plotcurves(handles);
if get(handles.diff, 'value') == 1
    set(handles.axtitle(1), 'String', ['Normalized depth difference(s) vs. ' handles.core{handles.masterno} ' depth']);
else
    set(handles.axtitle(1), 'String', ['Depth vs. ' handles.core{handles.masterno} ' depth']);
end;

%---

function cursor_Callback(hObject, handles) % Change cursor type from crosshair to fullcrosshair and back
if strcmp(get(handles.fig, 'pointer'), 'crosshair')
    set(handles.fig, 'pointer', 'fullcrosshair');
else
    set(handles.fig, 'pointer', 'crosshair');
end;

%---

function exit_Callback(hObject, handles);
matchmaker('evaluate_Callback',hObject,[],handles, 'closing_evaluate'); % THIS LINE SHOULD BE UNCOMMENTED !!!!
figure(handles.matchmakerfighandle)
delete(handles.fig)

%---

function inc_Callback(hObject, handles);
value = str2num(get(handles.inc, 'string'));
if ~isreal(value) | length(value) ~= 1
    set(handles.inc, 'string', '1');
elseif value < 0
    set(handles.inc, 'string', num2str(abs(value)));
end

%---
    
function xlim_Callback(hObject, handles, lowhigh);
lowmp = str2num(get(handles.lowmp, 'string'));
highmp = str2num(get(handles.highmp, 'string'));
oldmplim = handles.current_mp;
if mod(lowmp,1)==0 & isreal(lowmp) & length(lowmp) == 1 & mod(highmp,1)==0 & isreal(highmp) & length(highmp)==1 % check if input is valid
    for i = 1:handles.N
        mp = handles.mp{i};
        Nmpi134(i) = length(mp(find(mp(:,2)==1 | mp(:,2)==3 | mp(:,2)==4),1));
    end;
    Nmp = min(Nmpi134);
    lowmp = max(1, lowmp);
    lowmp = min(Nmp-1, lowmp);
    highmp = min(Nmp, highmp);
    highmp = max(2, highmp);
    if lowmp>highmp
        if lowhigh == -1
            highmp = min(Nmp, lowmp+(diff(oldmplim)-1));
        else
            lowmp = max(1, highmp-(diff(oldmplim)-1));
        end;
        lowmp = max(1, lowmp);
        lowmp = min(Nmp-1, lowmp);
        highmp = min(Nmp, highmp);
        highmp = max(2, highmp);
    end;
    set(handles.highmp, 'string', num2str(highmp));
    set(handles.lowmp, 'string', num2str(lowmp));
    handles.current_mp = [lowmp highmp];
    guidata(handles.fig, handles);
    mp = handles.mp{handles.masterno};
    mp134 = mp(find(mp(:,2)==1 | mp(:,2)==3 | mp(:,2)==4),1);
    xlim = [mp134(lowmp)-0.1 mp134(highmp)+0.1];
    set(handles.ax, 'xlim', xlim);
    plotcurves(handles);
else
    set(handles.lowmp, 'string', handles.current_mp(1));
    set(handles.highmp, 'string', handles.current_mp(2));
end;
%---

function plotcurves(handles) 
lowmp = str2num(get(handles.lowmp, 'string'));
highmp = str2num(get(handles.highmp, 'string'));
idx1 = [lowmp:highmp];
cla(handles.ax(1));
cla(handles.ax(2));
if length(idx1) < 2
    return
end;
colours = [{'b'} {'g'} {'r'} {'c'} {'m'} {'k'}];

mpmaster = handles.mp{handles.masterno};
mpmaster134 = mpmaster(find(mpmaster(:,2)==1 | mpmaster(:,2)==3 | mpmaster(:,2)==4),:);
mpmaster25 = mpmaster(find((mpmaster(:,2)==2 | mpmaster(:,2)==5) & mpmaster(:,1)>=mpmaster134(lowmp) & mpmaster(:,1)<=mpmaster134(highmp)),:);

for i = setdiff([1:handles.N], handles.masterno);
    mp = handles.mp{i};
    if length(mp(:,1))>=2
        mp134 = mp(find(mp(:,2)==1 | mp(:,2)==3 | mp(:,2)==4),:);
        mp25 = mp(find((mp(:,2)==2 | mp(:,2)==5) & mp(:,1)>=mp134(lowmp) & mp(:,1)<=mp134(highmp)),:);
        idx_13 = intersect(find(mpmaster134(1:idx1(end),2)<4 & mp134(1:idx1(end),2)<4), idx1);
        if length(idx_13) == 1
            disp('Not enough first order matchpoints on screen (blue type 4 matchpoints do not count)');
        else
            if length(idx_13) == 2
                deltadepth{i} = [mpmaster134(idx_13(1), 1) diff(mp134(idx_13, 1)); mp134(idx_13(2), 1) diff(mp134(idx_13, 1))];
            else
                deltadepth{i} = stepit([mpmaster134(idx_13(2:end), 1) diff(mp134(idx_13, 1))./diff(mpmaster134(idx_13, 1))]);
                deltadepth{i}(1,1) = mpmaster134(idx_13(1), 1);
            end;
            plot(deltadepth{i}(:,1), deltadepth{i}(:,2), 'color', colours{i}, 'linewidth', 2, 'parent', handles.ax(2));
%        plot(deltadepth{i}([1 end],1), [1 1]*mean(deltadepth{i}(:,2)), ':', 'color', colours{i}, 'linewidth', 1, 'parent', handles.ax(2));
        
            plotdiff = get(handles.diff, 'Value');
            offset = mean(mp134(idx1, 1)-mpmaster134(idx1, 1));

            plot(mpmaster134(idx1, 1), mp134(idx1, 1)-plotdiff*(mpmaster134(idx1, 1)+offset), '--', 'color', colours{i}, 'Linewidth', 1, 'parent', handles.ax(1), 'hittest', 'off');
            plot(mpmaster134(idx_13, 1), mp134(idx_13, 1)-plotdiff*(mpmaster134(idx_13, 1)+offset), '-', 'color', colours{i}, 'Linewidth', 2, 'parent', handles.ax(1), 'hittest', 'off');
            if size(mpmaster25, 1)==size(mp25, 1)
    % faster !, but without callbacks : plot(mpmaster25, mp25, '-', 'color', colours{i}, 'Linewidth', 0.5, 'parent', handles.ax(1), 'hittest', 'off');
                for j = 1:size(mp25,1)
                    if mp25(j,2) == 2 & mpmaster25(j,2) == 2
                        plot(mpmaster25(j, 1), mp25(j, 1)-plotdiff*(mpmaster25(j, 1)+offset), '.', 'color', colours{i}, 'Markersize', 8, 'parent', handles.ax(1), 'ButtonDownFcn', 'matchmaker_evaluate(''mpclick_Callback'',gcbo,[],guidata(gcbo))', 'Tag', num2str(mpmaster25(j, 1)));
                    else
                        plot(mpmaster25(j, 1), mp25(j, 1)-plotdiff*(mpmaster25(j, 1)+offset), 'o', 'color', colours{i}, 'Markersize', 2, 'parent', handles.ax(1), 'ButtonDownFcn', 'matchmaker_evaluate(''mpclick_Callback'',gcbo,[],guidata(gcbo))', 'Tag', num2str(mpmaster25(j, 1)));
                    end;
                end;
            end;
            for j = 1:length(idx1)
                if (mp134(idx1(j),2) < 1.5) & (mpmaster134(idx1(j),2) < 1.5)
                    plot(mpmaster134(idx1(j), 1), mp134(idx1(j), 1)-plotdiff*(mpmaster134(idx1(j), 1)+offset), '.', 'color', colours{i}, 'Markersize', 20, 'parent', handles.ax(1), 'ButtonDownFcn', 'matchmaker_evaluate(''mpclick_Callback'',gcbo,[],guidata(gcbo))', 'Tag', num2str(idx1(j)));
                else
                    plot(mpmaster134(idx1(j), 1), mp134(idx1(j), 1)-plotdiff*(mpmaster134(idx1(j), 1)+offset), 'o', 'color', colours{i}, 'Markersize', 5, 'parent', handles.ax(1), 'ButtonDownFcn', 'matchmaker_evaluate(''mpclick_Callback'',gcbo,[],guidata(gcbo))', 'Tag', num2str(idx1(j)));
                end;
            end;
        end;
    end;
end;
axes(handles.ax(2))
h_lgd = legend(handles.ax(2), handles.core(setdiff([1:handles.N], handles.masterno)));
set(h_lgd, 'box', 'off');
    

function mpclick_Callback(hObject, handles);
id = get(hObject, 'tag');
set(handles.info, 'string', id);
mpmaster = handles.mp{handles.masterno};
mpmaster134 = mpmaster(find(mpmaster(:,2)==1 | mpmaster(:,2)==3 | mpmaster(:,2)==4),1);
mpmaster25 = mpmaster(find(mpmaster(:,2)==2 | mpmaster(:,2)==5), 1);
pos = str2num(id);
if isfield(handles, 'indicatorline')
    if ishandle(handles.indicatorline)
        delete(handles.indicatorline);
    end;
end;
if ismember(pos, mpmaster25)
    handles.indicatorline = plot([pos pos], get(handles.ax(2), 'ylim')+[1e-3 -1e-3], ':k', 'parent', handles.ax(2));
else
    pos = mpmaster134(pos);
    handles.indicatorline = plot([pos pos], get(handles.ax(2), 'ylim')+[1e-3 -1e-3], ':k', 'parent', handles.ax(2));
end;
guidata(handles.fig, handles);