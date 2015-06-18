function matchmaker(varargin)
% ------------------------------------------------------------------------
% MATCHMAKER version 2.03, 22 April 2010, Sune Olander Rasmussen
% Copyright (C) Sune Olander Rasmussen
% MATCHMAKER syntax : matchmaker(filenames, datafileID, numberofpanels);
%    filenames        Name of .m file containing a list of data and matchpoint file names (try 'files_main')
%    datafileID       Vector of length N indicating which datafiles from "filenames" that should be used
%    numberofpanels   Vector of length N determining the number of data sub-panels in each data window
%
% Example : matchmaker('files_main', [2 3], [1 4]); 
%    opens GRIP (no. 2 in "files_main") with only one data panel, and
%    NGRIP1 (no. 3 in "files_main") with four data sub-panels.
% ------------------------------------------------------------------------
% If the program fails to start and you have problems getting rid of the
% figure window, issue the following commands in the command line :
%    set(0, 'showhiddenhandles', 'on')
%    delete(gcf)
%    set(0, 'showhiddenhandles', 'off')

%try
switch nargin
    case {0 1 2} % Only for development purposes
        delete(gcf);
        clc;
        disp('MATCHMAKER needs input arguments! Take a look at the users'' guide.');
    case 3
        eval('open_request(varargin{1:3})');
    case {4, 5, 6, 7}
        eval([char(varargin{1}) '(varargin{[2 4:nargin]})']);
    otherwise
        disp('Wrong number of arguments when calling matchmaker.m');
 end;
%catch
%    disp(lasterr)
%    dummy = input('ERROR : Log the error and press enter to continue');
%    disp(' ');
%end;

%---

function open_request(datafiles, fileno, NN)
    try
        eval(datafiles);
    catch
        disp('ERROR in Matchmaker : File name file could not be read.');
        return;
    end;
    load matchmaker_sett
    if size(fileno) ~= size(NN)
        disp('ERROR in Matchmaker : 2nd and 3rd argument must have same dimension.');
        return;
    end;        
    scrsize = get(0, 'screensize');
    hgt = scrsize(4);
    handles.fig = figure('position', [1 60 scrsize(3) scrsize(4)-100], 'name', 'Matchmaker', 'CloseRequestFcn', 'matchmaker(''exit_Callback'',gcbo,[],guidata(gcbo))', 'KeyPressFcn', 'matchmaker(''keypressed_Callback'',gcbo,[],guidata(gcbo))', 'nextplot', 'add', 'color', 0.9*[1 1 1], 'pointer', 'crosshair', 'Interruptible', 'off', 'Numbertitle', 'off', 'tag', 'matchmakermainwindow', 'integerhandle', 'off');
    
    N = length(fileno);
    handles.fileno = fileno;
    handles.N = N;
    handles.NN = NN;
    h_wait = waitbar(0.05, 'Please be patient ... loading'); % Wait window is launched
    for i = 1:N
        waitbar(0.2*i, h_wait, ['Please be patient ... loading data file no. ' num2str(i)]); % Wait window is launched
        try
            load(['Subroutines' filesep 'matchmaker' filesep 'data' filesep files.datafile{fileno(i)}]);
        catch
            disp(['Data file of record ' num2str(fileno(i)) ' not found']);
            delete(handles.fig);
            delete(h_wait);
            return
        end;
        handles.data{i} = data;
        handles.depth{i} = depth;
        handles.depth_no{i} = depth_no;
        handles.species{i} = species;
        handles.colours{i} = colours;
        waitbar(0.2*i+0.1, h_wait, ['Please be patient ... loading matchpoint file no. ' num2str(i)]); % Wait window is updated
        try
            load(files.matchfile{fileno(i)});
        catch
            disp(['Matchpoint file ' files.matchfile{fileno(i)} ' not found']);
            delete(handles.fig);
            delete(h_wait);
            return
        end;
        handles.matchfile{i} = files.matchfile{fileno(i)};
        handles.core{i} = files.core{fileno(i)};
        handles.mp{i} = mp;
        if length(sett.specs{fileno(i)}) == NN(i)
            handles.selectedspecs{i} = sett.specs{fileno(i)};
        else
            if length(sett.specs{fileno(i)}) < NN(i)
                handles.selectedspecs{i} = [sett.specs{fileno(i)} (length(sett.specs{fileno(i)})+1):NN(i)];
            else
                handles.selectedspecs{i} = sett.specs{fileno(i)}(1:NN(i));
            end;
        end;

%       if max(sett.specs{fileno(i)}) > NN(i) CHANGED 22/4 2007 to solve
%       problem that resets specis to 1:5 all the time

%       Version replaced 22/4 2010 to solve problem when switching to a
%       core with less available species
%       if max(sett.specs{fileno(i)}) > length(depth_no)
%           handles.selectedspecs{i} = 1:NN(i); disp('bing'); disp(sett.specs{fileno(i)}); disp(NN(i));
%       end;
        handles.selectedspecs{i} = min(handles.selectedspecs{i}, length(depth_no));
        
    end;
    close(h_wait)
        
    font1 = 8;
    font2 = 10;
    
    y0 = 0.03;
    y1 = 0.02;
    dy = 0.015;
    yh = (1-y0-(N+2)*dy)/N;
    yoverlap = 0.45; % fraction
    yhsub = (yh-y1)./(NN-(NN-1)*yoverlap);
    
    x0 = 0.05;
    x1 = 0.08;
    dx = 0.005;
    xw = 1-x0-2*x1-4*dx;
    
    eh = 0.029*710/hgt;
    dh = 0.002*710/hgt;
    
    dummyax = axes('position', [0 0 1 1], 'xlim', [0 1], 'ylim', [0 1], 'visible', 'off', 'nextplot', 'add');

    %handles.title = text(dx, 0.4*y0, 'MATCHMAKER', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Fontsize', font2, 'fontweight', 'bold', 'parent', dummyax);
    handles.save = uicontrol('units', 'normalized', 'position', [x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Save', 'style', 'pushbutton', 'callback', ['matchmaker(''save_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center', 'enable', 'off');
    handles.accordianize = uicontrol('units', 'normalized', 'position', [3*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Accordianize', 'style', 'pushbutton', 'callback', ['matchmaker(''accordianize_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
    handles.masterno = uicontrol('units', 'normalized', 'position', [5*x0+x1+3*dx 0.4*y0 0.6*x0 eh], 'string', '1', 'style', 'edit', 'callback', ['matchmaker(''masterno_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
    handles.evaluate = uicontrol('units', 'normalized', 'position', [6*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Evaluate', 'style', 'pushbutton', 'callback', ['matchmaker(''evaluate_Callback'',gcbo,[],guidata(gcbo), ''button'')'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
    handles.mark = uicontrol('units', 'normalized', 'position', [8*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Mark ?', 'style', 'togglebutton', 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center', 'value', 0);
    handles.dummymark = uicontrol('units', 'normalized', 'position', [10*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Dummy ?', 'style', 'togglebutton', 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center', 'value', 0);
    handles.plotmp2 = uicontrol('units', 'normalized', 'position', [12*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', '2nd order ?', 'style', 'togglebutton', 'callback', ['matchmaker(''plotmp2_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center', 'value', 1, 'Selected', 'off');
    handles.check = uicontrol('units', 'normalized', 'position', [14*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Check mps', 'style', 'pushbutton', 'callback', ['matchmaker(''check_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
    handles.exit = uicontrol('units', 'normalized', 'position', [16*x0+x1+3*dx 0.4*y0 1.6*x0 eh], 'string', 'Exit', 'style', 'pushbutton', 'callback', ['matchmaker(''exit_Callback'',gcbo,[],guidata(gcbo))'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'bold', 'horizontalalignment', 'center');
        
    for i = 1:N
        handles.bigax(i) = axes('position', [x0+x1+3*dx y0+2*dy+y1+(i-1)*(yh+dy) xw yh-y1], 'nextplot', 'add', 'ylim', [0 1], 'ytick', [], 'box', 'off', 'fontsize', font1, 'ycolor', 'w', 'xcolor', 'w', 'ButtonDownFcn', ['matchmaker(''axesclick_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ')']);
        handles.bigax2(i) = axes('position', [x0+x1+3*dx y0+2*dy+y1+(i-1)*(yh+dy) xw yh-y1], 'nextplot', 'add', 'visible', 'off', 'hittest', 'off', 'clipping', 'on', 'ylim', [0 1]);
        handles.name(i) = text(dx+x0/2, y0+dy+i*(yh+dy), files.core{fileno(i)}, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Fontsize', font2, 'fontweight', 'bold', 'parent', dummyax,'interpreter','none');
        handles.minx(i) = uicontrol('units', 'normalized', 'position', [dx y0+dy+i*(yh+dy)-2*(eh+dh) x0 eh], 'string', '0', 'style', 'edit', 'callback', ['matchmaker(''xscale_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ', 1)'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'right');
        handles.maxx(i) = uicontrol('units', 'normalized', 'position', [dx y0+dy+i*(yh+dy)-3*(eh+dh) x0 eh], 'string', '1', 'style', 'edit', 'callback', ['matchmaker(''xscale_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ', 2)'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'right');
        handles.back(i) = uicontrol('units', 'normalized', 'position', [dx y0+dy+i*(yh+dy)-5*(eh+dh) x0 eh], 'string', '<--', 'style', 'pushbutton', 'callback', ['matchmaker(''move_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ', -1)'], 'fontname', 'default', 'fontsize', font2, 'fontweight', 'bold', 'horizontalalignment', 'center');
        handles.fwd(i) = uicontrol('units', 'normalized', 'position', [dx y0+dy+i*(yh+dy)-6*(eh+dh) x0 eh], 'string', '-->', 'style', 'pushbutton', 'callback', ['matchmaker(''move_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ', 1)'], 'fontname', 'default', 'fontsize', font2, 'fontweight', 'bold', 'horizontalalignment', 'center');
        handles.incx(i) = uicontrol('units', 'normalized', 'position', [dx y0+dy+i*(yh+dy)-7*(eh+dh) x0 eh], 'string', '1', 'style', 'edit', 'callback', ['matchmaker(''incx_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ')'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'right');
    end;
    sides = [{'right'} {'left'}];
    for i = 1:N
        handles.tickax{i,1} = axes('position', [x0+x1+3*dx y0+2*dy+y1+(i-1)*(yh+dy) xw yhsub(i)], 'nextplot', 'add', 'color', 'none', 'xcolor', 'k', 'ylim', [0 1], 'box', 'off', 'fontsize', font1, 'yaxislocation', 'left', 'xaxislocation', 'bottom', 'hittest', 'off', 'fontweight', 'bold');
        handles.plotax{i,1} = axes('position', [x0+x1+3*dx y0+2*dy+y1+(i-1)*(yh+dy) xw yhsub(i)], 'nextplot', 'replacechildren', 'visible', 'off', 'hittest', 'off');
        for j = 1:NN(i)
            if j > 1
                handles.tickax{i,j} = axes('position', [x0+x1+3*dx y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i) xw yhsub(i)], 'nextplot', 'add', 'ycolor', handles.colours{i}(handles.selectedspecs{i}(j),:), 'color', 'none', 'xcolor', 'w', 'xtick', [], 'ylim', [0 1], 'box', 'off', 'fontsize', font1, 'yaxislocation', sides{mod(j,2)+1}, 'xaxislocation', 'top', 'hittest', 'off', 'fontweight', 'bold');
                handles.plotax{i,j} = axes('position', [x0+x1+3*dx y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i) xw yhsub(i)], 'nextplot', 'replacechildren', 'visible', 'off', 'hittest', 'off');
            end;
            handles.spec{i,j}  = uicontrol('units', 'normalized', 'position', [(mod(j,2)==1)*(x0+3*dx)+(mod(j,2)==0)*(1-x0-dx) y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i)+yhsub(i)-(eh-2*dh) x0 eh], 'string', handles.species{i}, 'value', handles.selectedspecs{i}(j), 'style', 'popupmenu', 'callback', ['matchmaker(''spec_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ',' num2str(j) ')'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'center');
            handles.offset{i,j} = uicontrol('units', 'normalized', 'position', [(mod(j,2)==1)*(x0+3*dx)+(mod(j,2)==0)*(1-x0-dx) y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i)+yhsub(i)-2*(eh+dh) x0 eh], 'string', 0, 'style', 'edit', 'callback', ['matchmaker(''offset_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ',' num2str(j) ')'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'right');
            handles.autoy{i,j} = uicontrol('units', 'normalized', 'position', [(mod(j,2)==1)*(x0+3*dx)+(mod(j,2)==0)*(1-x0-dx) y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i)+yhsub(i)-3*(eh+dh) (x0-dh)/2 eh], 'string', 'Aut', 'style', 'togglebutton', 'callback', ['matchmaker(''autoy_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ',' num2str(j) ')'], 'value', 1, 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'center');
            handles.logy{i,j}  = uicontrol('units', 'normalized', 'position', [(mod(j,2)==1)*(x0+3*dx)+(mod(j,2)==0)*(1-x0-dx)+(x0+dh)/2 y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i)+yhsub(i)-3*(eh+dh) (x0-dh)/2 eh], 'string', 'Log', 'style', 'togglebutton', 'callback', ['matchmaker(''logy_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ',' num2str(j) ')'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'center');
            handles.miny{i,j} = uicontrol('units', 'normalized', 'position', [(mod(j,2)==1)*(x0+3*dx)+(mod(j,2)==0)*(1-x0-dx) y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i)+yhsub(i)-4*(eh+dh) x0 eh], 'string', 0, 'style', 'edit', 'callback', ['matchmaker(''yscale_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ',' num2str(j) ', 1)'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'right');
            handles.maxy{i,j} = uicontrol('units', 'normalized', 'position', [(mod(j,2)==1)*(x0+3*dx)+(mod(j,2)==0)*(1-x0-dx) y0+2*dy+y1+(i-1)*(yh+dy)+(j-1)*(1-yoverlap)*yhsub(i)+yhsub(i)-5*(eh+dh) x0 eh], 'string', 1, 'style', 'edit', 'callback', ['matchmaker(''yscale_Callback'',gcbo,[],guidata(gcbo),' num2str(i) ',' num2str(j) ', 2)'], 'fontname', 'default', 'fontsize', font1, 'fontweight', 'normal', 'horizontalalignment', 'right');
        end;
        set([handles.tickax{i,:} handles.plotax{i,:} handles.bigax(i) handles.bigax2(i)], 'xlim', sett.xlim(fileno(i),:)); % Set the x limits according to the last values
        set(handles.minx(i), 'string', num2str(sett.xlim(fileno(i),1)));
        set(handles.maxx(i), 'string', num2str(sett.xlim(fileno(i),2)));
    end;
    for i = 1:N
        axes(handles.bigax2(i));
    end;
    handles.storeidx = [];
    guidata(handles.fig, handles);
    for i = 1:N
        for j = 1:NN(i)
            axes(handles.plotax{i,j});
            plotcurve(handles, i, j);
        end;
        handles = plotmp(handles, i);
        update_yminmax(handles.fig, handles, i, 0); % Update the y-scaling
    end;
    set(handles.fig, 'HandleVisibility', 'callback');
    guidata(handles.fig, handles);
    warning('off', 'MATLAB:Axes:NegativeDataInLogAxis');
    
%-----------------------

function varargout = move_Callback(hObject, handles, no, direc);
incX = str2num(get(handles.incx(no), 'string')); % Get increment value. Validity check is done elsewhere
oldlimX = get(handles.bigax(no), 'xlim'); % Get depth limits
set([handles.bigax(no) handles.bigax2(no) handles.plotax{no,:} handles.tickax{no,:}], 'xlim', oldlimX+direc*incX); % Add/Subtract increment value ...
set(handles.minx(no), 'String', num2str(oldlimX(1)+direc*incX)); % and update the edit-boxes
set(handles.maxx(no), 'String', num2str(oldlimX(2)+direc*incX));
handles = plotmp(handles, no);
plotcurve(handles, no, 0);
update_yminmax(hObject, handles, no, 0); % The view has changed, so update the y-sacle edit-boxes
guidata(handles.fig, handles)
if nargout == 1
    varargout{1} = handles;
end;

%---

function offset_Callback(hObject, handles, no1, no2)
value = str2num(get(handles.offset{no1, no2}, 'string')); % get input
if ~isreal(value) | length(value) ~= 1 % check if input is valid
    set(handles.offset{no1, no2}, 'string', '0');
end;
plotcurve(handles, no1, no2);
update_yminmax(hObject, handles, no1, no2); % The view has changed, so update the y-sacle edit-boxes

%---

function masterno_Callback(hObject, handles)
value = str2num(get(handles.masterno, 'string')); % get input
if ~isreal(value) | length(value) ~= 1 | mod(value,1)~=0 | value>handles.N | value<1% check if input is valid
    set(handles.masterno, 'string', '1');
    value = 1;
end;
idx1 = handles.mp1_depth{value};
if length(idx1)<2
    set(handles.evaluate, 'Enable', 'off');
else
    set(handles.evaluate, 'Enable', 'on');
end;

%---

function incx_Callback(hObject, handles, no) % Depth increment edit-box callback. The validity of the input is checked, the actual value is used elsewhere
value = str2num(get(handles.incx(no), 'string')); % get input
if ~isreal(value) | length(value) ~= 1 % check if input is valid
    set(handles.incx(no), 'string', '1');
end;
if value <= 0
    set(handles.incx(no), 'string', num2str(-value)); % Negative values are not logically acceptable
end;

%---

function xscale_Callback(hObject, handles, no, limit) % X-scaling edit-boxes callback
if limit == 1
    value = get(handles.minx(no), 'string'); % Get input
else;
    value = get(handles.maxx(no), 'string'); % Get input
end;
if length(value)>0 & value(1) == 'b'
    if limit == 1
        numvalue = (str2num(value(2:end))-1)*0.55;
        set(handles.minx(no), 'string', num2str(numvalue)); % ... use old lower limit
    else
        numvalue = str2num(value(2:end))*0.55;
        set(handles.maxx(no), 'string', num2str(numvalue)); % ... use old upper limit
    end;
else
    numvalue = str2num(value);
end;
Xlim = get(handles.bigax(no), 'xlim'); % Get old limits
if ~isreal(numvalue) | length(numvalue) ~= 1 % If input is not a valid number ...
    set(handles.minx(no), 'string', num2str(Xlim(1))); % ... use old lower limit
    set(handles.maxx(no), 'string', num2str(Xlim(2))); % ... use old upper limit
    return
end;
if limit == 1
    if numvalue >= Xlim(2) % Check if low < high limit
        set([handles.bigax(no) handles.bigax2(no) handles.plotax{no,:} handles.tickax{no,:}], 'xlim', [numvalue, numvalue + Xlim(2)-Xlim(1)]); % If so, move high limit up as well
        set(handles.maxx(no), 'string', num2str(numvalue + Xlim(2)-Xlim(1))); % and update the corresponding edit-box
    else
        set([handles.bigax(no) handles.bigax2(no) handles.plotax{no,:} handles.tickax{no,:}], 'xlim', [numvalue, Xlim(2)]); % If everything is OK, use the input value
    end;
else
    if numvalue <= Xlim(1)
        set([handles.bigax(no) handles.bigax2(no) handles.plotax{no,:} handles.tickax{no,:}], 'xlim', [numvalue - Xlim(2) + Xlim(1), numvalue]);
        set(handles.minx(no), 'string', num2str(numvalue - Xlim(2) + Xlim(1)));
    else
        set([handles.bigax(no) handles.bigax2(no) handles.plotax{no,:} handles.tickax{no,:}], 'xlim', [Xlim(1), numvalue]);
    end;
end;
plotcurve(handles, no, 0);
handles = plotmp(handles, no);
update_yminmax(hObject, handles, no, 0); % Update the y-scaling
guidata(handles.fig, handles)

%---

function yscale_Callback(hObject, handles, no1, no2, limit) % Y-scaling edit-boxes callback
if limit == 1
    value = str2num(get(handles.miny{no1, no2}, 'string')); % Get input
else;
    value = str2num(get(handles.maxy{no1, no2}, 'string')); % Get input
end;
Ylim = get(handles.tickax{no1, no2}, 'ylim'); % Get old limits
if ~isreal(value) | length(value) ~= 1 % If input is not a valid number ...
    set(handles.miny{no1, no2}, 'string', num2str(Ylim(1))); % ... use old lower limit
    set(handles.maxy{no1, no2}, 'string', num2str(Ylim(2))); % ... use old upper limit
    return
end;
if limit == 1
    if value >= Ylim(2) % Check if low < high limit
        set([handles.plotax{no1, no2} handles.tickax{no1, no2}], 'ylim', [value, value + Ylim(2)-Ylim(1)]); % If so, move high limit up as well
        set(handles.maxy{no1, no2}, 'string', num2str(value + Ylim(2)-Ylim(1))); % and update the corresponding edit-box
    else
        set([handles.plotax{no1, no2} handles.tickax{no1, no2}], 'ylim', [value, Ylim(2)]); % If everything is OK, use the input value
    end;
else
    if value <= Ylim(1)
        set([handles.plotax{no1, no2} handles.tickax{no1, no2}], 'ylim', [value - Ylim(2) + Ylim(1), value]);
        set(handles.miny{no1, no2}, 'string', num2str(value - Ylim(2) + Ylim(1)));
    else
        set([handles.plotax{no1, no2} handles.tickax{no1, no2}], 'ylim', [Ylim(1), value]);
    end;
end;
set(handles.autoy{no1, no2}, 'Value', 0);

%---

function autoy_Callback(hObject, handles, no1, no2);
if get(handles.autoy{no1, no2}, 'value') == 1
    set(handles.plotax{no1, no2}, 'ylimmode', 'auto');
    update_yminmax(hObject, handles, no1, no2);
else
    set(handles.plotax{no1, no2}, 'ylimmode', 'manual');
end;
set(handles.tickax{no1, no2}, 'ylim', get(handles.plotax{no1, no2}, 'ylim'));

%---

function logy_Callback(hObject, handles, no1, no2);
if get(handles.logy{no1, no2}, 'value') == 1
    set([handles.plotax{no1, no2} handles.tickax{no1, no2}], 'yscale', 'log');
else
    set([handles.plotax{no1, no2} handles.tickax{no1, no2}], 'yscale', 'linear');
end;
set(handles.tickax{no1, no2}, 'ylim', get(handles.plotax{no1, no2}, 'ylim'));
update_yminmax(hObject, handles, no1, no2);

%---

function spec_Callback(hObject, handles, no1, no2);
handles.selectedspecs{no1}(no2) = get(handles.spec{no1, no2}, 'Value');
plotcurve(handles, no1, no2);
set(handles.autoy{no1, no2}, 'Value', 1);
autoy_Callback(hObject, handles, no1, no2);
update_yminmax(hObject, handles, no1, no2);
guidata(handles.fig, handles)

%---

function axesclick_Callback(hObject, handles, no1);
if get(handles.mark, 'Value') == 1
    pos = get(handles.bigax(no1), 'currentpoint');
    pos = 0.001*round(1000*pos(1,1));
    but = get(handles.fig, 'selectiontype');
    dummy = get(handles.dummymark, 'Value');
    if strcmp(but, 'normal')
        if dummy
            mptype = 4;
        else
            mptype = 1;
        end;
    elseif strcmp(but, 'alt')
        if dummy
            mptype = 5;
        else
            mptype = 2;
        end;
        if get(handles.plotmp2, 'value') == 0
            return
        end;
    elseif strcmp(but, 'extend')
        mptype = 3;
    else
        return;
    end
    mp = handles.mp{no1};
    mp = [mp; pos mptype];
    handles.mp{no1} = sortrows(mp);
    handles = plotmp(handles, no1);
    set(handles.save, 'enable', 'on');
    guidata(hObject, handles);
end;

%---

function mpclick_Callback(hObject, handles, no1);
if get(handles.mark, 'Value') == 1
    pos = get(hObject, 'xdata');
    pos = pos(1);
    mp = handles.mp{no1};
    delindx = find(mp(:,1)==pos);
    mp = mp(setdiff(1:length(mp(:,1)), delindx),:);
    handles.mp{no1} = mp;
    handles = plotmp(handles, no1);
    set(handles.save, 'enable', 'on');
    guidata(handles.fig, handles);
end;

%---

function update_yminmax(hObject, handles, no1, no2) % Updates the minY and maxY edit-boxes with the current Y-axis limits. Is called whenever the view changes.
if no2 == 0;
    no2 = 1:handles.NN(no1);
end;
for j = 1:length(no2)
    limY = get(handles.plotax{no1, no2(j)}, 'ylim');
    set(handles.miny{no1, no2(j)}, 'string', num2str(limY(1)));
    set(handles.maxy{no1, no2(j)}, 'string', num2str(limY(2)));
    set(handles.tickax{no1, no2(j)}, 'ylim', limY);
end;

%---

function plotcurve(handles, no1, no2) 
if no2 == 0;
    no2 = 1:handles.NN(no1);
end;
xlim = get(handles.bigax(no1), 'xlim');
for j = no2
    offset = str2num(get(handles.offset{no1, j}, 'String'));
    offsetxlim = xlim - offset;
    depth = handles.depth{no1}{handles.depth_no{no1}(handles.selectedspecs{no1}(j))};
    data = handles.data{no1}{handles.selectedspecs{no1}(j)};
    idx = [max(1, min(find(depth>=offsetxlim(1)))-1) : min(length(data), max(find(depth<=offsetxlim(2)))+1)];
    if isempty(idx)
        plot(0, 0, 'parent', handles.plotax{no1,j}, 'hittest', 'off');
    else
        plotdepth = depth(idx);
        plotdata = data(idx);
        plot(offset+plotdepth, plotdata, 'color', handles.colours{no1}(handles.selectedspecs{no1}(j),:), 'linewidth', 1.5, 'parent', handles.plotax{no1,j}, 'hittest', 'off');
    end;
    set([handles.spec{no1, j} handles.offset{no1, j} handles.autoy{no1, j} handles.logy{no1, j} handles.miny{no1, j} handles.maxy{no1, j}], 'Backgroundcolor', 'w', 'Foregroundcolor', handles.colours{no1}(handles.selectedspecs{no1}(j),:));
    set(handles.tickax{no1, j}, 'ycolor', handles.colours{no1}(handles.selectedspecs{no1}(j),:)); 
end;
%guidata(handles.fig, handles); % MAYBE THIS SHOULD BE INCLUDED ????

%---

function handles = plotmp(handles, no1)
cla(handles.bigax2(no1));
greytone = 0.85*[1 1 1];
redtone = 0.85*[1 0 0];
bluetone = 0.85*[0 0 1];
mp = handles.mp{no1};
xlim = get(handles.bigax(no1), 'xlim');
if isempty(mp)
    return
else
    mp1 = mp(find(mp(:,2)==1),1);
    mp3 = mp(find(mp(:,2)==3),1);
    mp4 = mp(find(mp(:,2)==4),1);
    mp134 = mp(find(mp(:,2)==1 | mp(:,2)==3| mp(:,2)==4),1);
    if ~isempty([mp1' mp3' mp4'])
        idx1 = find(mp1>=xlim(1) & mp1<=xlim(2));
        idx3 = find(mp3>=xlim(1) & mp3<=xlim(2));
        idx4 = find(mp4>=xlim(1) & mp4<=xlim(2));
        idx134 = find(mp134>=xlim(1) & mp134<=xlim(2));
        if ~isempty(idx1)
            plot((mp1(idx1)*[1 1])', repmat([0.01 0.93]', 1, length(idx1)), 'linewidth', 6, 'color', greytone, 'parent', handles.bigax2(no1), 'ButtonDownFcn', ['matchmaker(''mpclick_Callback'',gcbo,[],guidata(gcbo),' num2str(no1) ')']);
        end;
        if ~isempty(idx3)
            plot((mp3(idx3)*[1 1])', repmat([0.01 0.93]', 1, length(idx3)), 'linewidth', 6, 'color', redtone, 'parent', handles.bigax2(no1), 'ButtonDownFcn', ['matchmaker(''mpclick_Callback'',gcbo,[],guidata(gcbo),' num2str(no1) ')']);
        end;
        if ~isempty(idx4)
            plot((mp4(idx4)*[1 1])', repmat([0.01 0.93]', 1, length(idx4)), 'linewidth', 6, 'color', bluetone, 'parent', handles.bigax2(no1), 'ButtonDownFcn', ['matchmaker(''mpclick_Callback'',gcbo,[],guidata(gcbo),' num2str(no1) ')']);
        end;
        text(mp134(idx134), 0.93*ones(length(idx134),1), num2str(idx134), 'parent', handles.bigax2(no1), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'fontsize', get(handles.tickax{1,1}, 'fontsize'), 'color', 'k');
        handles.mp1_idx{no1} = idx134;
        handles.mp1_depth{no1} = mp134(idx134);
    else
        handles.mp1_idx{no1} = 0;
        handles.mp1_depth{no1} = [];
    end;
    idx2 = find(mp(:,2)==2 & mp(:,1)>=xlim(1) & mp(:,1)<=xlim(2));
    idx5 = find(mp(:,2)==5 & mp(:,1)>=xlim(1) & mp(:,1)<=xlim(2));
    if ~isempty(idx2) & get(handles.plotmp2, 'Value') == 1
        plot((mp(idx2,1)*[1 1])', repmat([0.05 0.89]', 1, length(idx2)), 'linewidth', 4, 'color', greytone, 'parent', handles.bigax2(no1), 'ButtonDownFcn', ['matchmaker(''mpclick_Callback'',gcbo,[],guidata(gcbo),' num2str(no1) ')']);
    end;
    if ~isempty(idx5) & get(handles.plotmp2, 'Value') == 1
        plot((mp(idx5,1)*[1 1])', repmat([0.05 0.89]', 1, length(idx5)), 'linewidth', 4, 'color', bluetone, 'parent', handles.bigax2(no1), 'ButtonDownFcn', ['matchmaker(''mpclick_Callback'',gcbo,[],guidata(gcbo),' num2str(no1) ')']);
    end;
end;

if length(handles.mp1_idx{str2num(get(handles.masterno, 'string'))})<2
    set(handles.evaluate, 'Enable', 'off');
else
    set(handles.evaluate, 'Enable', 'on');
end;
%---

function keypressed_Callback(hObject, handles) % Translate keypress to appropriate button actions.
key = double(get(handles.fig, 'currentcharacter'));
if ~isempty(key)
    switch key
    case 28   %<-
        handles = move_Callback(hObject, handles, str2num(get(handles.masterno, 'string')), -1);
        accordianize_Callback(hObject, handles);
    case 29  %->
        handles = move_Callback(hObject, handles, str2num(get(handles.masterno, 'string')), 1);
        accordianize_Callback(hObject, handles);
    case {80, 112}   %p, P
        set(handles.fig, 'InvertHardcopy', 'on', 'paperunits', 'centimeters', 'paperorientation', 'landscape', 'papertype', 'A4', 'paperposition', [1 1 27.7 19], 'renderer', 'painters');
%        print(handles.fig, '-dpsc2', '-r300', 'matchmaker.eps');
        print(handles.fig, '-dpsc2', '-noui', '-r300', 'matchmaker.eps');
    case {97, 65}   %a, A
        accordianize_Callback(hObject, handles);
    case {99, 67}   %c, C
        cursor_Callback(hObject, handles);
    case {109, 77}  %m, M
        set(handles.mark, 'value', get(handles.mark, 'value')==0);
    case {115, 83}  %s, S
        if strcmp(get(handles.save, 'enable'), 'on');
            save_Callback(hObject, handles);
        end;
    case {120, 88}  %x, X
        exit_Callback(hObject, handles);
    case 50  %2
        set(handles.plotmp2, 'value', get(handles.plotmp2, 'value')==0);
        plotmp2_Callback(hObject, handles);
    otherwise % If key not defined, show info window
         disp(['Just for information : MATCHMAKER undefined key callback : ' num2str(key)]);
         h_help = helpdlg({...
             'Available keyboard commands:';
             '<- = move selected core one frame back and accordianize'
             '-> = move selected core one frame forward and accordianize'
             'A  = Accordianize'
             'C  = change Cursor type'
             'M  = Toggle "Mark ?" button on/off'
             'S  = Save matchpoint files'
             'X  = eXit'
             '2  = Toggle "2nd order" button on/off'});
    end;
end;

%---

function cursor_Callback(hObject, handles) % Change cursor type from crosshair to fullcrosshair and back
if strcmp(get(handles.fig, 'pointer'), 'crosshair')
    set(gcf, 'pointer', 'fullcrosshair');
else
    set(gcf, 'pointer', 'crosshair');
end;

%---

function accordianize_Callback(hObject, handles);
masterno = str2num(get(handles.masterno, 'string'));
mastermp = handles.mp1_idx{masterno};
if isempty(mastermp)
    return
end;
xlim = get(handles.bigax(masterno), 'xlim');
mp_depth = handles.mp1_depth{masterno}([1 end]);
frac = (mp_depth-xlim(1))./(xlim(2)-xlim(1));
if length(mastermp) == 1
    for i = setdiff(1:handles.N, masterno)
        mpi = handles.mp{i};
        if length(mpi)>=max(mastermp)
            mpi134 = mpi(find(mpi(:,2)==1 | mpi(:,2)==3 | mpi(:,2)==4),1);
            newxlim = round(1000*(mpi134(mastermp)+[-frac(1) (1-frac(1))]*(xlim(2)-xlim(1))))/1000;
            set([handles.bigax(i) handles.bigax2(i) handles.plotax{i,:} handles.tickax{i,:}], 'xlim', newxlim);
            plotcurve(handles, i, 0);
            set(handles.minx(i), 'string', num2str(newxlim(1)));
            set(handles.maxx(i), 'string', num2str(newxlim(2)));
            handles = plotmp(handles, i);
        end;
    end;
else
    for i = setdiff(1:handles.N, masterno)
        mpi = handles.mp{i};
        if length(mpi)>=max(mastermp)
            mpi134 = mpi(find(mpi(:,2)==1 | mpi(:,2)==3 | mpi(:,2)==4),1);
            mpi134_depth = mpi134(mastermp([1 end]));
            newwidth = diff(mpi134_depth)/(frac(2)-frac(1));
            newxlim = round(100*(mpi134_depth(1)-frac(1)*newwidth + [0 newwidth]))/100;
            set([handles.bigax(i) handles.bigax2(i) handles.plotax{i,:} handles.tickax{i,:}], 'xlim', newxlim);
            plotcurve(handles, i, 0);
            handles = plotmp(handles, i);
            set(handles.minx(i), 'string', num2str(newxlim(1)));
            set(handles.maxx(i), 'string', num2str(newxlim(2)));
        end;
    end;
end;
guidata(handles.fig, handles); % Is it OK to have it outside the loop ???

%---

function evaluate_Callback(hObject, handles, identify);
if strcmp(identify, 'button')
    masterno = str2num(get(handles.masterno, 'string'));
    if isfield(handles, 'evaluatefigurehandle')
        matchmaker_evaluate('evalreuse', handles.evaluatefigurehandle, [], handles.mp, handles.core, masterno, handles.mp1_idx{masterno});
    else
        matchmaker_evaluate('evalopen', handles.fig, [], handles.mp, handles.core, masterno, handles.mp1_idx{masterno});
    end;
elseif strcmp(identify, 'opening_evaluate')
    evalhandles = handles; % the input argument 'handles' is the handles of the evaluate window. 
    handles = guidata(evalhandles.matchmakerfighandle); %This line, redefines 'handles' to contain the handles of the matchmaker window
    handles.evaluatefigurehandle = evalhandles.fig;
    guidata(handles.fig, handles);
elseif strcmp(identify, 'closing_evaluate')
    evalhandles = handles; % the input argument 'handles' is the handles of the evaluate window. 
    handles = guidata(evalhandles.matchmakerfighandle); %This line, redefines 'handles' to contain the handles of the matchmaker window
    handles = rmfield(handles, 'evaluatefigurehandle');
    guidata(handles.fig, handles);
else
   disp(['Error in MATCHMAKER.m, evaluate_Callback : ' identify]);
end;
%mp = handles.mp; core = handles.core; mp1 = handles.mp1_idx{masterno}; save('testevaldata.mat', 'mp', 'core', 'masterno', 'mp1');

%---

function check_Callback(hObject, handles);
clc
names = [];
for i = 1:handles.N
    names = [names char(handles.core{i}) ' '];
end;
disp(['MATCHMAKER DIAGNOSTICS            for matchpoints from the cores : ' names]);
disp(' ');
for i = 1:handles.N
    mpi = handles.mp{i};
    mpi1 = mpi(find(mpi(:,2)==1),1);
    mpi25 = mpi(find(mpi(:,2)==2 | mpi(:,2)==5),1);
    mpi134 = mpi(find(mpi(:,2)==1 | mpi(:,2)==3 | mpi(:,2)==4),1);
    Nmp1(i) = length(mpi1);
    Nmp134(i) = length(mpi134);
    for j = 1:Nmp134(i)-1
        Nmp25(i,j) = sum(mpi25>=mpi134(j) & mpi25<mpi134(j+1));
    end;
    if length(mpi134)<2
        mindist134(i) = 0;
    else
        mindist134(i) = min(diff(mpi134));
    end;
    if length(mpi25)<2
        mindist25(i) = 0;
    else
        mindist25(i) = min(diff(mpi25));
    end;
end;
if isempty(setdiff(Nmp1, Nmp1(1)))
    disp('The number of type 1 matchpoints is the same for all cores.');
else
    disp(['The number of type 1 matchpoints is not the same for all cores   : ' num2str(Nmp1, 2)]);
end;
if isempty(setdiff(Nmp134, Nmp134(1)))
    disp('The total number of type 1/3/4 matchpoints is the same for all cores.');
else
    disp(['The total number of type 1/3/4 matchpoints is not the same for all cores   : ' num2str(Nmp134, 2)]);
end;

disp(' ')
diffidx = [];
for j = 1:max(Nmp134)-1
    if ~isempty(setdiff(Nmp25(:,j), Nmp25(1,j)))
        diffidx = [diffidx j];
    end;
end;
if isempty(diffidx)
    disp('The number of type 2/5 matchpoints between each set of type 1/3/4 matchpoints is the same for all cores.');
else
    disp('The number of type 2/5 matchpoints between each set of type 1/3/4 matchpoints is not the same for all cores :');
    for k = 1:length(diffidx)
        disp(['Number of type 2/5 matchpoints between type 1/3/4 matchpoint ' num2str(diffidx(k)) ' and ' num2str(diffidx(k)+1) ' : ' num2str(Nmp25(:, diffidx(k))', 2)]);
    end;
end;
disp(' ');
disp(['The minimum distance between two adjacent type 1/3/4 matchpoints is  : ' num2str(mindist134, 2)]);
disp(['The minimum distance between two adjacent type 2/5 matchpoints is  : ' num2str(mindist25, 2)]);
disp('(rounded off to nearest centimeter value)');
disp(' ');
answer = input('(Press ENTER to return)');
disp(' ');
disp(' ');
figure(handles.fig);

%---

function save_Callback(hObject, handles);
for i = 1:handles.N
    status = copyfile(handles.matchfile{i},[handles.matchfile{i} '.backup']);
    if status == 0
        disp(['MATCHMAKER.m warning : Could not back up matchpoint file ' handles.matchfile{i} ' before saving']);
        disp('   (this error does not influence the saving procedure itself)');
    end;
    mp = handles.mp{i};
    save(handles.matchfile{i}, 'mp', '-MAT');
end;
set(handles.save, 'enable', 'off');

%---

function plotmp2_Callback(hObject, handles);
for i = 1:handles.N
    handles = plotmp(handles, i);
end;
guidata(hObject, handles);

%---

function exit_Callback(hObject, handles);
if strcmp(get(handles.save, 'enable'), 'on'); % Ask about saving marks if something has been changed, otherwise confirm and exit
    answer = questdlg([{'The matchpoints have been changed.'} {'Do you want to save the matchpoints before leaving Matchmaker ?'}], 'Save marks ?', 'Yes', 'No', 'Yes');
    switch answer
        case 'Yes'
            save_Callback(hObject, handles);
        case 'No'
            answer = questdlg('Do you want to exit without saving ?', 'Exit ?', 'Yes', 'No', 'No');
            switch answer
                case {'No' ''}
                    return;
            end;
        case ''
            return;
    end;
end;
fileno = handles.fileno;
load(['Subroutines' filesep 'matchmaker' filesep 'matchmaker_sett.mat']);
for i = 1:length(fileno)
    sett.xlim(fileno(i),:) = get(handles.bigax(i), 'xlim');
    sett.specs{fileno(i)} = handles.selectedspecs{i};
end;
save(['Subroutines' filesep 'matchmaker' filesep 'matchmaker_sett.mat'], 'sett');
if isfield(handles, 'evaluatefigurehandle')
    delete(handles.evaluatefigurehandle)
end;
delete(handles.fig)
warning('on', 'MATLAB:Axes:NegativeDataInLogAxis');