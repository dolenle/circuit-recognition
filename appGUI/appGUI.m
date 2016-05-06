function varargout = appGUI(varargin)
%APPGUI M-file for appGUI.fig
%      APPGUI, by itself, creates a new APPGUI or raises the existing
%      singleton*.
%
%      H = APPGUI returns the handle to a new APPGUI or the handle to
%      the existing singleton*.
%
%      APPGUI('Property','Value',...) creates a new APPGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to appGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      APPGUI('CALLBACK') and APPGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in APPGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help appGUI

% Last Modified by GUIDE v2.5 06-May-2016 09:53:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @appGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @appGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before appGUI is made visible.
function appGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for appGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes appGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
set(hObject,'Resize','off');


% --- Outputs from this function are returned to the command line.
function varargout = appGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.source_dir = uigetdir();
set(handles.edit3, 'String', handles.source_dir);
addpath(genpath(handles.source_dir));
guidata(hObject, handles);


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dest_dir = uigetdir();
set(handles.edit4, 'String', handles.dest_dir);
guidata(hObject, handles);


function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
labeled_feature_mtx(handles.source_dir,handles.dest_dir, handles);

% --- Executes on button press in checkbox1 (Select loaded .MAT).
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mat_source_dir = handles.dest_dir;
set(handles.edit5, 'String', handles.mat_source_dir);
addpath(genpath(handles.mat_source_dir));
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox1



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6 (select .MAT)
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mat_source_dir = uigetdir();
set(handles.edit5, 'String', handles.mat_source_dir);
addpath(genpath(handles.mat_source_dir));
guidata(hObject, handles);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path = get(handles.edit5,'String');
if(exist(path))
    dir_struct = dir(fullfile(path,'*mat'));
    if(numel(dir_struct) == 0)
        set(handles.text11, 'String', 'Files not found.');
    else  
        set(handles.text11, 'String', 'Feature/Label Matrix Loaded.');
    end
else
    set(handles.text11, 'String', 'Invalid Directory.');
end


function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8 (select ckt)
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path, filter] = uigetfile('*.jpg;*.JPG;*.jpeg;*.JPEG;*.png;*.bmp', 'Schematic Images');
handles.ckt_source = fullfile(path,file);
set(handles.edit6, 'String', handles.ckt_source);
addpath(genpath(path));
axes(handles.axes6);
imshow(~im2bw(imread(handles.ckt_source)));
guidata(hObject, handles);


function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uiputfile({'*.net'},'Save as');
handles.net_dest = fullfile(path,file);
set(handles.edit7, 'String', handles.net_dest);
guidata(hObject, handles);


% --- Executes on button press in pushbutton11. (RUN)
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
minLength = str2num(get(handles.edit8, 'String'));
run(handles.ckt_source, handles.mat_source_dir, handles.mat_source_dir, handles.net_dest, minLength);
guidata(hObject, handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
axes(hObject);
imshow('blocks.png');
% Hint: place code in OpeningFcn to populate axes5


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton6.
function pushbutton6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function axes5_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on pushbutton7 and none of its controls.
function pushbutton7_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Read In Scanned Image & Crop Border
minLength = str2num(get(handles.edit8, 'String'));
scanned_image = imread(handles.ckt_source);
cropped_scan = imcrop(scanned_image, [10 10 size(scanned_image,2)-20 size(scanned_image,1)-20]); %crop 10px edges
% Resize scan image and invert image
resized_scan = imresize(cropped_scan,[4096 4096]);
resized_scan = ~im2bw(resized_scan,0.6);
axes(handles.axes6);
% imshow(resized_scan); hold on;
% Thin Scanned Image 
dilate_resized_scan = imdilate(resized_scan, strel('disk', 10));
thin_scan = bwmorph(dilate_resized_scan, 'thin', Inf);
    %figure;imshow(thin_scan); title('thin_scan');
% Denoising by removing island size
    % Remove islands smaller than 1/10th of average
    %avgIslandSize = mean(cellfun(@length, CC.PixelIdxList))
removed_specks = 0;
CC = bwconncomp(thin_scan);
denoised_scan = thin_scan;
for aa = 1:CC.NumObjects
    %if(length(CC.PixelIdxList{aa}) < avgIslandSize/10)
    if(length(CC.PixelIdxList{aa}) < 500)
        denoised_scan(CC.PixelIdxList{aa}) = 0;
        removed_specks = removed_specks+1;
    end
end
less_denoised_scan = thin_scan;
for aa = 1:CC.NumObjects
    %if(length(CC.PixelIdxList{aa}) < avgIslandSize/10)
    if(length(CC.PixelIdxList{aa}) < 50)
        less_denoised_scan(CC.PixelIdxList{aa}) = 0;
    end
end
removed_specks;
    %figure;imshow(denoised_scan); title('denoised scan');
    %figure;imshow(less_denoised_scan); title('less denoised scan');

% Thicken Images and Transpose for Hough Transform
thick_denoised_scan = imdilate(denoised_scan, strel('rectangle', [1,20])); % Used 'disk' @ 930PM
thick_denoised_scan_rot = imdilate(denoised_scan', strel('rectangle', [1,20]));
%thick_denoised_scan = thick_denoised_scan_rot;
    %figure; imshow(thick_denoised_scan); title('thick_denoised_scan');

% Compute Hough Transform for Line Detection
[H,theta,rho] = hough(thick_denoised_scan);
[H_rot,theta_rot,rho_rot] = hough(thick_denoised_scan_rot);
% Visualize Hough Transform
    % figure;imshow(imadjust(mat2gray(H)),'XData',theta,'YDATA',rho,'InitialMagnification','fit');
    % axis on;axis normal,hold on; colormap(hot);
peaks = houghpeaks(H,100,'threshold',ceil(0.3*max(H(:))));
peaks_rot = houghpeaks(H_rot,100,'threshold',ceil(0.3*max(H_rot(:))));
    % plot(theta(peaks(:,2)),rho(peaks(:,1)),'s','color','black');

% Find Lines of Hough Transform
lines = houghlines(thick_denoised_scan,theta,rho,peaks,'FillGap',3,'MinLength',minLength);
pre_deleted_lines = lines;
lines_rot = houghlines(thick_denoised_scan_rot,theta_rot,rho_rot,peaks_rot,'FillGap',3,'MinLength',minLength);
pre_deleted_lines_rot = lines_rot;
    %hold off;

%figure; imshow(~imdilate(denoised_scan, strel('disk', 5))); title('not thick denoised scan'); hold on;
for ii = 1:2
    if (ii == 2)
        lines = lines_rot;
    end
    % Remove angled lines (i.e. BJT) and any lines which intersect them=
    angles = [lines.theta];
    slantInd = find(angles < 50 & angles > 40 | angles > -50 & angles < -40);
    slantLines = horzcat(vec2mat([lines(slantInd).point1],2),vec2mat([lines(slantInd).point2],2)); %End points of angled lines
    row_index = [];
    output_intersections = lineSegmentIntersect(horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2)), slantLines);
    [row_index,col_index] = find(output_intersections.intAdjacencyMatrix);
    
    slantInd = unique([slantInd row_index']);
    % Delete Slanted & Intersected Lines
    lines(slantInd) = [];
    linesRemoved = length(slantInd);
    lineVec = horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2));
    %disp 'Joining nearby points'
    allpts = vertcat(lineVec(:,1:2),lineVec(:,3:4)); % everything x two
    ptsJoined = 0;
    thresh = 25;
    for k = 1:size(allpts, 1)
        x = allpts(k,:);
        y = vertcat(allpts(1:k-1,:),allpts(k+1:end,:));
        distances = sqrt((y(:,1)-x(1)).^2 + (y(:,2)-x(2)).^2);
        pts = vertcat(x, y(find(distances < thresh),:));
        if(size(pts,1) > 1)
            pt = round(mean(pts));
            idx = find(distances < thresh);
            idx = horzcat(idx(idx<k)', k, idx(idx>k)');
            idx(idx>k) = idx(idx>k) + 1;
            allpts(idx,:) = repmat(pt,length(idx),1);
            ptsJoined = ptsJoined + length(pts)-1;
        end
    end
    %ptsJoined
    
    % Remove duplicate lines
    %disp 'Removing duplicate lines'
    line_mtx_no_duplicates = unique(horzcat(allpts(1:size(allpts,1)/2,:),allpts(size(allpts,1)/2+1:end,:)),'rows');
    %linesRemoved = size(allpts,1)/2-size(line_mtx_no_duplicates,1)
    if (ii == 1)
        line_mtx = line_mtx_no_duplicates;
    else
        swaped_xy = zeros(size(line_mtx_no_duplicates));
        swaped_xy(:,1) = line_mtx_no_duplicates(:,2);
        swaped_xy(:,2) = line_mtx_no_duplicates(:,1);
        swaped_xy(:,3) = line_mtx_no_duplicates(:,4);
        swaped_xy(:,4) = line_mtx_no_duplicates(:,3);
        line_mtx = [line_mtx; swaped_xy];
    end
end

% Erasing Lines From Original Image
%ckt_no_line = denoised_scan;
ckt_no_line = less_denoised_scan;
%ckt_no_line = imdilate(denoised_scan, strel('disk', 5));
point_mtx = line_mtx;
erase_thickness = 30;
erase_limit = 10;
for ii = 1:size(point_mtx,1)
    if (abs(point_mtx(ii,1)-point_mtx(ii,3))< abs(point_mtx(ii,2)-point_mtx(ii,4)))
        x1 = point_mtx(ii,2); y1 = point_mtx(ii,1);
        x2 = point_mtx(ii,4); y2 = point_mtx(ii,3);
        horiz = 0;
    else
        x1 = point_mtx(ii,1);y1 = point_mtx(ii,2);
        x2 = point_mtx(ii,3);y2 = point_mtx(ii,4);
        horiz = 1;
    end
    mm = (y2-y1)/(x2-x1);
    bb = y1 - mm*x1;

    for xx = (min(x1,x2)+erase_limit):((max(x1,x2))-erase_limit)
        yy = round(mm*xx+bb);
        yy_min = yy-erase_thickness;
        yy_max = yy+erase_thickness;
        if (horiz == 1)
            ckt_no_line([yy_min:yy_max],xx)=0;
        else
            ckt_no_line(xx,[yy_min:yy_max])=0;
        end
    end
end
prev_img = imdilate(ckt_no_line, strel('disk', 10));
% prev_img_bb = regionprops(prev_img, 'BoundingBox');
% prev_img_bb = prev_img_bb.BoundingBox;
xS = find(sum(prev_img) > 0);
yS = find(sum(prev_img, 2) > 0);
imshow(imcrop(prev_img, [yS(1) xS(1) yS(end)-yS(1) xS(end)-xS(1)]));
% imshow(prev_img);
