function varargout = Cholesterol(varargin)
% CHOLESTEROL M-file for Cholesterol.fig
%      CHOLESTEROL, by itself, creates a new CHOLESTEROL or raises the existing
%      singleton*.
%
%      H = CHOLESTEROL returns the handle to a new CHOLESTEROL or the handle to
%      the existing singleton*.
%
%      CHOLESTEROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOLESTEROL.M with the given input arguments.
%
%      CHOLESTEROL('Property','Value',...) creates a new CHOLESTEROL or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before Cholesterol_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Cholesterol_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Cholesterol

% Last Modified by GUIDE v2.5 03-Mar-2016 14:45:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cholesterol_OpeningFcn, ...
                   'gui_OutputFcn',  @Cholesterol_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before Cholesterol is made visible.
function Cholesterol_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Cholesterol (see VARARGIN)

% Choose default command line output for Cholesterol
handles.output = hObject;
handles.flag=0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Cholesterol wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Cholesterol_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_fileOpen.
function btn_fileOpen_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fileOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
    {'*.bmp;*.tif;*.tiff;*.jpg;*.jpeg','All Picture Files (*.bmp;*.tif;*.tiff;*.jpg;*.jpeg)'; ...
    '*.bmp',  'Bitmap File (*.bmp)'; ...
    '*.tif;*tiff','TIFF (*.tif;*.tiff)'; ...
    '*.jpg;*.jpeg','JPEG (*.jpg;*.jpeg)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Insert Pre-process Latex Glove Sample Image');
if isequal(filename,0) || isequal(pathname,0)
    return
end

Im = imread([pathname,filename]);

handles.Im = Im;
handles.flag = 1;
guidata(hObject,handles)

axes(handles.Pic1);
imshow(handles.Im);
title('Localization Iris and Pupil');


function txt_out_Callback(hObject, eventdata, handles)
% hObject    handle to txt_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_out as text
%        str2double(get(hObject,'String')) returns contents of txt_out as a double


% --- Executes during object creation, after setting all properties.
function txt_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_close.
function btn_close_Callback(hObject, eventdata, handles)
% hObject    handle to btn_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close;

% --- Executes on button press in btn_calculation.
function btn_calculation_Callback(hObject, eventdata, handles)
% hObject    handle to btn_calculation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( handles.flag==1)
        Im=handles.Im;
        Im1 = im2double(Im); % need to convert image to double for histogram to work
        if(size(Im1,1)>=1000)
            Im1=imresize(Im1,200/size(Im1,1));
        end    
        
        ImR = Im1(:,:,1);    % Red
        ImG = Im1(:,:,2);    %Green
        ImB = Im1(:,:,3);    %Blue
        ImL = rgb2gray(Im1); %Gray grayscale


        [counts,x] = imhist(ImB); % Histogram of image data
            thresh = 200/255;
            for x = 1: size(ImL,1)
                for y = 1: size(ImB,2)
                    if ImG(x,y) < thresh
                        h(x,y) = 0;
                    else
                        h(x,y) = ImB(x,y);
                    end;
                end
            end

         BW = edge(h,'canny',0.6); %specifies sensitivity thresholds for the Canny method. 0.7
         radii = 45:1:50; 
         H = circle_hough(BW, radii, 'same', 'normalise');

        peaks = circle_houghpeaks(H, radii, 'nhoodxy', 15, 'nhoodr', 21, 'npeaks', 1);
    %%%%%    
        guidata(hObject,handles)
        axes(handles.Pic2);
        imshow(ImL);
        title('Iris and Pupil area');
        hold on;
        for peak = peaks
            [x, y] = circlepoints(peak(3));
            plot(x+peak(1), y+peak(2), 'g-','LineWidth',2);
            circlepoint = [x+peak(1); y+peak(2)];
        end


        xcentre = nanmedian(circlepoint(1,:));
        ycentre = nanmedian(circlepoint(2,:));

        radius = pdist([xcentre, circlepoint(1,1); ycentre, circlepoint(2,1)],'euclidean')*0.95;
  %%
        mCount=0;mthred=0;
        for i=xcentre-7:xcentre+7
            for j=ycentre-7:ycentre+7
                if(ImL(j,i)<=0.15 )
                    mCount=mCount+1;
                    mthred=mthred+ImL(j,i); 
                end
            end
        end
        mthred=mthred/mCount;

        for i=ycentre+radius:-1:ycentre
            if(ImL(int16(i),int16(xcentre))<=mthred*1.1)
                radiusPupil=i-ycentre;
                break;
            end
        end

        num=0;
        for t=0:0.1:8*atan(1)+0.1
            num=num+1;
            xx(num)=xcentre+radiusPupil*sin(t);
            yy(num)=ycentre+radiusPupil*cos(t);
        end
        
        plot(xx,yy,'r-','LineWidth',2);
        hold off;
            
        NormalizedImag=zeros(100,170);
        numx=0;
        for t=0:8*atan(1)/170:8*atan(1)
            numx=numx+1;numy=0;
            for r=0:radius/100:radius        
                numy=numy+1;
                x1=xcentre+r*sin(t);
                y1=ycentre+r*cos(t);
                NormalizedImag(numy,numx)=ImL(int16(y1),int16(x1));
            end
        end
        guidata(hObject,handles)
        axes(handles.Pic3);
        imshow(NormalizedImag);
        title('100% Normalization Polar & rec');
        
        NormalizedImag30=NormalizedImag(int16(numy*0.7):numy,:);        
        guidata(hObject,handles)
        axes(handles.Pic4);
        imshow(NormalizedImag30);     
        title('30% Normalization Polar & rec');
%%
        Nmax=1000;
        R=5;
        for n=1:Nmax

            r2(n) = R*sqrt(rand(1,1));
            % and theta as before:
            theta2(n)=2*pi*rand(1,1);   
            % convert to cartesian
            x2(n)=r2(n)*cos(theta2(n));
            y2(n)=r2(n)*sin(theta2(n));
        end


        for y = 1: size(ImL,1)
            for x = 1: size(ImB,2)
                pnt = sqrt(power((x - xcentre),2) +  power((y - ycentre),2));
                if pnt <= radius                    
                    ImIris(y,x) = ImL(y,x);
                else 
                    ImIris(y,x) = 0;
                end;
            end
        end


        [counts,x2] = imhist(ImIris); % Histogram of image data
        counts(1) = 0;

        px=[1/9 1/9 1/9;1/9 1/9 1/9;1/9 1/9 1/9]; % high pass filter
        px=[-1 0 1;-1 0 1;-1 0 1]; % Low pass filter
        icx=filter2(px,ImIris);
        %figure (11), imshow(icx);
        py=px';
        icy=filter2(py,ImIris);
        %figure (12), imshow(icx);
        pedge=sqrt(icx.^2 + icy.^2);

        %%
        %handles.pedge = pedge;
        %guidata(hObject,handles)
        %axes(handles.Pic2);
        %imshow(handles.pedge); 

        %%
        Hitogram=imhist(NormalizedImag30);
%         Hitogram=imhist(pedge);
%         Hitogram(1)=0;Hitogram(256)=0;
%         for y = 1: size(ImL,1)
%             for x = 1: size(ImB,2)
%                 pnt = sqrt(power((x - xcentre),2) +  power((y - ycentre),2));
%                 if pnt <= radius
%                     if(pedge(y,x)==0)
%                         Hitogram(1)=Hitogram(1)+1;
%                     else if(pedge(y,x)==1)
%                         Hitogram(256)=Hitogram(256)+1;
%                         end
%                     end
%                 end;
%             end
%         end

        guidata(hObject,handles)
        axes(handles.Pic5);
        imhist(NormalizedImag30);
        title('Histogram');
%         for i=1:256
%             line([i,i],[0,Hitogram(i)],'Color','b');
%         end

        %%
        for t=1:255
            mu1(t)=0.0; mu2(t)=0.0;
            sigma1(t)=0.0; sigma2(t)=0.0;
        end

        totalC=sum(Hitogram);
        for t=1:255
            q1(t)=sum(Hitogram(1:t))/totalC;
            q2(t)=sum(Hitogram(t+1:256))/totalC;
            for i=1:255
                if(i<=t)
                    if(q1(t)~=0) 
                        mu1(t)=mu1(t)+i*double(Hitogram(i)/totalC)/q1(t);
                    end            
                else
                    if(q2(t)~=0) 
                        mu2(t)=mu2(t)+i*double(Hitogram(i)/totalC)/q2(t);
                    end
                end
            end
            for i=1:256
                if(i<=t)
                    if(q1(t)~=0) 
                        sigma1(t)=sigma1(t)+(i-mu1(t))^2*double(Hitogram(i)/totalC)/q1(t);
                    end
                else
                     if(q2(t)~=0) 
                        sigma2(t)=sigma2(t)+(i-mu2(t))^2*double(Hitogram(i)/totalC)/q2(t);
                     end
                end
            end    
            sigmaW(t)=q1(t)*sigma1(t)+q2(t)*sigma2(t);  
            sigmaT(t)=q1(t)*(1-q1(t))*(mu1(t)-mu2(t))^2;
        end
        clear totalC Hitogram sigma1 sigma2 q1 q2 mu1 mu2;

        %%
        [mmax mposition]=max(sigmaT);
        guidata(hObject,handles)
        axes(handles.Pic6);        
        plot(sigmaT,'-b','LineWidth',2);
        title('OTSU thredhold');
        line([mposition,mposition],[0,mmax],'Color','r','LineWidth',2,'Marker','.','LineStyle','-');
        set(handles.txt_out,'String',mposition);
end
 
