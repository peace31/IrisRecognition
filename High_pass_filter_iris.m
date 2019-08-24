ic = imread('NewEye.jpg');
Im1 = im2double(ic); % need to convert image to double for histogram to work

 %% Step 3. Convert RGB image or colormap to intensity image
ImR = Im1(:,:,1);    % Red
ImG = Im1(:,:,2);    %Green
ImB = Im1(:,:,3);    %Blue
ImL = rgb2gray(Im1); %Gray grayscale

px=[-1 0 1;-1 0 1;-1 0 1];
icx=filter2(px,ImL);
figure (1), imshow(icx);
py=px';
icy=filter2(py,ImL);
figure (2), imshow(icx);
pedge=sqrt(icx.^2 + icy.^2);
figure (3), imshow(pedge);
fe=im2bw(pedge,0.25);
figure (4),imshow(fe);

