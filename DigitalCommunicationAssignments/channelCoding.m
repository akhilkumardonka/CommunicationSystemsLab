clc;
clear all;
close all;

%% Image Input and Framing

image = imread('./car.jpg');
imshow(image);
bw = rgb2gray(image);
[height,width]=size(bw);

% Converting to bitstream and k=4 bits sized grouping
image_1D=reshape(bw.',1,[]);
inputMatrix=de2bi(image_1D,'left-msb');
im_uncoded_1D = reshape(inputMatrix.',1,[]);
Image_4bit=reshape(im_uncoded_1D,4,[]).'; 

%% (7,4) LBC

P=[1 1 1;1 0 1;0 1 1;1 1 0];
G=[eye(4),P];
H=[P',eye(3)];

%% Adding Errors/Noise to Image Bits and Uncoded Txn

p = [0.001, 0.01, 0.1, 0, 0.5, 0.9, 0.99, 0.999];
error_matrix_uncoded = uint8(zeros(length(p),8*height*width));
noisy_im_uncoded = uint8(zeros(length(p),8*height*width));
recon_image1 = uint8(zeros(height*width,8,length(p)));
recon_image2 = uint8(zeros(height*width,1,length(p)));
recon_image_uncoded_final = uint8(zeros(height,width,length(p)));

for i=1:length(p) 
    error_matrix_uncoded(i,:) = uint8(randsrc(1,8*height*width,[1,0; p(i),(1-p(i))]));
    noisy_im_uncoded(i,:) = uint8(mod(im_uncoded_1D+error_matrix_uncoded(i,:),2));
    recon_image1(:,:,i) = reshape(noisy_im_uncoded(i,:),8,[]).';
    recon_image2(:,:,i) = uint8(bi2de(recon_image1(:,:,i),'left-msb'));          
    recon_image_uncoded_final(:,:,i) = reshape(recon_image2(:,:,i).',width,[]).';
end

%% With Channel Coding

imCoded=uint8(mod(double(Image_4bit)*G,2));
imCoded_1D=reshape(imCoded.',1,[]);                           

error_matrix = uint8(zeros(length(p),14*height*width));
noisy_coded = uint8(zeros(length(p),14*height*width));

for i=1:length(p)               
    error_matrix(i,:) = uint8(randsrc(1,14*height*width,[1,0;p(i),(1-p(i))]));
    noisy_coded(i,:) = mod(double(error_matrix(i,:)+imCoded_1D),2);
end

%% Reconstruction after Coded Txn

imRxed = uint8(zeros(2*height*width,7,length(p)));
Syndrome = uint8(zeros(2*height*width,3,length(p)));
Index = uint8(zeros(2*height*width,1,length(p)));
corrected_coded_image = uint8(zeros(2*height*width,7,length(p)));
decoded_image = uint8(zeros(2*height*width,4,length(p)));
decoded_image1 = uint8(zeros(1,8*height*width,length(p)));
decoded_image_reshaped = uint8(zeros(height*width,8,length(p)));
decoded_image_gray_1D = uint8(zeros(height*width,1,length(p)));
decoded_image_final = uint8(zeros(height,width,length(p)));

% since it is systematic form, selecting first 4 bits
decoding_matrix=[1,0,0,0,0,0,0;0,1,0,0,0,0,0;0,0,1,0,0,0,0;0,0,0,1,0,0,0].';

for i=1:length(p)
    imRxed(:,:,i)=reshape(noisy_coded(i,:),7,[]).';
    corrected_coded_image(:,:,i)=imRxed(:,:,i);
    Syndrome(:,:,i)=uint8(mod(double(imRxed(:,:,i))*H',2));
    Index(:,:,i)=uint8(bi2de(Syndrome(:,:,i),'right-msb'));        
    for j=1:2*height*width
        if Index(j,1,i)~=0             
           corrected_coded_image(j,Index(j,1,i),i)=uint8(mod(corrected_coded_image(j,Index(j,1,i),i)+1,2));
        end
    end
    decoded_image(:,:,i)=uint8(mod(double(corrected_coded_image(:,:,i))*decoding_matrix,2));
    decoded_image1(:,:,i) = reshape(decoded_image(:,:,i).',1,[]);
    decoded_image_reshaped(:,:,i)=reshape(decoded_image1(:,:,i),8,[]).';
    decoded_image_gray_1D(:,:,i)=uint8(bi2de(decoded_image_reshaped(:,:,i),'left-msb'));
    decoded_image_final(:,:,i)=reshape(decoded_image_gray_1D(:,:,i).',width,[]).';
end

%% Plots

f1=figure;
subplot(2,4,1);
imshow(image),title('(Original colour image)');
subplot(2,4,5);
imshow(bw),title('(Original Grayscale image)');
subplot(2,4,2);
imshow(uint8(recon_image_uncoded_final(:,:,1))),title('Reconstucted Uncoded image (Pe=0.001)');
subplot(2,4,3);
imshow(uint8(recon_image_uncoded_final(:,:,2))),title('Reconstucted Uncoded image (Pe=0.01)');
subplot(2,4,4);
imshow(uint8(recon_image_uncoded_final(:,:,3))),title('Reconstucted Uncoded image (Pe=0.1)');
subplot(2,4,6);
imshow(uint8(decoded_image_final(:,:,1))),title('Reconstucted Coded image (Pe=0.001)');
subplot(2,4,7);
imshow(uint8(decoded_image_final(:,:,2))),title('Reconstucted Coded image (Pe=0.01)');
subplot(2,4,8);
imshow(uint8(decoded_image_final(:,:,3))),title('Reconstucted Coded image (Pe=0.1)');

f2=figure(2);
subplot(2,3,1),imshow(image),title('(Original colour image)');
subplot(2,3,4),imshow(bw),title('(Original grayscale image)');
subplot(2,3,2),imshow(uint8(recon_image_uncoded_final(:,:,4))),title('Reconstucted Uncoded image (Pe=0)');
subplot(2,3,3),imshow(uint8(recon_image_uncoded_final(:,:,5))),title('Reconstucted Uncoded image (Pb=0.5)');
subplot(2,3,5),imshow(uint8(decoded_image_final(:,:,4))),title('Reconstucted coded image (Pb=0)');
subplot(2,3,6),imshow(uint8(decoded_image_final(:,:,5))),title('Reconstucted coded image (Pb=0.5)');

f3=figure(3);
subplot(2,4,1);imshow(image),title('(Original colour image)');
subplot(2,4,5);imshow(bw),title('(Original grayscale image)');
subplot(2,4,2);imshow(uint8(recon_image_uncoded_final(:,:,6))),title('Reconstucted Uncoded image (Pe=0.9)');
subplot(2,4,3);imshow(uint8(recon_image_uncoded_final(:,:,7))),title('Reconstucted Uncoded image (Pe=0.99)');
subplot(2,4,4);imshow(uint8(recon_image_uncoded_final(:,:,8))),title('Reconstucted Uncoded image (Pe=0.999)');
subplot(2,4,6);imshow(uint8(decoded_image_final(:,:,6))),title('Reconstucted Coded image (Pe=0.9)');
subplot(2,4,7);imshow(uint8(decoded_image_final(:,:,7))),title('Reconstucted Coded image (Pe=0.99)');
subplot(2,4,8);imshow(uint8(decoded_image_final(:,:,8))),title('Reconstucted Coded image (Pe=0.999)');
