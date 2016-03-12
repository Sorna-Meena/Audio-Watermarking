% Robust Blind Digital Audio Watermarking                                   
% Implementation Based on Watermark Scrambling Algorithm-Arnold Transform   
% Discrete Wavelet Transform,Discrete Cosine Transform,Erro Correcting Code 
% Coded by Subir, Undergraduate Sudent, ECE, LNMIIT, Jaipur                 
% Under Guidance of Prof. Amit M. Joshi, MNIT, Jaipur                       

clc
clear all
close all
warning off

%% audio input
ok_classical=0;   % 1 selected 0 not-selected
ok_jazz=0;
ok_pop=1;
ok_looney=0;
 if ok_classical
[y,Fs] = audioread('classical.wav');
 end
 if ok_jazz
[y,Fs] = audioread('jazz.wav');
 end
 if ok_pop
[y,Fs] = audioread('pop.wav');
 end
 if ok_looney
[y,Fs] = audioread('loopyMusic.wav');
 end
x1=y(1:262144);
x1=reshape(x1,1,262144);
x=reshape(x1,512,512);
Key=7;
[L1, L2]=size(x1);
%% image input 

img  = imread('CW32.jpg'); %Get the input image
I = rgb2gray(img);
original_img=im2bw(I);
w1=im2bw(I);%Convert to grayscale image

%% 3 level of DWT decomposition with HAAR wavelet

[L,H,V,D]=dwt2(x,'haar','sym');
[CL,CH,CV,CD]=dwt2(H,'haar','sym');
[CA1,CA2,CA3,CA4]=dwt2(CH,'haar','sym');
for i=1:16
    a(i)=4;
end
for i=1:16
    b(i)=4;
end
%% DCT of sub-band block CH->CA2

newCA2=mat2cell(CA2,a,b);
d=zeros(64,64);
nCA2=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nCA2{i,j}=dct2(newCA2{i,j});
    end
end
%% DCT of sub-band block CV->CV2

[CV1,CV2,CV3,CV4]=dwt2(CV,'haar','sym');

newCV2=mat2cell(CV2,a,b);
d=zeros(64,64);
nCV2=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nCV2{i,j}=dct2(newCV2{i,j});
    end
end
%% DCT of sub-band block VH->VH2

[VL,VH,VV,VD]=dwt2(V,'haar','sym');
[VH1,VH2,VH3,VH4]=dwt2(VH,'haar','sym');
newVH2=mat2cell(VH2,a,b);
d=zeros(64,64);
nVH2=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nVH2{i,j}=dct2(newVH2{i,j});
    end
end
%% dct OF SUB-band VV->VV2

[VV1,VV2,VV3,VV4]=dwt2(VV,'haar','sym');
newVV2=mat2cell(VV2,a,b);
d=zeros(64,64);
nVV2=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nVV2{i,j}=dct2(newVV2{i,j});
    end
end
%% UNCORRELATED PSEUDO RANDOM SEQUENCE GENERATOR

G=Key;
%Generation of first m-sequence using generator polynomial [45]
sd1 =[0 0 0 0 1];             % Initial state of Shift register
PN1=[];                       % First m-sequence
for j=1:G
    PN1=[PN1 sd1(5)];
    if sd1(1)==sd1(4)
        temp1=0;
    else temp1=1;
    end
    sd1(1)=sd1(2);
    sd1(2)=sd1(3);
    sd1(3)=sd1(4);
    sd1(4)=sd1(5);
    sd1(5)=temp1;
end
sd2 =[0 0 0 0 1];             % Initial state of Shift register
PN2=[];                       % Second m-sequence
for j=1:G
    PN2=[PN2 sd2(5)];
    if sd2(1)==sd2(2)
        temp1=0;
    else temp1=1;
    end
    if sd2(4)==temp1
        temp2=0;
    else temp2=1;
    end
    if sd2(5)==temp2
        temp3=0;
    else temp3=1;
    end
    sd2(1)=sd2(2);
    sd2(2)=sd2(3);
    sd2(3)=sd2(4);
    sd2(4)=sd2(5);
    sd2(5)=temp3;
end
%% WATERMARKING COEFFICIENT OR WEIGHT ALPHA

alpha=0.1;
%% WATERMARK BIT ENCODING ALGORITHM IN SELECTED/MID-BAND COEFFICIENTS OF DCT TRANSFORMED BLOCKS ONLY
%% FOR BLOCK CA2

for p=1:16
    for q=1:16
        if w1(p,q)==0
            nCA2{p,q}(1,3)=nCA2{p,q}(1,3)+alpha*PN1(1);
            nCA2{p,q}(1,4)=nCA2{p,q}(1,4)+alpha*PN1(2);
            nCA2{p,q}(2,2)=nCA2{p,q}(2,2)+alpha*PN1(3);
            nCA2{p,q}(2,3)=nCA2{p,q}(2,3)+alpha*PN1(4);
            nCA2{p,q}(3,1)=nCA2{p,q}(3,1)+alpha*PN1(5);
            nCA2{p,q}(3,2)=nCA2{p,q}(3,2)+alpha*PN1(6);
            nCA2{p,q}(4,1)=nCA2{p,q}(4,1)+alpha*PN1(7);
        else
            nCA2{p,q}(1,3)=nCA2{p,q}(1,3)+alpha*PN2(1);
            nCA2{p,q}(1,4)=nCA2{p,q}(1,4)+alpha*PN2(2);
            nCA2{p,q}(2,2)=nCA2{p,q}(2,2)+alpha*PN2(3);
            nCA2{p,q}(2,3)=nCA2{p,q}(2,3)+alpha*PN2(4);
            nCA2{p,q}(3,1)=nCA2{p,q}(3,1)+alpha*PN2(5);
            nCA2{p,q}(3,2)=nCA2{p,q}(3,2)+alpha*PN2(6);
            nCA2{p,q}(4,1)=nCA2{p,q}(4,1)+alpha*PN2(7);
        end
    end
end

for i=1:16
    for j=1:16
        newCA2{i,j}=idct2(nCA2{i,j});
    end
end
CA2=cell2mat(newCA2);
%% FOR BLOCK CV2

for p=1:16
    for q=1:16
        if w1(p,16+q)==0
            nCV2{p,q}(1,3)=nCV2{p,q}(1,3)+alpha*PN1(1);
            nCV2{p,q}(1,4)=nCV2{p,q}(1,4)+alpha*PN1(2);
            nCV2{p,q}(2,2)=nCV2{p,q}(2,2)+alpha*PN1(3);
            nCV2{p,q}(2,3)=nCV2{p,q}(2,3)+alpha*PN1(4);
            nCV2{p,q}(3,1)=nCV2{p,q}(3,1)+alpha*PN1(5);
            nCV2{p,q}(3,2)=nCV2{p,q}(3,2)+alpha*PN1(6);
            nCV2{p,q}(4,1)=nCV2{p,q}(4,1)+alpha*PN1(7);
        else
            nCV2{p,q}(1,3)=nCV2{p,q}(1,3)+alpha*PN2(1);
            nCV2{p,q}(1,4)=nCV2{p,q}(1,4)+alpha*PN2(2);
            nCV2{p,q}(2,2)=nCV2{p,q}(2,2)+alpha*PN2(3);
            nCV2{p,q}(2,3)=nCV2{p,q}(2,3)+alpha*PN2(4);
            nCV2{p,q}(3,1)=nCV2{p,q}(3,1)+alpha*PN2(5);
            nCV2{p,q}(3,2)=nCV2{p,q}(3,2)+alpha*PN2(6);
            nCV2{p,q}(4,1)=nCV2{p,q}(4,1)+alpha*PN2(7);
        end
    end
end

for i=1:16
    for j=1:16
        newCV2{i,j}=idct2(nCV2{i,j});
    end
end
CV2=cell2mat(newCV2);
%% FOR BLOCK VH2

for p=1:16
    for q=1:16
        if w1(p+16,q)==0
            nVH2{p,q}(1,3)=nVH2{p,q}(1,3)+alpha*PN1(1);
            nVH2{p,q}(1,4)=nVH2{p,q}(1,4)+alpha*PN1(2);
            nVH2{p,q}(2,2)=nVH2{p,q}(2,2)+alpha*PN1(3);
            nVH2{p,q}(2,3)=nVH2{p,q}(2,3)+alpha*PN1(4);
            nVH2{p,q}(3,1)=nVH2{p,q}(3,1)+alpha*PN1(5);
            nVH2{p,q}(3,2)=nVH2{p,q}(3,2)+alpha*PN1(6);
            nVH2{p,q}(4,1)=nVH2{p,q}(4,1)+alpha*PN1(7);
        else
            nVH2{p,q}(1,3)=nVH2{p,q}(1,3)+alpha*PN2(1);
            nVH2{p,q}(1,4)=nVH2{p,q}(1,4)+alpha*PN2(2);
            nVH2{p,q}(2,2)=nVH2{p,q}(2,2)+alpha*PN2(3);
            nVH2{p,q}(2,3)=nVH2{p,q}(2,3)+alpha*PN2(4);
            nVH2{p,q}(3,1)=nVH2{p,q}(3,1)+alpha*PN2(5);
            nVH2{p,q}(3,2)=nVH2{p,q}(3,2)+alpha*PN2(6);
            nVH2{p,q}(4,1)=nVH2{p,q}(4,1)+alpha*PN2(7);
        end
    end
end

for i=1:16
    for j=1:16
        newVH2{i,j}=idct2(nVH2{i,j});
    end
end
VH2=cell2mat(newVH2);
%% FOR BLOCK VV2

for p=1:16
    for q=1:16
        if w1(p+16,q+16)==0
            nVV2{p,q}(1,3)=nVV2{p,q}(1,3)+alpha*PN1(1);
            nVV2{p,q}(1,4)=nVV2{p,q}(1,4)+alpha*PN1(2);
            nVV2{p,q}(2,2)=nVV2{p,q}(2,2)+alpha*PN1(3);
            nVV2{p,q}(2,3)=nVV2{p,q}(2,3)+alpha*PN1(4);
            nVV2{p,q}(3,1)=nVV2{p,q}(3,1)+alpha*PN1(5);
            nVV2{p,q}(3,2)=nVV2{p,q}(3,2)+alpha*PN1(6);
            nVV2{p,q}(4,1)=nVV2{p,q}(4,1)+alpha*PN1(7);
        else
            nVH2{p,q}(1,3)=nVV2{p,q}(1,3)+alpha*PN2(1);
            nVV2{p,q}(1,4)=nVV2{p,q}(1,4)+alpha*PN2(2);
            nVV2{p,q}(2,2)=nVV2{p,q}(2,2)+alpha*PN2(3);
            nVV2{p,q}(2,3)=nVV2{p,q}(2,3)+alpha*PN2(4);
            nVV2{p,q}(3,1)=nVV2{p,q}(3,1)+alpha*PN2(5);
            nVV2{p,q}(3,2)=nVV2{p,q}(3,2)+alpha*PN2(6);
            nVV2{p,q}(4,1)=nVV2{p,q}(4,1)+alpha*PN2(7);
        end
    end
end

for i=1:16
    for j=1:16
        newVV2{i,j}=idct2(nVV2{i,j});
    end
end
VV2=cell2mat(newVV2);
%% IDWT AND IDCT TO GET MODIFIED COEFFICIENTS AND FORM THE WATERMARKED AUDIO 

CH=idwt2(CA1,CA2,CA3,CA4,'haar','sym');
CV=idwt2(CV1,CV2,CV3,CV4,'haar','sym');
H=idwt2(CL,CH,CV,CD,'haar','sym');
VH=idwt2(VH1,VH2,VH3,VH4,'haar','sym');
VV=idwt2(VV1,VV2,VV3,VV4,'haar','sym');
V=idwt2(VL,VH,VV,VD,'haar','sym');
newx=idwt2(L,H,V,D,'haar','sym');
y1=reshape(newx,1,512*512);
SNR=snr(y1,x1)
y2=reshape(newx,1,512*512);
L=length(y1);
%% DISTURBANCE ADDED TO THE WATERMARKED AUDIO

ok_noise =0;        % select 1 for adding the particular disturbance to host audio
ok_filtering=0;
ok_cropping=0;
ok_resampling=0;
ok_requantization=0;

if ok_noise
    % Additional noise
    y1 = awgn(y1,10,'measured');
end

if ok_filtering
    % Filtering
    myfilter = ones(512,1);
    myfilter = myfilter/sum(myfilter);
    y1   = filter(myfilter,1,y1);
end

if ok_cropping
    % Cropping
    Lmin = round(L/11);
    Lmax = round(5*L/11);
    y1(1:1,1:Lmin)   = 0;
    y1(1:1,Lmax:end) = 0;
end
if ok_resampling
    % Resampling
    Fs_0 = Fs;
    Fs_1 = round(Fs/9);
    y1 = resample(y1,Fs_1,Fs_0);
    y1 = resample(y1,Fs_0,Fs_1);
    if length(y1)<L
        y1(end+L) = 0;
    end
    if length(y1)>L
        y1 = y1(1:L);
    end
end

if ok_requantization
    % Requantization
    bits_new = 8;
    wavwrite(y1,Fs,bits_new,'requantized_sound.wav');
    y1 = audioread('requantized_sound.wav');
end
%% WATERMARK EXTRACTION ALGORITHM

newx=reshape(y1,512,512);
[A,B,C,D]=dwt2(newx,'haar','sym');
[B1,B2,B3,B4]=dwt2(B,'haar','sym');
[B21,B22,B23,B24]=dwt2(B2,'haar','sym');
[B31,B32,B33,B34]=dwt2(B3,'haar','sym');
[C1,C2,C3,C4]=dwt2(C,'haar','sym');
[C21,C22,C23,C24]=dwt2(C2,'haar','sym');
[C31,C32,C33,C34]=dwt2(C3,'haar','sym');
%% DCT OF SUB-BANDS

for i=1:16
    a(i)=4;
end
for i=1:16
    b(i)=4;
end
newB22=mat2cell(B22,a,b);
d=zeros(64,64);
nB22=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nB22{i,j}=dct2(newB22{i,j});
    end
end
%% EXTRACTION OF WATERMARK BIT BY COMPARISON OF CORRELATION

for p=1:16
    for q=1:16
        if corr([nB22{p,q}(1,3) nB22{p,q}(1,4) nB22{p,q}(2,2) nB22{p,q}(2,3) nB22{p,q}(3,1) nB22{p,q}(3,2) nB22{p,q}(4,1)]',PN1(1:7)')>=corr([nB22{p,q}(1,3) nB22{p,q}(1,4) nB22{p,q}(2,2) nB22{p,q}(2,3) nB22{p,q}(3,1) nB22{p,q}(3,2) nB22{p,q}(4,1)]',PN2(1:7)')
            w2(p,q)=0;
        else
            w2(p,q)=1;
        end
    end
end
%% DCT OF SUB-BANDS

for i=1:16
    a(i)=4;
end
for i=1:16
    b(i)=4;
end
newB32=mat2cell(B32,a,b);
d=zeros(64,64);
nB32=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nB32{i,j}=dct2(newB32{i,j});
    end
end
%% EXTRACTION OF WATERMARK BIT BY COMPARISON OF CORRELATION

for p=1:16
    for q=1:16
        if corr([nB32{p,q}(1,3) nB32{p,q}(1,4) nB32{p,q}(2,2) nB32{p,q}(2,3) nB32{p,q}(3,1) nB32{p,q}(3,2) nB32{p,q}(4,1)]',PN1(1:7)')>=corr([nB32{p,q}(1,3) nB32{p,q}(1,4) nB32{p,q}(2,2) nB32{p,q}(2,3) nB32{p,q}(3,1) nB32{p,q}(3,2) nB32{p,q}(4,1)]',PN2(1:7)')
            w2(p,q+16)=0;
        else
            w2(p,q+16)=1;
        end
    end
end

%% DCT OF SUB-BAND

for i=1:16
    a(i)=4;
end
for i=1:16
    b(i)=4;
end
newC22=mat2cell(C22,a,b);
d=zeros(64,64);
nC22=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nC22{i,j}=dct2(newC22{i,j});
    end
end
%% EXTRACTION OF WATERMARK BIT BY COMPARISON OF CORRELATION

for p=1:16
    for q=1:16
        if corr([nC22{p,q}(1,3) nC22{p,q}(1,4) nC22{p,q}(2,2) nC22{p,q}(2,3) nC22{p,q}(3,1) nC22{p,q}(3,2) nC22{p,q}(4,1)]',PN1(1:7)')>=corr([nC22{p,q}(1,3) nC22{p,q}(1,4) nC22{p,q}(2,2) nC22{p,q}(2,3) nC22{p,q}(3,1) nC22{p,q}(3,2) nC22{p,q}(4,1)]',PN2(1:7)')
            w2(p+16,q)=0;
        else
            w2(p+16,q)=1;
        end
    end
end
%% DCT OF SUB-BAND

for i=1:16
    a(i)=4;
end
for i=1:16
    b(i)=4;
end
newC32=mat2cell(C32,a,b);
d=zeros(64,64);
nC32=mat2cell(d,a,b);
for i=1:16
    for j=1:16
        nC32{i,j}=dct2(newC32{i,j});
    end
end
%% EXTRACTION OF WATERMARK BIT BY COMPARISON OF CORRELATION

for p=1:16
    for q=1:16
        if corr([nC32{p,q}(1,3) nC32{p,q}(1,4) nC32{p,q}(2,2) nC32{p,q}(2,3) nC32{p,q}(3,1) nC32{p,q}(3,2) nC32{p,q}(4,1)]',PN1(1:7)')>=corr([nC32{p,q}(1,3) nC32{p,q}(1,4) nC32{p,q}(2,2) nC32{p,q}(2,3) nC32{p,q}(3,1) nC32{p,q}(3,2) nC32{p,q}(4,1)]',PN2(1:7)')
            w2(p+16,q+16)=0;
        else
            w2(p+16,q+16)=1;
        end
    end
end
sbr=w2;
figure
imshow(logical(w2))
%% Quality check
q=max(imabsdiff(y2,reshape(x,1,512*512)))
BER=biterr(original_img,sbr)
ssimval = ssim(uint8(sbr),uint8(original_img))
ps=psnr(uint8(original_img),uint8(sbr)) % Must be greater than 35 dB
ssimValues = zeros(1,10);
qualityFactor = 10:10:100;
for i = 1:10
    imwrite(I,'sbr.jpg','jpg','quality',qualityFactor(i));
    ssimValues(i) = ssim(imread('sbr.jpg'),I);
end
figure
plot(qualityFactor,ssimValues,'r--*')
xlabel(' Copression Quality Factor');
ylabel(' SSIM value');
%% END %%































