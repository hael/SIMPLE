function WriteImagic(map, basename)
% function WriteImagic(map, basename);
% Write a 2D image or a stack of 2D images contained in 'map'.
% We assume that all images in a stack are identical in size.  The data are
% written as reals.
% 

[n1 n2 nimages]=size(map);

[nx slen]=size(basename);
% strip any 'img' or 'hed' extension from the base name
if (slen>4) && (basename(slen-3)=='.') % has an extension
    if (basename(slen-2:slen)=='hed') || (basename(slen-2:slen)=='img')
        basename=basename(1:slen-4);  % remove the extension.
    end;
end;
hdrname=strcat(basename,'.hed');
hdr=fopen(hdrname,'wb','ieee-le');  % little-endian
if hdr<0
    error(['File could not be opened for writing: ' hdrname]);
end;

% Fill in the 'minimal' header information
t=clock;
idat=zeros(256,1);
% idat(1)=1;  % Image location number
% idat(2)=max(0,nimages-1); % number of images following.
idat(4)=1; % number of headers per images
idat(5)=t(3);  % day
idat(6)=t(2);  % month
idat(7)=t(1);  % year
idat(8)=t(4);  % hour
idat(9)=t(5);  % minute
idat(10)=round(t(6)); % second
% idat(12)=n1*n2;  % should be left at zero.
idat(13)=n1;  % x pixels
idat(14)=n2;  % y pixels
typeString='REAL';

idat(61)=1;  % number of planes
idat(69)=33686018;  % machine type

for i=1:nimages  % write multiple headers if necessary
    idat(1)=i;          % image location number
    idat(2)=nimages-i; % number of images following
    cnt=fwrite(hdr,idat(1:14),'int32');
    cnt=fwrite(hdr,typeString,'char');
    cnt=fwrite(hdr,idat(16:256),'int32');
end;

fclose(hdr);

% Write the data file
imgname=strcat(basename,'.img');
img=fopen(imgname,'wb','ieee-le');  % little-endian
cnt=fwrite(img,map,'float32');
fclose(img);
