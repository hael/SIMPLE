function [map, info]=ReadImagic(basename, startimage, numimages, debug)
% function [map info]=ReadImagic(basename, startimage, numimages, debug)
% Read a stack of 2D images from an Imagic file pair.
% We don't care if the name has an extension or not.
% All arguments except the first are optional.
% The first image is number 1; thus startimage by default is 1.
% If debug>0 extra information will be printed out.
%
% We assume all images are identical in size and data type.
% No complex data types are supported at present.
% The returned array 'map' is the same type (uint8, int16, float) as the
% original data.
%
% The info structure is returned even if numimages is set to zero.
% It contains the fields
% nx, ny, nim = nz (number of pixels in 1st and 2nd dimensions, no. of images)
% string (type string for fread function)
% 
% If numimages<0, the data file is opened but not read.  The file handle is
% returned as 'map'.  Thus to read individual images you can do the
% following:
%   [h s]=ReadImagic(name,1,-1);
%   for i=1:s.nz
%       m1=fread(h,s.nx*s.ny,s.string);
%       m1=reshape(m1,s.nx,s.ny);
%        --do something with m1--
%   end;
%   fclose(h);

% Modified to handle both big and little-endian data, and
% to pick up the EMAN CTF information into the vector CPars.
% fs 28 March 05
% Added debug option fs 1 Sept 05
% Added ability to read part of a file fs 15 Mar 08
% Changed to avoid idat(12) fs 19 Nov 09.
%
% this part is no longer supported:
% the field info.CPars is returned if an EMAN CTF string is found in the
% namestr field of the header.
% The CPars information is
% CPars(1) defocus in um (typically a negative value) e.g. -1.5
% CPars(2) B factor e.g. 200
% CPars(3) Structure amplitude scaler e.g. 1.6
% CPars(4) Amplitude contrast e.g. 0.1
% CPars(5) Noise sqrt exponent
% CPars(6) Noise linear exponent
% CPars(7) Noise squared exponent
% CPars(8) keV
% CPars(9) Cs, in mm
% CPars(10) Angstroms/pixel

if nargin<2 startimage=1; end;
if nargin<3 numimages=inf; end;
if nargin<4 debug=0; end;

% Strip the extension if present in basename.
[path basename ext]=fileparts(basename);
basename=fullfile(path,basename);

% [n1 slen]=size(basename);
% if (slen>4) && (basename(slen-3)=='.') % has an extension
%     if (basename(slen-2:slen)=='hed') | (basename(slen-2:slen)=='img')
%         basename=basename(1:slen-4);  % remove the extension.
%     end;
% end;
hdrname=[basename '.hed'];
t=[];

% Open the header
order='ieee-le';
hdr=fopen(hdrname,'r',order);  % try little-endian
if hdr <= 0
    error(['File could not be opened: ' hdrname]);
end;

% Try reading it
idat=fread(hdr,14,'int32');
if abs(idat(4))>6e4 % idat(4) is headers/image, typically = 1.
    fclose(hdr);  % wrong byte order, start over.
    order='ieee-be';
    hdr=fopen(hdrname,'r',order);  % big-endian
    idat=fread(hdr,14,'int32');
end;

typeString=char(fread(hdr,4,'uchar')');
idat(16:29)=fread(hdr,14,'int32');
namestr=char(fread(hdr,80,'uchar'));

% namestr'

% Pick up the EMAN style CTF parameters
if (namestr(1)=='!') && (namestr(2)=='-') % Contains ctf parameters
    info.CPars=sscanf(namestr(3:78),'%f');
end;

totalImages=idat(2)+1;
info.nim=totalImages;
info.nz=totalImages;  % compatibility with ReadMRC
info.nx=idat(14);
info.ny=idat(13);

% Here we would read the rest of the headers, but we don't...
% t.idat=idat;
% for i=2:nimages
%     [t(i).idat,cnt]=fread(hdr,256,'int32');
% end;
fclose(hdr);

if debug
    disp('Debug output: ');
    totalImages
    typeString
    LengthTypeString=size(typeString);
    LengthTypeString
    namestr'
    idat
end;

map=[];
imgname=[basename '.img'];
img=fopen(imgname,'r',order);

% We assume that all the images are identical in type and size.
% We leave the data in the original numeric type.
% typeString
switch typeString
    case 'PACK'
        string='uint8=>uint8';
        bytes=1;
    case 'INTG'
        string='int16=>int16';
        bytes=2;
    case 'REAL'
        string='float32=>float32';
        bytes=4;
    otherwise
        disp(strcat('ReadImagic: unsupported data mode: ',typeString));
        string='???';
end;
info.string=string;

if numimages<0
    map=img;   % just return the file handle
    return
    
elseif numimages>0  % Read the image file.
    
    %!!!! added line.
    idat(12)=idat(13)*idat(14); % Doesn't seem always to contain the number of words/image.
    
    imagebytes=bytes*double(idat(12));
    
    numimages=double(max(0,min(totalImages-startimage+1,numimages))); % Make sure we don't read beyond end of file.
    ok=1;
    
    if startimage>1
        ok=~fseek(img,double(startimage-1)*imagebytes,'bof');
    end;
    if ok
        [map,cnt]=fread(img,double(numimages*double(idat(12))),string);
        if cnt < numimages*idat(12)
            error(['End of file reached after reading ' num2str(cnt/idat(12)) ' images.']);
        end;
        map=reshape(map,idat(14),idat(13),numimages);
    end;
end;
fclose(img);
