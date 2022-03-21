function HdArray=ReadImagicHeader(basename)
% function HdArray=ReadImagicHeader(basename)
% Read an Imagic header file and return its contents as a structure of 2D arrays.  Numeric contents for
% image i are returned as doubles in HdArray.Vals(:,i).  The strings
% HdArray.Type(i,:), Name and Strings (the History array) are
% also returned.  See MakeImagicHeader for more details.
% We don't care if the name has an extension or not.
%

maxheaders=1e6;  % maximum number of images we can read (100 MB of data)

% Ignore the extension, and construct the filename *.hed
[n1 slen]=size(basename);
if (slen>4) && (basename(slen-3)=='.') % has an extension
    if strcmp(basename(slen-2:slen),'hed') || strcmp(basename(slen-2:slen),'img')
        basename=basename(1:slen-4);  % remove the extension.
    end;
end;
hdrname=strcat(basename,'.hed');

q=typecast(int32(1),'uint8');
machineLE=(q(1)==1);  % true for little-endian machine

% Open the header
fileLE=1;
order='ieee-le';
hdr=fopen(hdrname,'r',order);  % try little-endian
% Try reading it
idat=fread(hdr,4,'int32');
if abs(idat(4))>6e4 % idat(4) is headers/image, typically = 1.
    fclose(hdr);  % wrong byte order, start over.
    fileLE=0;
    order='ieee-be';
    hdr=fopen(hdrname,'r',order);  % big-endian
end;

% fileLE

fseek(hdr,0,'bof');  % move back to the beginning

[IH,count]=fread(hdr,[256 maxheaders+1],'int32');
fclose(hdr);

if count >= 256*(maxheaders+1)
    error('File too long, maxheaders exceeded');
end;

[n nim]=size(IH);

HdArray.Vals=zeros(199,nim);
HdArray.Name=char(zeros(nim,80));
HdArray.Strings=char(zeros(nim,228));
HdArray.Type=char(zeros(1,4));
HdArray.ByteOrder=order;

% rd is a an array showing which elements are reals
rd=zeros(199,1);
rd(18:23)=1;
rd(25:29)=1;
rd(50)=1;
rd(65:67)=1;
rd(100:101)=1;
rd(104:106)=1;
rd(108:120)=1;
rd(123:199)=1;

TheReals=find(rd);
id=1-rd;
% Exclude the characters
id(15)=0;
id(30:49)=0;
id(200:256)=0;
TheInts=find(id);

% And the characters are stored in elements 200 to 256.


for i=1:nim
    HdArray.Vals(TheInts,i)=IH(TheInts,i);
    HdArray.Vals(TheReals,i)=typecast(int32(IH(TheReals,i)),'single');
    if machineLE == fileLE  % File byte-order matches machine
        HdArray.Type(1,:)=char(typecast(int32(IH(15)),'uint8'));
        HdArray.Name(i,:)=char(typecast(int32(IH(30:49,i)),'uint8'));
        HdArray.Strings(i,:)=char(typecast(int32(IH(200:256,i)),'uint8'));
    else  % Big-endian machine with LE file, requires byte-order swap.
        HdArray.Type(1,4:-1:1)=char(typecast(int32(IH(15)),'uint8'));
        HdArray.Name(i,80:-1:1)=char(typecast(int32(IH(49:-1:30,i)),'uint8'));
        HdArray.Strings(i,228:-1:1)=char(typecast(int32(IH(256:-1:200,i)),'uint8'));
    end;
end;
HdArray.ByteOrder=order;
