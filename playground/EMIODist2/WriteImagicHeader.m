function ok=WriteImagicHeader(HdrArray, name)
% function ok=WriteImagicHeader(HdrArray, name)
%
% Write out an Imagic header file.
% We don't care if the name has an extension or not.
% See MakeImagicHeader() for information on the HdrArray structure.

[nn nim]=size(HdrArray.Vals);

IH=int32(zeros(256,nim));  % image of data to be written

% Ignore the extension, and construct the filename *.hed
[n1 slen]=size(name);
if (slen>4) && (name(slen-3)=='.') % has an extension
    if strcmp(name(slen-2:slen),'hed') || strcmp(name(slen-2:slen),'img')
        name=name(1:slen-4);  % remove the extension.
    end;
end;
hdrname=strcat(name,'.hed');

% Check to see if we're little-endian
q=typecast(int32(1),'uint8');
machineLE=(q(1)==1);  % true for little-endian machine

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
    IH(TheInts,i)=int32(HdrArray.Vals(TheInts,i));
    IH(TheReals,i)=typecast(single(HdrArray.Vals(TheReals,i)),'int32');
    if machineLE
        IH(15,i)=typecast(uint8(HdrArray.Type(1,:)),'int32');
    else
        IH(15,i)=typecast(uint8(HdrArray.Type(1,4:-1:1)),'int32');
    end;
    IH(30:49,i)=typecast(uint8(HdrArray.Name(i,:)),'int32');
    IH(200:256,i)=typecast(uint8(HdrArray.Strings(i,:)),'int32');
end;

% Open the header file and write it
count=0;

hdr=fopen(hdrname,'w',HdrArray.ByteOrder); % We use the same byte order as was read.
count=fwrite(hdr,IH,'int32');
fclose(hdr);

ok=(count==numel(IH));
