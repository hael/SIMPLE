function [m,s,a]=ReadTiffStack(filename,startSlice,numSlices)
% [m,s,a]=ReadTiffStack(filename,startSlice,numSlices)
% Returns the stack of images m.  The structure s is returned regardless of
% the number of slices requested.  The pixelsize s.pixA is read from the
% imageDescription field of files that we write with WriteTiffStack;
% otherwise it's returned as zero.  The returned array of structs s are the
% image info returned by iminfo().

if nargin<2
    startSlice=1;
end;
if nargin<3
    numSlices=inf;
end;
a=imfinfo(filename);
nz=numel(a);  % number of slices
start=min(nz,startSlice);
num=max(0,min(nz-startSlice+1,numSlices));
s=struct;
% in files we have created, we have Matlab code 
doEval=isfield(a(1),'ImageDescription');
nRot=0;
if doEval
    try
        eval(a(1).ImageDescription);  % assign s.pixA, nx, ny and nz variables.
    catch
        doEval=0;
    end;
end;
if doEval % we created this file, s values are evaluated already.
    flipImages=0;
else
    if strcmp(a(1).ResolutionUnit,'Inch') && (numel(a(1).XResolution)>0) && (a(1).XResolution ~= 0)  % SerialEM fill this in
        s.pixA=25.4e7/a(1).XResolution;
    else
        s.pixA=0; % unknown value, we have a DM TIFF
    end;
%     standard k2 movie is 3838x3710 = (Height,Width) = nx,ny
    flipImages=(a(1).Width>a(1).Height); % Check for non-rotated SerialEM option.
%     k2 data have flipImages=0.
end;
    if flipImages || ~doEval % k2 data
        s.nx=a(1).Height;
        s.ny=a(1).Width;
    else
        nRot=3;
        s.nx=a(1).Width;  % doEval data: s is already already assigned from string.
        s.ny=a(1).Height;
    end;
    s.nz=nz;

% Fill in fields for compatibility with ReadMRC.
s.mx=s.nx;
s.my=s.ny;
s.mz=s.nz;
s.rez=s.pixA*s.mz;

if a(1).BitDepth>8
    m=zeros(s.nx,s.ny,num,'uint16');
else
    m=zeros(s.nx,s.ny,num,'uint8');
end;
for i=1:num
    if flipImages
        m(:,:,i)=flip(imread(filename,'tiff','index',i+start-1,'info',a),1);
    else
        m(:,:,i)=rot90(imread(filename,'tiff','index',i+start-1,'info',a),nRot);
    end;
end;
