function [m,s,ok]=ReadMovie(name,start,numSlices)
% function [m,s,ok]=ReadMovie(name,start,numSlices)
%   simplest:  m=ReadMovie;  % put up file selector.
% File reader for k2 or f2 camera movies, .mrc, .tif, z.tif, .stk
% Based on the file extension, read the image and return the struct s
% which has fields nx,ny,nz and pixA.  When no pixA value
% is available, as in generic TIFF files, s.pixA is returned as zero.
% previously: handled gzipped files with .gz following the true extension.

m=[]; % default values
s=struct;
ok=true;
showInfo=false;

if nargin<1 || numel(name)<1
    fprintf('Opening the movie file: ');
    typeStr='*.mrc;*.mrcs;*.st;*.tif;*.stk';
    [name,path]=uigetfile(typeStr,'Open a file');
    if isnumeric(name)
        disp(' no file selected.');
        return
    else
        disp([path name]);
        showInfo=true;
    end;
    cd(path);
end;
if nargin<2
    startSlice=1;
end;
if nargin<3
    numSlices=inf;
end;
someData=numSlices>0;

if nargin<2
    start=1;
end;
if nargin<3
    numSlices=inf;
end;

name=InterpretFilename(name);
[~,~,ex]=fileparts(name);
% % Special case for gzipped files:
% isGzFile=strcmp(ex,'.gz');
% if isGzFile
%     gzName=name;
%     name=[AddSlash(pa) nm];
%     [~,~,ex]=fileparts(name);
%     if ~exist(name,'file')  % not already unzipped
%         zipErr=system(['gunzip -kf ' gzName]);
%         if zipErr
%             ok=false;
%         return
%         end;
%     end;
% end;

% Special case of ztiff files
if numel(name)>4 && strcmp(name(end-4:end),'z.tif')
    ex='z.tif';
end;

switch lower(ex)
    case {'.mrc','.st','.mrcs'}
        [m, s]=ReadMRC(name,start,numSlices);
        ok=numel(m)>0;
    case '.tif'
        [m, s]=ReadTiffStack(name,start,numSlices);
    case 'z.tif'
        [m,s,ok]=ReadZTiff(name,start,numSlices);
    case '.stk'  % primitive reader for metamorph stacks
        q=bfopen(name);
        q=q{1};
        m=q{1,1};
        for i=2:size(q,1)
            m(:,:,i)=q{i,1};
        end;
        m=flipud(m);
    otherwise
        if nargout<3
            error(['Unknown file type: ' ex]);
        else
            ok=false;
        end;
end;
if showInfo && numel(m)>0
    whos m
end;
if nargout<1
    clear m % Don't list the huge array.
end;

% if isGzFile % We'll delete the unzipped file
%     disp(['deleting ' name]);
% %     delete(name);
% end;
