function WriteStarFileStruct(s,dataName,fileName,activeFlags)
% function WriteStarFileStruct(s,dataName,fileName,activeFlags)
% function WriteStarFileStruct(s,FID)
% Write a STAR file from the fields of the struct s.  These fields are
% assumed to be numeric or cell arrays (of strings).  The arrays must all
% have the same number of elements, which we'll call nLines.
% The optional nLines x 1 boolean array activeFlags determines which lines
% actually get writen out rather than skipped.
% 
% Give a file handle FID instead of fileName if you want to write several
% structs to a single file.  In this case you'd call
% FID=fopen('file.star','w'); WriteStarFileStruct(..,FID); fclose(FID);
% dataName is the string written after 'data_' on the first line of the
% file.

if ~isa(s,'struct')
    error('Input is not a struct.');
end;

fnames=fieldnames(s);
nFields=numel(fnames);
if nFields<1
    error('Input structure has no fields');
end;
nLines=numel(s.(fnames{1}));

if nargin<4
    activeFlags=true(nLines,1);
end;
%%

if isnumeric(fileName)
    fi=fileName;
else
    fi=fopen(fileName,'w');
end;

fprintf(fi,'\n');
fprintf(fi,'data_%s\n',dataName);
fprintf(fi,'\n');

fprintf(fi,'loop_\n');
for i=1:nFields
    fprintf(fi,'_%s # %d\n',fnames{i},i);
end;

for iLine=1:nLines
    if activeFlags(iLine)  % write the line
        for iField=1:nFields
            x=s.(fnames{iField});
            if isa(x,'numeric')
                fprintf(fi,'%12.6g ',x(iLine));
            elseif isa(x,'char')
                fprintf(fi,'%s ',x);
                if iLine>1
                    error([fnames{iField} ' is not a cell array.']);
                end;
            elseif isa(x,'cell')
                fprintf(fi,'%s ',x{iLine});
            else
                x
                error('Unrecognized datatype');
            end;
        end;
        fprintf(fi,'\n');
    end;
end;
fprintf(fi,'\n');

if ~isnumeric(fileName)
    fclose(fi);
end;
