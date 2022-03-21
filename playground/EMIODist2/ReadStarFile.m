function [blockNames,blockData,ok]=ReadStarFile(name)
% function [blockNames,blockData,ok]=ReadStarFile(name)
% Read a Relion .star file and put its contents into two cell arrays.
% blockNames is a cell array of the blockNames.  For example here are
% values from reading a post_run.star file.
% 
% >> blockNames{1} =
%   1×12 char array
% data_general
% >> blockData{1} =
%   struct with fields:
%                rlnFinalResolution: 7.0800
%       rlnBfactorUsedForSharpening: -215.7768
%         rlnFittedSlopeGuinierPlot: -53.9442
% >> blockNames{3} =
%   1×12 char array
% data_guinier
%  >> blockData{3} =
%   struct with fields:
%             rlnResolutionSquared: [51×1 double]
%         rlnLogAmplitudesOriginal: [51×1 double]
%                  rlnParticleName: [51×1 cell] % cell array of strings

% cd('/Users/fred/EMWork/relion13_tutorial/betagal/PrecalculatedResults/Refine3D')
% 
% name='post_run1.star';

blockNames={};
blockData={};

commentMarkers={'#'};

ok=exist(name,'file');
if ~ok
    error(['STAR file not found: ' name]);
end;
fi=fopen(name);

nLines=0;
C=cell(1,1);

% Load the whole file into the cell array C, handling comments
disp('loading...');
while ~feof(fi)
    line=fgetl(fi);
    p=strfind(line,commentMarkers);
    hasComment = numel(p)>0;
    if hasComment
        line(p(1):end)=[];
    end;
    % Lines that only contain comments are treated as blank lines.
    nLines=nLines+1;  % Count this line
    if numel(line)<1
        C{nLines}={{}}; % Count as blank
    else
        C{nLines}=textscan(line,'%s');
    end;
end;
fclose(fi);
%% ------------------
% C is a cell array {1,1} (a single '%s' is picked up) containing a cell
% array {nc,1} where nc is the number of tokens in the line.

nBlocks=0;
P=1;  % line pointer

disp('scanning...');
while P<=numel(C) % loop through all the entries
    
    % skip blank lines
    while numel(C{P})<1 || numel(C{P}{1})<1
        P=P+1;
        if P>numel(C)
            return  % exit the function.
        end;
    end;
    
    % Get the block name
    % It starts with 'data_'
    if strncmp(C{P}{1}{1},'data_',5)  % data block
        nBlocks=nBlocks+1;
        blockNames{nBlocks,1}=C{P}{1}{1};  % whole string data_xxx
    else
        error(['''data_'' expected at line ' num2str(nLines)])
    end;
    P=P+1;
    
    % skip blank lines
    while numel(C{P})<1 || numel(C{P}{1})<1
        P=P+1;
    end;
    if P>numel(C)  % reached the end of the file
        error(['End of file. Expected data at line ' num2str(nLines)]);
    end;
    
    loopMode=strcmpi(C{P}{1}{1},'loop_');    
    if loopMode
        P=P+1;
        % skip blank lines
        while numel(C{P})<1 || numel(C{P}{1})<1
            P=P+1;
        end;
        if P>numel(C)
            error(['End of file.  Expected field name at line ' num2str(nLines)]);
        end;
    end;
    
    %  pick up fieldnames
    nFields=0;
    fieldNames=cell(0,1);
    fieldVals=cell(1,0);
    while numel(C{P})>0 && numel(C{P}{1})>0 && numel(C{P}{1}{1})>0 && C{P}{1}{1}(1)=='_' % begins with underscore
        nFields=nFields+1;
        fieldNames{nFields}=C{P}{1}{1}(2:end);
        if ~loopMode
            if numel(C{P})>0 && numel(C{P}{1})>1
                fieldVals{1,nFields}=C{P}{1}{2};
            else
                fieldVals{1,nFields}='';
            end;
        end;
        P=P+1;
    end;
    %%
    nRows=1;
    if loopMode  % Now the values follow immediately after the fieldnames
        nRows=0;
        fieldVals=cell(0,nFields);
        while numel(C{P}{1})>=nFields
            nRows=nRows+1;
            fieldVals(nRows,:)=C{P}{1}';
            P=P+1;
        end;
    end;
    
%     Convert fieldVals to numeric when possible
    q=struct;
    for i=1:nFields
        fn=fieldNames{i};
        numericFVs=str2double(fieldVals(:,i));
        if all(~isnan(numericFVs))
            q.(fn)=numericFVs;
        else
            if nRows==1
                q.(fn)=fieldVals{1,i}; % field is a string
            else
                q.(fn)=fieldVals(:,i); % field is a cell array
            end;
        end;
    end;
   
    blockData{nBlocks,1}=q;
end;
disp('done.');
return;

