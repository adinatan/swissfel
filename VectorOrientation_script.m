% VectorOrientation_script
% turns all vectors into columns similar to a=a(:)
%clc
ToScreen=false;
vars=who;
try % switch to prevent text in command window
    if ~ToScreen
        for ii=1:numel(vars)
            if eval(['(isnumeric(' num2str(vars{ii}) ') | islogical(' num2str(vars{ii}) '))']) && eval(['numel(' num2str(vars{ii}) ')>1']) && eval(['ndims(' num2str(vars{ii}) ')==2'])  % ensure it is numeric 1 or 2-D variable
                if eval(['length(' num2str(vars{ii}) ')==size(' num2str(vars{ii}) ',2)']) && eval(['size(' num2str(vars{ii}) ',1)==1']) % isrow
                    % Convert all Row vectors to colums
                    %
                    eval([num2str(vars{ii}) '=' num2str(vars{ii}) ''';']);
                    %
                elseif eval(['length(' num2str(vars{ii}) ')==size(' num2str(vars{ii}) ',1)']) && eval(['size(' num2str(vars{ii}) ',2)==1']) % iscolumn
                    % Convert all coulumn vectors into rows
                    % disp([num2str(vars{ii}) ' is a column vector'])
                    % eval([num2str(vars{ii}) '=' num2str(vars{ii}) ''';']);
                end
            end
        end
    else %  lots of notifications in command window
        disp(' ')
        for ii=1:numel(vars)
            if eval(['(isnumeric(' num2str(vars{ii}) ') | islogical(' num2str(vars{ii}) '))']) && eval(['numel(' num2str(vars{ii}) ')>1']) && eval(['ndims(' num2str(vars{ii}) ')==2'])  % ensure it is numeric 1 or 2-D variable
                if eval(['length(' num2str(vars{ii}) ')==size(' num2str(vars{ii}) ',2)']) && eval(['size(' num2str(vars{ii}) ',1)==1']) % isrow
                    % Convert all Row vectors to colums
                    disp([num2str(vars{ii}) ' is a row vector'])
                    eval([num2str(vars{ii}) '=' num2str(vars{ii}) ''';']);
                    disp([num2str(vars{ii}) ' is now a column vector'])
                elseif eval(['length(' num2str(vars{ii}) ')==size(' num2str(vars{ii}) ',1)']) && eval(['size(' num2str(vars{ii}) ',2)==1']) % iscolumn
                    % Convert all coulumn vectors into rows
                    disp([num2str(vars{ii}) ' is a column vector'])
                    %eval([num2str(vars{ii}) '=' num2str(vars{ii}) ''';']);
                end
            end
        end
    end
catch
    disp('Something failed in vector orientation')
end

clear vars ToScreen;