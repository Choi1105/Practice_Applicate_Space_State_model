%**************************************************************************
% By Martin M. Andreasen. 
% This program display's the content in the matrix AA(1:N1,1:N2,1:N3,1:N4)

NAME                = 'fyp';        %The name for the matrix in the output file

AA                  = eval(NAME); %The matrix of focus

DISP_matlab         = 1;          %Output for Matlab

DISP_fortran        = 0;          %Output for Fortran
NUM_pr_line         = 130;        %The number of characters pr. linie in fortran - max is 132

% The results are saved in the files Matlab_NAME.dat and/or Fortran_NAME_dat

% Log:
% 23.10 2011 so it can deal with double digits in the matrices
% 18.02 2012 corrected a bug in the loop for the "remaining characters" by 
%            replacing Num_char-1 by Num_char-4 (around line 124)
%**************************************************************************
% Clears the command windue
clc;

if ndims(AA) == 2
    N1 = size(AA,1);
    N2 = size(AA,2);
    N3 = 1;
    N4 = 1;
elseif ndims(AA) == 3
    N1 = size(AA,1);
    N2 = size(AA,2);
    N3 = size(AA,3); 
    N4 = 1;
    order_app = 2;
elseif ndims(AA) == 4
    N1 = size(AA,1);
    N2 = size(AA,2);
    N3 = size(AA,3);
    N4 = size(AA,4);
    order_app = 3;    
end
 
% For Matlab
if DISP_matlab == 1

% The results are saved in the file SAVE_matlab - results are appended to this file
SAVE_matlab = ['Matlab_' NAME '.txt']; 
diary(SAVE_matlab);

disp(['% START DISPLAYING ' NAME])   
for g=1:N4
    for k=1:N3
        for j=1:N2
            for i=1:N1
                if AA(i,j,k,g) ~= 0
                    if ndims(AA) == 2
                        disp([NAME '(' num2str(i) ',' num2str(j) ')=' char(AA(i,j,k,g)) ';']);
                    elseif ndims(AA) == 3
                        disp([NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ')=' char(AA(i,j,k,g)) ';']);
                    elseif ndims(AA) == 4
                        disp([NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ',' num2str(g) ')=' char(AA(i,j,k,g)) ';']);
                    end
                end
            end
        end
    end
end
disp(['% END DISPLAYING ' NAME])
diary off;

end



% For Fortran 
if DISP_fortran == 1

% The results are saved in the file SAVE_fortran - results are appended to this file
SAVE_fortran =['Fortran_' NAME '.txt'];
diary(SAVE_fortran);
    
disp(['! START DISPLAYING ' NAME])   
Index_start = 110;
for g=1:N4
    for k=1:N3
        for j=1:N2
            for i=1:N1
                if AA(i,j,k,g) ~= 0
                    % Convert the symbolic expression in AA() into a character string
                    TEMP = char(AA(i,j,k,g));

                    % Replacing ^ with **,exp with EXP and log with LOG
                    TEMP = strrep(TEMP,'^','**');
                    TEMP = strrep(TEMP,'exp','EXP');
                    TEMP = strrep(TEMP,'log','LOG');

                    % Changing double digits from 20 to 20._8
                    for m1=1:9
                        for m2=0:9
                            TEMP = strrep(TEMP,[num2str(m1) num2str(m2)],[[num2str(m1) num2str(m2)] '._8']);
                        end
                    end
                    
                    % Furthermore, all constants are changed from 5 to 5._8
                    % where 8 is double precision in fortran. This must be done without changing
                    % for instance PHI2 to PHI2._8
                    Num_char = size(TEMP,2);
                    Index = 1;
                    clear TEMP2;
                    TEMP_new = TEMP;
                    % For the first character in TEMP
                    for m = 1:9
                        if TEMP(1) == num2str(m)
                            Index = Index + 3;
                            TEMP2(1:Index) = [TEMP_new(1) '._8'];
                            TEMP_new = [TEMP2(1:Index) TEMP(2:end)];
                        end
                    end
                    
                    % For the remaining characters in TEMP
                    for n = 2:Num_char
                        Index = Index + 1;
                        for m = 1:9
                            if n < Num_char-4
                                if TEMP(n)== num2str(m) && TEMP(n+1) ~= '0' && TEMP(n+1) ~= '1' && TEMP(n+1) ~= '2' ...
                                        && TEMP(n+1) ~= '3' && TEMP(n+1) ~= '4' && TEMP(n+1) ~= '5' && TEMP(n+1) ~= '6'...
                                        && TEMP(n+1) ~= '7' && TEMP(n+1) ~= '8' && TEMP(n+1) ~= '9'...
                                        && TEMP(n+3) ~= '8' && TEMP(n+4) ~= '8' ...
                                        && TEMP(n-1) ~= 'A' && TEMP(n-1) ~= 'I' && TEMP(n-1) ~= 'F' ...
                                        && TEMP(n-1) ~= 'X' && TEMP(n-1) ~= 'S' && TEMP(n-1) ~= 'p' ...
                                        && TEMP(n-1) ~= 'a' && TEMP(n-1) ~= 'i' && TEMP(n-1) ~= 'f' ...
                                        && TEMP(n-1) ~= 'x' && TEMP(n-1) ~= 's' && TEMP(n-1) ~= '_'
                                    Index = Index + 3;
                                    TEMP2(1:Index) = [TEMP_new(1:Index-3) '._8'];
                                    TEMP_new = [TEMP2(1:Index) TEMP(n+1:end)];
                                end
                            else
                                if TEMP(n)== num2str(m) ...
                                        && TEMP(n-1) ~= 'A' && TEMP(n-1) ~= 'I' && TEMP(n-1) ~= 'F' ...
                                        && TEMP(n-1) ~= 'X' && TEMP(n-1) ~= 'S' && TEMP(n-1) ~= 'p' ...
                                        && TEMP(n-1) ~= 'a' && TEMP(n-1) ~= 'i' && TEMP(n-1) ~= 'f' ...
                                        && TEMP(n-1) ~= 'x' && TEMP(n-1) ~= 's' && TEMP(n-1) ~= '_'
                                    Index = Index + 3;
                                    TEMP2(1:Index) = [TEMP_new(1:Index-3) '._8'];
                                    TEMP_new = [TEMP2(1:Index) TEMP(n+1:end)];
                                end
                            end
                        end
                    end

                    % Udating the character expression
                    TEMP = TEMP_new;

                    % Recounting the number of characters for the i,j,k'th element
                    Num_char = size(TEMP,2);

                    % If the number of characters for the i,j,k'th element is large than
                    % 120 then the characters are broken up into more lines by using &
                    if Num_char > 120
                        if ndims(AA) == 2
                            disp([NAME '(' num2str(i) ',' num2str(j) ')=' TEMP(1:110) '&']);
                        elseif ndims(AA) == 3
                            disp([NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ')=' TEMP(1:110) '&']);
                        elseif ndims(AA) == 4
                            disp([NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ',' num2str(g) ')=' TEMP(1:110) '&']);                            
                        end

                        % The number of lines required
                        Num_loop = floor((Num_char-110)/NUM_pr_line);
                        if Num_loop*NUM_pr_line == Num_char-110
                            Num_loop = Num_loop - 1;
                        end
                        Index_start = 110;

                        for l = 1:Num_loop
                            Index_slut = Index_start + NUM_pr_line;
                            disp(['&' TEMP(Index_start+1:Index_slut) '&' '    !' num2str(l)]);
                            Index_start = Index_start + NUM_pr_line;
                        end
                        disp(['&' TEMP(Index_start+1:Num_char)]);
                        disp('  ');
                    else
                        if ndims(AA) == 2
                            disp([NAME '(' num2str(i) ',' num2str(j) ')=' TEMP]);
                        elseif ndims(AA) == 3
                            disp([NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ')=' TEMP]);
                        elseif ndims(AA) == 4
                            disp([NAME '(' num2str(i) ',' num2str(j) ',' num2str(k) ',' num2str(g) ')=' TEMP]);                            
                        end

                    end

                    % Resetting the index
                    Index_start = NUM_pr_line;
                end
            end
        end
    end
end
disp(['! END DISPLAYING ' NAME])
diary off;
end

    