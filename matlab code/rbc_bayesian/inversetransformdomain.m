function parameter=inversetransformdomain(x)

% to restrict paramters into the reasonable domain

% parameter([1 3 5 7]) = log(x([1 3 5 7]));       % to make them positive
parameter([1 3 5 6 7]) = log(x([1 3 5 6 7]));       % to make them positive
parameter([2 4])     = log((1-x([2 4]))./x([2 4])); % to make them lie between 0 and 1
% parameter(6)         = x(6);