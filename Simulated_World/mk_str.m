function [string_out] = mk_str(string_in)
%mk_str Replaces "_" in strings with "\_"
%   Detailed explanation goes here
string_out = [];
for i = 1:length(string_in),
    if string_in(i) == '_',
        string_out = [ string_out '\_' ];
    else
        string_out = [ string_out string_in(i) ];
    end
end

