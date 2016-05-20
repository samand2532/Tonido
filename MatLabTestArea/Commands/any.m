%%%% any command
clc; clear all;

% A = rand(1,15)*3
% 
% if any(A < 0.1)
%     disp('yes')
% else disp('no')
% end


A = [ 1 2 5 2 6 88 9]



for CoutA = 1:5
    if any(A == -1)
        disp('less 0')
        continue
    else A = A -1
    end
end
