%%% continue command
clc; clear all;
% for a = 1:16
%     if mod(a,5) == 0 %%%% this line says that a, divide by 5, remainder 0
%         continue %%%% skips THIS ITERATION OF LOOP ONLY
%     end
%     
%     a
% end



for b = 1:15
    if b*b > 100
        continue
    end
    b
    b*b
end

%%%% Happy with continue statement.