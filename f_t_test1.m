% Tests the null hypothesis that two groups with n and m functions do not differ
% Calculates maximum t-values of all possible permutations
%
%_______________________________________________________________________________
% This program performs all possible permutations. This can take a long time for 
% large groups (for example, more than 1 day if both groups consist of 
% more than 10 functions).
% Therefore, there is an alternative program that only performs the
% functional t-test for a limited number of permutations (f_t_test2).
%_______________________________________________________________________________
% Name this file f_t_test1.m
% this file has to be accessible from MATLAB search path; e.g., it can be saved in
% MATLAB Startup Folder: userhome/Documents/MATLAB.
%
% the file f_t_test_T_max.m has to be accessible from MATLAB search path
%_______________________________________________________________________________
% runs after the command
% f_t_test1(f_1,f_2,sampler)
%_______________________________________________________________________________
% f_1 - function group 1
% f_2 - function group 2
% sampler - x axis of the fuctions (e.g. [0:0.1:12])
%
% function groups have to be inserted as anonymus functions
% these ananymus functions have to be crated before the functional t-test
% is startet
%
% Example for the creation of a function group (consisting of gamma functions)
% e.g test1=[10 7 2;10 8 2;4 18 5]
%    for n = 1:length(test1(:,1))
%    f_1{n} =@(t)((test(1)*t.^(test(2)-1).*exp(-t.*test(2)) .* test(2)^test(3) / gamma(test(2))));
%    end
%_______________________________________________________________________________
% following temporary results are specified in the command line:
%   ->      done                - percentage of how many permutations
%                               have been completed
%   ->      p-value             - intermediate result of p-value
%
%
% result: a vector named 'p'
%   -> p(1)     -     p-value of the functional t-test
%   -> p(2)     -     T_org_max
%   -> p(3:end) -     T_max of the computed permutations
%_______________________________________________________________________________
% 
% Martin Segeroth, 2018


function p = f_t_test1(f_1,f_2,sampler)

f_1_length = length(f_1);
f_2_length = length(f_2);
k = f_1_length; 

length_permut_line = f_1_length + f_2_length;
n = length_permut_line;

f = [f_1 f_2];


% compute T of the original groups
orginal_line = zeros(1,n);
orginal_line(1:f_1_length) = 1;

select_f_1 = find(orginal_line);
org_m_f_1 = @(x) 0;
for i = 1:f_1_length
    org_m_f_1 = @(x) ( org_m_f_1(x) + f{select_f_1(i)}(x) );
end
org_m_f_1 = @(x) (1/f_1_length)*( org_m_f_1(x) );

org_v_f_1 = @(x) 0;
for i = 1:f_1_length
    org_v_f_1 = @(x) ( org_v_f_1(x) + (f{select_f_1(i)}(x)-org_m_f_1(x)).^2 );
end
org_v_f_1 = @(x) (1/(f_1_length-1))*( org_v_f_1(x) );

org_f_2 = ones(1,length_permut_line) - orginal_line;
select_f_2 = find(org_f_2);
org_m_f_2 = @(x) 0;
for i = 1:f_2_length
    org_m_f_2 = @(x) ( org_m_f_2(x) + f{select_f_2(i)}(x) );
end
org_m_f_2 = @(x) (1/f_2_length)*( org_m_f_2(x) );

org_v_f_2 = @(x) 0;
for i = 1:f_2_length
    org_v_f_2 = @(x) ( org_v_f_2(x) + (f{select_f_2(i)}(x)-org_m_f_2(x)).^2 );
end
org_v_f_2 = @(x) (1/(f_2_length-1))*( org_v_f_2(x) );

T_org = @(x) ((abs( org_m_f_1(x) - org_m_f_2(x) ))./( sqrt( (1/f_1_length)*org_v_f_1(x) + (1/f_2_length)*org_v_f_2(x) ) ));
        
y = sampler;
R = T_org(y);
T_org_max = max(R);






v = ones(1,n);

str_v = [num2str(v(1))];
for i = 2:n
    str_v = [str_v , num2str(v(i))];
end

str_v_comp = [num2str(1)];
for i = 2:n
    str_v_comp = [str_v_comp , num2str(0)];
end

T_max = [0];

while sum(str_v == str_v_comp) ~= n
    str_v = dec2bin(bin2dec(str_v)-bin2dec('1'));
    
    v = [str2num(str_v(1))];
    for i = 2:n
        v = [v , str2num(str_v(i))];
    end
    
    v_u = v - [1 zeros(1,n-1)];
    
    if (sum(v) == k)
        
        T_max = [ T_max , f_t_test_T_max(f_1,f_2,v,sampler) ];
        
        p_value = sum((T_max > T_org_max))/(length(T_max)-2)
        length(T_max)/(factorial(n)/(factorial(k)*factorial(n-k)))
    
    elseif (sum(v_u) == k)
        
        T_max = [ T_max , f_t_test_T_max(f_1,f_2,v_u,sampler)];
        
        p_value = sum((T_max > T_org_max))/(length(T_max)-2)
        done=length(T_max)/(factorial(n)/(factorial(k)*factorial(n-k)))
    
    end
    
end


p = [p_value, T_org_max, T_max]