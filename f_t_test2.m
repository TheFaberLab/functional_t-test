% Tests the null hypothesis that two groups with n and m functions do not differ
% Calculates maximum t-values of a selected number of permutations
%
%_______________________________________________________________________________
% Since it takes a long time to calculate the maximum t-values for all
% permuations if there are large groups (for example, more than 1 day if both 
% groups consist of more than 10 functions), this program offers the possibility 
% to perform the functional t-test for a choosen number of permutations.
%
% There is an alternative program that performs all
% permutations (f_t_test1).
%_______________________________________________________________________________
% name this f_t_test2.m
% this file has to be accessible from MATLAB search path; e.g., it can be saved
% in MATLAB Startup Folder: userhome/Documents/MATLAB.
%_______________________________________________________________________________
% runs after the command
% f_t_test2(f_1,f_2,permutations,max_same_permutation,sampler)
%_______________________________________________________________________________
%                   f_1 - function group 1
%                   f_2 - function group 2
%          permutations - number of permutations to be performed
% max_same_permutation  - Since the permutations are randomly selected, 
%                         a permutation can occur more than once.
%                         This value limits how often a permutation may occur.
%               sampler - x axis of the fuctions (e.g. [0:0.1:12])
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
%   -> already_found            - appears, if a permutation should be 
%                               calculated more often than max_same_permutation
%   -> counter_of_already_found - counts how many permutations have been
%                               choosen more often than max_same_permutation
%   -> T_max_current            - maximum t-value of the current permutation
%
%
% final result: a vector named 'p'
%   ->      p(1)     -     p-value of the functional t-test
%   ->      p(2)     -     T_org_max
%   ->      p(3:end) -     T_max of the computed permutations
%_______________________________________________________________________________
%
% 
% Martin Segeroth, 2018


function p = f_t_test2(f_1,f_2,permutations,max_same_permutation,sampler)

f_1_length = length(f_1);
f_2_length = length(f_2);

length_permut_line = f_1_length + f_2_length;

f = [f_1 f_2];


% compute T of the original groups
orginal_line = zeros(1,length_permut_line);
orginal_line(1:f_1_length) = 1;

select_1 = find(orginal_line);
org_m_1 = @(x) 0;
for i = 1:f_1_length
    org_m_1 = @(x) ( org_m_1(x) + f{select_1(i)}(x) );
end
org_m_1 = @(x) (1/f_1_length)*( org_m_1(x) );

org_v_1 = @(x) 0;
for i = 1:f_1_length
    org_v_1 = @(x) ( org_v_1(x) + (f{select_1(i)}(x)-org_m_1(x)).^2 );
end
org_v_1 = @(x) (1/(f_1_length-1))*( org_v_1(x) );

org_2 = ones(1,length_permut_line) - orginal_line;
select_2 = find(org_2);
org_m_2 = @(x) 0;
for i = 1:f_2_length
    org_m_2 = @(x) ( org_m_2(x) + f{select_2(i)}(x) );
end
org_m_2 = @(x) (1/f_2_length)*( org_m_2(x) );

org_v_2 = @(x) 0;
for i = 1:f_2_length
    org_v_2 = @(x) ( org_v_2(x) + (f{select_2(i)}(x)-org_m_2(x)).^2 );
end
org_v_2 = @(x) (1/(f_2_length-1))*( org_v_2(x) );

T_org = @(x) ((abs( org_m_1(x) - org_m_2(x) ))./( sqrt( (1/f_1_length)*org_v_1(x) + (1/f_2_length)*org_v_2(x) ) ));
        
y = sampler;
R = T_org(y);
T_org_max = max(R);


permut_line = zeros(1,length_permut_line);
zv_x = randperm(length_permut_line);
permut_line(zv_x(1:f_1_length)) = 1;

counter = 1;
select_1 = find(permut_line);
perm_m_1 = @(x) 0;
for i = 1:f_1_length
    perm_m_1 = @(x) ( perm_m_1(x) + f{select_1(i)}(x) );
end
perm_m_1 = @(x) (1/f_1_length)*( perm_m_1(x) );

perm_v_1 = @(x) 0;
for i = 1:f_1_length
    perm_v_1 = @(x) ( perm_v_1(x) + (f{select_1(i)}(x)-perm_m_1(x)).^2 );
end
perm_v_1 = @(x) (1/(f_1_length-1))*( perm_v_1(x) );

perm_2 = ones(1,length_permut_line) - permut_line;
select_2 = find(perm_2);
perm_m_2 = @(x) 0;
for i = 1:f_2_length
    perm_m_2 = @(x) ( perm_m_2(x) + f{select_2(i)}(x) );
end
perm_m_2 = @(x) (1/f_2_length)*( perm_m_2(x) );

perm_v_2 = @(x) 0;
for i = 1:f_2_length
    perm_v_2 = @(x) ( perm_v_2(x) + (f{select_2(i)}(x)-perm_m_2(x)).^2 );
end
perm_v_2 = @(x) (1/(f_2_length-1))*( perm_v_2(x) );

T = @(x) ((abs( perm_m_1(x) - perm_m_2(x) ))./( sqrt( (1/f_1_length)*perm_v_1(x) + (1/f_2_length)*perm_v_2(x) ) ));
        
y = sampler;
R = T(y);
T_max(counter) = max(R);

counter = counter + 1;


str_line = [];
for i = 1:length_permut_line
    str_line = [str_line , num2str(permut_line(i))];
end

values_of_line = 0;
values_of_line = [values_of_line ; bin2dec(str_line) ];
values_of_line = [values_of_line ; 2^length_permut_line ];

permut_matrix = zeros(1,length_permut_line);
permut_matrix = [permut_matrix ; permut_line ];
permut_matrix = [permut_matrix ; ones(1,length_permut_line)];

counter_of_already_found = 0;
k = 1;
while k <= permutations
    permut_line = zeros(1,length_permut_line);
    zv_x = randperm(length_permut_line);
    permut_line(zv_x(1:f_1_length)) = 1;

    str_line = [];
    for i = 1:length_permut_line
        str_line = [str_line , num2str(permut_line(i))];
    end
    
    v_line = bin2dec(str_line);
    
    how_many_equal = sum(values_of_line == v_line);
    
    if ((how_many_equal == 0)+(counter_of_already_found > max_same_permutation)) == 1
        
        select_1 = find(permut_line);
        perm_m_1 = @(x) 0;
        for i = 1:f_1_length
            perm_m_1 = @(x) ( perm_m_1(x) + f{select_1(i)}(x) );
        end
        perm_m_1 = @(x) (1/f_1_length)*( perm_m_1(x) );

        perm_v_1 = @(x) 0;
        for i = 1:f_1_length
            perm_v_1 = @(x) ( perm_v_1(x) + (f{select_1(i)}(x)-perm_m_1(x)).^2 );
        end
        perm_v_1 = @(x) (1/(f_1_length-1))*( perm_v_1(x) );

        perm_2 = ones(1,length_permut_line) - permut_line;
        select_2 = find(perm_2);
        perm_m_2 = @(x) 0;
        for i = 1:f_2_length
            perm_m_2 = @(x) ( perm_m_2(x) + f{select_2(i)}(x) );
        end
        perm_m_2 = @(x) (1/f_2_length)*( perm_m_2(x) );

        perm_v_2 = @(x) 0;
        for i = 1:f_2_length
            perm_v_2 = @(x) ( perm_v_2(x) + (f{select_2(i)}(x)-perm_m_2(x)).^2 );
        end
        perm_v_2 = @(x) (1/(f_2_length-1))*( perm_v_2(x) );
        
        T = @(x) ((abs( perm_m_1(x) - perm_m_2(x) ))./( sqrt( (1/f_1_length)*perm_v_1(x) + (1/f_2_length)*perm_v_2(x) ) ));
        
        y = sampler;
        R = T(y);
        T_max(counter) = max(R);

        
        counter = counter + 1;
        
        T_max_current=T_max(counter-1)
        p_value = sum((T_max > T_org_max))/(counter-1)
        
        
        values_smaller = values_of_line <= v_line;
        where = max(find(values_smaller));

        split_oben = permut_matrix(1:where,:);
        split_unten = permut_matrix((where+1):k+2,:);

        permut_matrix = [ split_oben ; permut_line ; split_unten ];

        split_oben_values = values_of_line(1:where);
        split_unten_values = values_of_line((where+1):k+2);

        values_of_line = [ split_oben_values ; v_line ; split_unten_values ];
        
        counter_of_already_found = 0;
    else
        k = k-1;
        
        already_found = strcat('permutation alredy found for',{' '},num2str(max_same_permutation),' times')
        counter_of_already_found = counter_of_already_found +1
    end
    
    
    
    done=k/permutations
    k = k+1;
end

p = [p_value, T_org_max, T_max]
