% Is used by the Matlab function f_t_test1.m
% Determines the maximal t-value of the actual permutation
%_______________________________________________________________________________
% Name this file f_t_test_T_max.m
% this file has to be accessible from MATLAB search path; e.g., it can be saved in
% MATLAB Startup Folder: userhome/Documents/MATLAB.
%_______________________________________________________________________________
% 
% Martin Segeroth, 2018

function T_max = f_t_test_T_max(f_1,f_2,permut_zeile,sampler)

f_1_length = length(f_1);
f_2_length = length(f_2);
k = f_1_length;

length_permut_line = f_1_length + f_2_length;
n = length_permut_line;

f = [f_1 f_2];

select_1 = find(permut_zeile);
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

perm_2 = ones(1,length_permut_line) - permut_zeile;
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
T_max = max(R);