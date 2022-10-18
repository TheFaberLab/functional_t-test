# Functional t-test
Tests the null hypothesis that two groups of functions do not differ.

# Overview
Two groups of functions, consisting of n1 and n2 functions can be compared using the permutation based functional t-test. For a detailed description of the functional t-test see: <br><br>
[1] Ramsay, J.O., Hooker, G., Spencer, G., 2009. Functional Data Analysis with R and MATLAB, first ed. Springer, Dordrecht, Heidelberg, London, New York <br> <br>
[2] Lambers et al Neuroimage. 2020 Mar;208:116446. doi: 10.1016/j.neuroimage.2019.116446.
https://www.sciencedirect.com/science/article/pii/S1053811919310377?via%3Dihub

# Software Requirements
The script was developed using MATLAB 2018b. 

# Details
There are two versions of the functional t-test (f_t_test1 und  f_t_test2). <br>
f_t_test1 performs all possible permutations. This can take a long time for large groups (for example, more than one day if both groups consist of more than 10 functions). Therefore, there is the program f_t_test2 that only performs the functional t-test for a limited number of permutations. <br>
If the script f_t_test1 is used, the script f_t_test_T_max.m must also be on the MATLAB path. <br><br>
Detailed instructions on how to use the scripts are provided in comments at the beginning of the scripts.

# How to cite the functional t-test
These scripts were created by your colleagues. Please the article in which the functional t-test was described:<br>
Lambers et al Neuroimage. 2020 Mar;208:116446. doi: 10.1016/j.neuroimage.2019.116446. <br>
https://www.sciencedirect.com/science/article/pii/S1053811919310377?via%3Dihub


