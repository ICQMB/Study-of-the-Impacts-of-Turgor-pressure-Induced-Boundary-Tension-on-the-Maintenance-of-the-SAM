The CSV_interpreter program is made to work with 
the batchmaker.cpp program, so compile by

g++ batchmaker.cpp -o batchGenerator.out
g++ CSV_interpreter.cpp -o CSV_interpreter.out

and run it with

./CSV_interpreter.out <csv file name>.csv

INPUT: First row index K; CSV file organized as follows: 

-batchmaker_flag_1,flag_value_1, ...,-batchmaker_flag_M,flag_value_M,-test,NAME
-parameter_flag_1,...,-parameter_flag_P
xK1,x12,...,xKP
x(K+1)1,x(k+2)2,...,x(K+1)P
...
x(K+N)1,X(K+N)2,...,X(K+N)P

OUTPUT: N simulations submitted with parameter vectors (Xk1,...,XkP)
submitted via slurm as NAME_K, NAME_(K+1), ... , NAME_(K+N). 

List of batchmaker_flags (Optional marked with *)
-test <NAME>
-hours <##> 
-cores <##>
-nodes <##>
-partition <short/batch/intel>
*-mem <#> (# is in GB, default 2)
*-nodes <#> (defaults to 1).
