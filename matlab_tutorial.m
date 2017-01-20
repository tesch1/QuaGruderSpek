% -*-matlab-*-
% Matlab Tutorial Examples
%
% type these examples into Matlab's command window...
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Command Line
1+1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Variables

% create a variable "a"
a = 1

% show value of variable
a

% change variable value, but dont show result (; makes silent)
a = 2;
a

% variables can have imaginary values
a = 2 + 4i

% variables are matrices

% create a row
x = [ 1 2 3 ]

% create a column (row transposed)
x = [ 1 2 3 ].'

% or, equivalently (use ; to separate rows)
A = [ 1 ; 2 ; 3 ]

% 2D matrices
A = [ 1 2 ; 3 4 ]

% variables can have any name, but dont use "i" !!!
i
i = 2

% accessing just part of a matrix (a sub-matrix)
A(1)

% create larger 2D matrix
A = [1 2 3 ; 4 5 6; 7 8 9]

% accessing sub-matrix
A(2:3, 2:3)

% same thing:
A([2 3], [2 3])

% notation for sequence using ':'
3:8

% a sequence is also a row vector
x = 3:8

% assign sequence to column vector
x = [3:8].'

% the help command
help
help zeros

% creating a square matrix with all zeros
zeros(4)

% creating a non-square matrix with all zeros: 3 rows, 4 columns
zeros(3, 4)

% all ones
ones(3)

% identity matrix
eye(3)
eye(4)

% change the shape of a matrix
A = reshape([1:9], 3, 3)

% Math with Matrices

% scalar addition and subtraction
A + 1
A - 1

% matrix element-wise operations
A - eye(3)

% matrix-matrix multiplication
x = [1 ; 1 ; 1 ]
y = x.';

% matrix-columnn
A*x
% row-matrix
y*A
% square matrix - square matrix
A*A

% element-wise operations (.*) instead of (*)
A .* eye(3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Looping

a = 0;
for n = 1:10
  a = a + n;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Functions

% inline-function definition, f(x) = x^2
f = @(x) (x * x);

% anti-example: f([2 2])

% better:
f = @(x) (x .* x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Plotting

% implicit x-axis
plot(30:40)

% explicit x-axis
plot(0:10, f(0:10))

% labelling the plot (DO THIS ALWAYS!)
xlabel('x-axis')
ylabel('x squared')
title('a good plot')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Numbers & Errors

% Special Numbers

% dividing by zero - Inf
1/0

% using Inf - NaN (Not-a-Number!)
Inf/Inf

% Limitations of numbers
sprintf('%f', 10^30)

% subtraction and addition of vastly different scales gives errors
10^30 - (10^30 + 1)

sprintf('%f', 10^30 + 1)

% also small numbers can be approximations
sprintf('%.30f', 0.1)

