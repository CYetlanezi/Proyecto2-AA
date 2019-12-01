clear; clc; clear all;
%Proyecto II - Análisis Aplicado
%Otoño 2019

% Parámetros centrales
%c1 = 10?4
%c2 = 0.99 
%? = 0.1 
%tol = 10?5
%?max = 1.5
%-- LMBGFS m ? {1, 3, 5, 17, 29}

tol = 1e-5;
radioMax = 1.5;

maxBFGS = 200;
maxRCsr1 = 200;
maxLMBGFS = 1000;

% Función Rosenbrock Extended
% Ref: CUTE UnconstrainedOptimizationTestFunctionsCollection.pdf
% x0 = x0 = [?1.2,1,?,?1.2,1] c=100. 
%n ? {2, 10, 100, 200, 1000}

display(" Función Rosenbrock Extended ");

[f, Df, ~] = fExtendRosenbrock();
x0_2 = [-1.2 ; 1]; 

%-----
display(" n = 2");

[x11, iter11] = lineBGFS( f, x0_2, tol, maxBFGS )
[x21, iter21] = mRCsr1( f, x0_2, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[x311, iter311] = lineLMBGFS( f, x0_2, tol, maxLMBGFS ,1)
display(" m = 3");
[x313, iter313] = lineLMBGFS( f, x0_2, tol, maxLMBGFS ,3)
display(" m = 5");
[x315, iter315] = lineLMBGFS( f, x0_2, tol, maxLMBGFS ,5)
display(" m = 17");
[x3117, iter3117] = lineLMBGFS( f, x0_2, tol, maxLMBGFS ,17)
display(" m = 29");
[x3129, iter3129] = lineLMBGFS( f, x0_2, tol, maxLMBGFS ,29)

%-----

display(" n = 10");
x0_10 = repmat(x0_2,5,1);
[x12, iter12] = lineBGFS( f, x0_10, tol, maxBFGS )
[x22, iter22] = mRCsr1( f, x0_10, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[x321, iter321] = lineLMBGFS( f, x0_10, tol, maxLMBGFS ,1)
display(" m = 3");
[x323, iter323] = lineLMBGFS( f, x0_10, tol, maxLMBGFS ,3)
display(" m = 5");
[x325, iter325] = lineLMBGFS( f, x0_10, tol, maxLMBGFS ,5)
display(" m = 17");
[x3217, iter3217] = lineLMBGFS( f, x0_10, tol, maxLMBGFS ,17)
display(" m = 29");
[x3229, iter3229] = lineLMBGFS( f, x0_10, tol, maxLMBGFS ,29)

%-----

display(" n = 100");
x0_100 = repmat(x0_2,50,1);
[x13, iter1] = lineBGFS( f, x0_100, tol, maxBFGS )
[x23, iter2] = mRCsr1( f, x0_100, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[x331, iter331] = lineLMBGFS( f, x0_100, tol, maxLMBGFS ,1)
display(" m = 3");
[x333, iter333] = lineLMBGFS( f, x0_100, tol, maxLMBGFS ,3)
display(" m = 5");
[x335, iter335] = lineLMBGFS( f, x0_100, tol, maxLMBGFS ,5)
display(" m = 17");
[x3317, iter3317] = lineLMBGFS( f, x0_100, tol, maxLMBGFS ,17)
display(" m = 29");
[x3329, iter3329] = lineLMBGFS( f, x0_100, tol, maxLMBGFS ,29)

%------

display(" n = 200");
x0_200 = repmat(x0_2,100,1);

[x14, iter14] = lineBGFS( f, x0_200, tol, maxBFGS )
[x24, iter24] = mRCsr1( f, x0_200, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[x341, iter341] = lineLMBGFS( f, x0_200, tol, maxLMBGFS ,1)
display(" m = 3");
[x343, iter343] = lineLMBGFS( f, x0_200, tol, maxLMBGFS ,3)
display(" m = 5");
[x345, iter345] = lineLMBGFS( f, x0_200, tol, maxLMBGFS ,5)
display(" m = 17");
[x3417, iter3417] = lineLMBGFS( f, x0_200, tol, maxLMBGFS ,17)
display(" m = 29");
[x3429, iter3429] = lineLMBGFS( f, x0_200, tol, maxLMBGFS ,29)

%------

display(" n = 1000");
x0_1000 = repmat(x0_2,500,1);
[x15, iter15] = lineBGFS( f, x0_1000, tol, maxBFGS )
[x25, iter25] = mRCsr1( f, x0_1000, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[x351, iter351] = lineLMBGFS( f, x0_1000, tol, maxLMBGFS ,1)
display(" m = 3");
[x353, iter353] = lineLMBGFS( f, x0_1000, tol, maxLMBGFS ,3)
display(" m = 5");
[x355, iter355] = lineLMBGFS( f, x0_1000, tol, maxLMBGFS ,5)
display(" m = 17");
[x3517, iter3517] = lineLMBGFS( f, x0_1000, tol, maxLMBGFS ,17)
display(" m = 29");
[x3529, iter3529] = lineLMBGFS( f, x0_1000, tol, maxLMBGFS ,29)


%-------------------------------------------------
%-------------------------------------------------

% Función DIXMAANA -- Letra del equipo G
% Ref: CUTE UnconstrainedOptimizationTestFunctionsCollection.pdf
% z0 = z0 = [2,2,?,2,2] 
% n ? {200, 1000} 

[f, Df, ~] = fDixmaanaG();
z0_2 = [2 ; 2];

display(" n = 200");
z0_200 = repmat(z0_2,100,1);
[z11, iter11] = lineBGFS( f, z0_200, tol, maxBFGS )
[z21, iter21] = mRCsr1( f, z0_200, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[z311, ziter311] = lineLMBGFS( f, z0_200, tol, maxLMBGFS ,1)
display(" m = 3");
[z313, ziter313] = lineLMBGFS( f, z0_200, tol, maxLMBGFS ,3)
display(" m = 5");
[z315, ziter315] = lineLMBGFS( f, z0_200, tol, maxLMBGFS ,5)
display(" m = 17");
[z317, ziter3117] = lineLMBGFS( f, z0_200, tol, maxLMBGFS ,17)
display(" m = 29");
[z3129, ziter1429] = lineLMBGFS( f, z0_200, tol, maxLMBGFS ,29)

%----

display(" n = 1000");
z0_1000 = repmat(z0_2,500,1);
[z12, ziter12] = lineBGFS( f, z0_1000, tol, maxBFGS )
[z22, ziter22] = mRCsr1( f, z0_1000, maxRCsr1, radioMax) %con la misma tolerancia

display(" m = 1");
[z321, ziter321] = lineLMBGFS( f, z0_1000, tol, maxLMBGFS ,1)
display(" m = 3");
[z323, ziter323] = lineLMBGFS( f, z0_1000, tol, maxLMBGFS ,3)
display(" m = 5");
[z325, ziter325] = lineLMBGFS( f, z0_1000, tol, maxLMBGFS ,5)
display(" m = 17");
[z327, ziter3217] = lineLMBGFS( f, z0_1000, tol, maxLMBGFS ,17)
display(" m = 29");
[z3229, ziter3229] = lineLMBGFS( f, z0_1000, tol, maxLMBGFS ,29)


