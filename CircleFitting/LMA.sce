clear all;
clc;

disp('Problem statement: Fit y= a*(1-exp(-b*x)) following data using LMA with a and b as parameters.  x 0.25 0.75 1.25 1.75 2.25  y 0.28 0.57 0.68 0.74 0.79');
 
disp('Use value of damping factor (when asked) as zero for Gauss Newton algorithm') ;

//Enter your function as y= f(x, [parameters])
function y=func(x,a,b)
    y=a*(1-exp(-b*x));
endfunction

//First column of jacobian matrix, partial derivative w.r.t. parameter a 
function Ja=JacobianA(x,b)
    Ja=1-exp(-b*x)    
endfunction

//Second column of jacobian matrix, partial derivative w.r.t. parameter b 
function Jb=JacobianB(x,a,b)    
    Jb= a*x*exp(-b*x)
endfunction

//Data to fit
xdata= [0.25 0.75 1.25 1.75 2.25]
ydata= [0.28 0.57 0.68 0.74 0.79]

//Initial guesses for parameters, damping factor and no of iterations from user 
a=input("Enter value of a");
b=input("Enter value of b");
lambda=input("Enter value of damping factor");
if lambda==0 then betaV=1
else    
betaV= input("Enter value of gain ratio");
end
numberOfIteration= input("Enter number of iterations to be performed");

IntialGuess= [a ; b]

//The following for loop iterates algorithm for new values of parameters and damping factor
for k=1:numberOfIteration   

xlength= length(xdata);

//Calculates jacobian matrix
for i=1: xlength   
   Z(i,1)=JacobianA(xdata(i),IntialGuess(2,1))      
   Z(i,2)=JacobianB(xdata(i),IntialGuess(1,1),IntialGuess(2,1))
end

// Calculates transpose of Z
ZT= Z'

ZTZ= ZT * Z;

lambdaI= lambda* eye(2,2);

ZTZPlusLI= ZTZ+lambdaI;

//calculates invere of ZTZ +Lambda*I
ZTZPlusLIInverse= inv(ZTZPlusLI);


for i=1: xlength
    D(i,1)= ydata(i)- func((xdata(i)),IntialGuess(1,1),IntialGuess(2,1))
end


ZTD= ZT*D;


Increments= ZTZPlusLIInverse* ZTD;


NewValues= IntialGuess + Increments;
disp(NewValues);

//Calculates D*D
for i=1:xlength
    yMinusyOld(i,1)=ydata(i)-func((xdata(i)),IntialGuess(1,1),IntialGuess(2,1))
    yMinusyOldSquare(i,1)=yMinusyOld(i,1)*yMinusyOld(i,1);
end

sumOld=sum(yMinusyOldSquare);

//Calculates Dnew*Dnew
for i=1:xlength
    yMinusyNew(i,1)=ydata(i)-func((xdata(i)),NewValues(1,1),NewValues(2,1))
    yMinusyNewSquare(i,1)=yMinusyNew(i,1)*yMinusyNew(i,1);
end

sumNew= sum(yMinusyNewSquare);



//Criterion for choosing new value of damping factor for next iteration
if sumOld > sumNew
    lambda= lambda/betaV
else
    lambda= lambda* betaV;
   end    
IntialGuess=[NewValues(1,1); NewValues(2,1)]

end
