clc;  
clear all;

 disp('Problem statement: Fit y=sqrt(R*R-(x-h)*(x-h))+k for following data using LMA with R, h and k as parameters. x [-2.4164 -1.55 -2.7251 -0.5035 3 3.031 1.7677 2.01525 0.5643 0 0 -2.192 -2.1213] y [1.692 2.6846 -0.9918 -2.8559 0 -1.75 -1.7677 2.01525 3.2006 3.5 -2.8 -2.192 2.1213 ]');

disp('Use value of damping factor (when asked) as zero for Gauss Newton algorithm') ;

disp('Use damping factor and gain ratio as 10 and 1-10 respectively for converging results') ;

//Data to fit
xdata= [-2.4164 -1.55 -2.7251 -0.5035 3 3.031 1.7677 2.01525 0.5643 0 0 -2.192 -2.1213]
ydata= [1.692 2.6846 -0.9918 -2.8559 0 -1.75 -1.7677 2.01525 3.2006 3.5 -2.8 -2.192 2.1213 ] 
ylength= length(ydata);

for h=1:ylength
    if ydata(h)==0
        signdecider(h)=1
    else
    signdecider(h)= ydata(h)/ abs(ydata(h))
    end
end

//Enter your function as y= f(x, [parameters])
function y=func(x,R,h,k)   
        y= abs(sqrt(R*R-(x-h)*(x-h)))+k     
endfunction

//First column of jacobian matrix, partial derivative w.r.t. parameter R 
function JR=JacobianR(x,R,h)
    JR= R/abs(sqrt(R*R-(x-h)*(x-h)))  
endfunction

//Second column of jacobian matrix, partial derivative w.r.t. parameter h 
function Jh=JacobianH(x,R,h)
    Jh= (x-h)/abs(sqrt(R*R-(x-h)*(x-h)))  
endfunction

//Third column of jacobian matrix, partial derivative w.r.t. parameter k 
function Jk=JacobianK()
    Jk= 1 
endfunction


//Initial guesses for parameters, damping factor and no of iterations from user 
R=input("Enter value of R");
h=input("Enter value of h");
k=input("Enter value of k");
lambda=input("Enter value of damping factor");
if lambda==0 then betaV=1
else 
    betaV= input("Enter value of beta");
end

numberOfIteration= input("Enter number of iterations to be performed");

IntialGuess= [R ; h; k]

//The following for loop iterates algorithm for new values of parameters and damping factor
for w=1:numberOfIteration   

xlength= length(xdata);

//Calculates jacobian matrix
for i=1: xlength   
   Z(i,1)= signdecider(i)*(JacobianR(xdata(i),IntialGuess(1,1),IntialGuess(2,1)))     
   Z(i,2)=signdecider(i)* (JacobianH(xdata(i),IntialGuess(1,1),IntialGuess(2,1)))
   Z(i,3)=signdecider(i)* (JacobianK())
end

// Calculates transpose of Z
ZT= Z'

ZTZ= ZT * Z;


lambdaI= lambda* eye(3,3);

ZTZPlusLI= ZTZ+lambdaI;

//calculates invere of ZTZ +Lambda*I
ZTZPlusLIInverse= inv(ZTZPlusLI);


for i=1: xlength
    functionvalue= func(xdata(i),IntialGuess(1,1),IntialGuess(2,1),IntialGuess(3,1))- IntialGuess(3,1)
    withsign= signdecider(i)* functionvalue
    value= withsign + IntialGuess(3,1)
    D(i,1)= ydata(i)- value
end


ZTD= ZT*D;

Increments= ZTZPlusLIInverse* ZTD;


NewValues= IntialGuess + Increments;
disp(NewValues);

//Calculates D*D
for i=1:xlength
    functionvalue= func(xdata(i),IntialGuess(1,1),IntialGuess(2,1),IntialGuess(3,1))- IntialGuess(3,1)
    withsign= signdecider(i)* functionvalue
    value= withsign + IntialGuess(3,1)
    
    yMinusyOld(i,1)=ydata(i)-value
    
    yMinusyOldSquare(i,1)=yMinusyOld(i,1)*yMinusyOld(i,1);
end

sumOld=sum(yMinusyOldSquare);

//Calculates Dnew*Dnew
for i=1:xlength
    
    functionvalue= func(xdata(i),NewValues(1,1),NewValues(2,1),NewValues(3,1))- NewValues(3,1)
    withsign= signdecider(i)* functionvalue
    value= withsign + NewValues(3,1)
    
    yMinusyNew(i,1)=ydata(i)-value
    yMinusyNewSquare(i,1)=yMinusyNew(i,1)*yMinusyNew(i,1);
end

sumNew= sum(yMinusyNewSquare);



//Criterion for choosing new value of damping factor for next iteration
if sumOld > sumNew
    lambda= lambda/betaV
else
    lambda= lambda* betaV;
   end    
IntialGuess=[NewValues(1,1); NewValues(2,1); NewValues(3,1)]

end
