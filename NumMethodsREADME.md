*********************************************************************************************************************************************************************************
Firstly let's tackle down all our options for numerical approximations.
  
  The root approximating algorithms are implemented here in two ways.
    Iteration Based: 
      The last parameter of the function is interpreted as the iteration count.
    Error Based: 
      The last parameter of the function is interpreted as the maximum error submissible for the final value,
      given in as the order to be base 0.1. 
      This allows one to directly control the error cap by specifying the number of decimals to be accurate upto.
      Error here is implemented as the difference of values from one iteration to the next.
  
  [
   Since the third parameter depends on what implementation is used(Iterative/Error-Prone), it is left out in the description
   The terminating conditions vary based off of this.
   In iteration based approximations, the terminating condition would be when there are no more iterations left/
   In error based approximations, the terminating condtion would be when the error from the last term falls under the maximum permitable error.
  ]
  
  1.BisectionMethod:
      It takes in the function in question and the initial interval.
      The mean of the interval(<x>) is computed and the functional projection of that value is realized.
        
          <x> = (lower bound + upper bound) / 2
      
      If f(<x>).f(lower bound) < 0, the interval is modified to have <x> as the upper bound.
      The corresponding interval is then Bisected again recursively.
      
      If f(<x>).f(upper bound) < 0, the interval is modified to have <x> as the lower bound.
      The corresponding interval is then Bisected again recursively.
      
      If both of these cases do not occur, the algorithm assumes that the accuracy is 100% and returns <x> and it's functional value.
      The reasoning for this is that if either one or both of products of functional values are zero,
      One or more of {<x>, lower bound, upper bound} are roots. If both products turn out to be positive,
      then the algorithm won't be able to move to another point and stops the reccursion then and there.
      
      It returns the resulting approximate root and it's functional value.
   
 2.NewtonsMethod:
    It takes in the function in question and the initial value of independent variable. Here the derivative is numerically approximated.
    It can also take the function in question, it's derivative, and the initial value of independent variable.
    
    If f'(last x) is not 0, 
    The next x value is directly computed by taking
       
       next x = last x - f(last x) / f'(last x)
    
    Otherwise the algorithm assumes 100% accuracy and terminates the recurrsion, returning the last x and it's functional value
    
    The resulting new value of independent variable is then passed in as the initial independent variable for the next recursion of NewtonsMethod
    
    It returns the resulting approximate root and it's functional value.
 
 3.SecantMethod:
      It takes in the function in question and two independent variables.
      
      If the functional values of both variables are different,
      the next approximation is computed.
      (Otherwise it assumes the 2nd x is the next x).
          
          next x = 2nd x - f(2nd x)(2nd x - 1st x)/(f(2nd x) - f(1st x))
      
      The new x is then reccursively parameterized as the 2nd x and the current 2nd x as the 1st x in SecantMethod
      
      It returns the resulting approximate root and it's functional value.
 
 4.RegulaFalsiMethod:
      It takes in the function in question and the initial interval.
      
      If the functional values of both variables are different,
      the next approximation is computed.
      (Otherwise it assumes the 2nd x is the next x).
          
          next x = upper bound - f(upper bound)(upper bound - lower bound)/(f(upper bound) - f(lower bound))
      
      If f(next x).f(lower bound) < 0, the interval is modified to have next x as the upper bound.
      The corresponding interval is then Bisected again recursively.
      
      If f(next x).f(upper bound) < 0, the interval is modified to have next x as the lower bound.
      The corresponding interval is then Bisected again recursively.
      
      If both of these cases do not occur, the algorithm assumes that the accuracy is 100% and returns next x and it's functional value.
      The reasoning for this is that if either one or both of products of functional values are zero,
      One or more of {next x, lower bound, upper bound} are roots. If both products turn out to be positive,
      then the algorithm won't be able to move to another point and stops the reccursion then and there.
      
      It returns the resulting approximate root and it's functional value.

5.FixedPointMethod:
      It takes in the function in question and the initial independent variable.
      It then computes the functional value at that point.
   
          f(x) = new x
      
      The functional value is then passed on as the initial independent variable in FixedPointMethod reccursively.
      
      It returns the resulting approximate root and it's functional value.
*********************************************************************************************************************************************************************************
Now onto the interpolation algorithms.

  1.NewtonsInterpolation:Equally spaced inputs
      It takes in a collection of values and the smallest interval of input variables.
      It then computes the class width of the dataset.
          
          h = (upper bound - lower bound) / (number of values - 1)

      It then proceeds to generate the difference table and extract only the upper row. Lets call this collection the DifferenceLadder.
      The resulting polynomial is given the first value in the DifferenceLadder(y(lower bound)).
      
      The first term is set to the polynomial p(x) = 1.
      
      In an interative loop, the term is updated to contain the root given by
         
         root = lower bound + h * (iteration number - 1)
          
      The value of the resulting polynomial is updated as the sum of itself with the ith term 
      multiplied by the corresponding DifferenceLadder and devided by the ith factorial.
      
      It return the resulting polynomial.
      
 2.NewtonsInterpolation:Ambigously spaced inputs
      It takes in a collection of 2d coordinates of the form (x, y).
    
      It then proceeds to generate the devided difference table and extract only the upper row. Lets call this collection the DifferenceLadder.
      The resulting polynomial is given the first value in the DifferenceLadder(y(lower bound)).
      
      The first term is set to the polynomial p(x) = 1.
      
      In an interative loop, the ith term is updated to contain the root given by
         
         root = x value of ith input coordinate
          
      The value of the resulting polynomial is updated as the sum of itself with the ith term
      multiplied by the corresponding DifferenceLadder.
      
      It return the resulting polynomial.      
      
 1.LegrangesInterpolation:Equally spaced inputs
      It takes in a collection of 2d coordinates of the form (x, y).
    
      It then proceeds to generate the values of (y of the ith input coordinate)/(x of ith input coordinate - x of jth input coordinate), where j is not i.
      Lets call this quantity the ith scale.
      The corresponding polynomial is set to p(x) = 0.
        
      In an interative loop, 
      polynomial containing all xs of jth input coordinate as roots, where j is not i, is generated.
      the ith term is then set to be the product of this polynomial and the ith scale
    
      The value of the resulting polynomial is updated as the sum of itself with the ith term 
      multiplied by the corresponding DifferenceLadder.
      
      It return the resulting polynomial.      
