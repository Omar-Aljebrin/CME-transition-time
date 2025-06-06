Python_fit2h_V1 & Python_fit2h_V2

Python_fit2h_V1 and Python_fit2h_V2 are two different versions of the fit2h script originally developed in IDL. The only difference between them is how the data is imported:

- Version 1 processes data from 2013 and 2011.
- Version 2 extracts data from 2016 and 2022.

To use the code, simply specify the correct file path and filename, and it should run as expected.

Important Notes:
This is still a rough translation of the original IDL version. It is not fully optimized and contains several issues, but it serves as a good starting point.

--------------------------------------------------------

Current Progress:

2013 Data:
- Successfully fitted the exponential-linear function and the power function.
- The other two functions have not yet converged.
- Plotted raw position vs. time and velocity vs. time for comparison.

In the code, I refer to different datasets using the prefix "ht" followed by a number.
- "ht" represents the sample data from 2013.
- "ht4" refers to the data from 2011.
- However, these names are specific to my setup, and you may need to adjust them based on your own file structure.

2011 Data:
- Currently unable to get the functions to fit properly.
- The next steps are:
  1. Accurately plot velocity data.
  2. Attempt to fit various functions, likely in this order:
     - Linear
     - Exponential
     - Quadratic
  3. Once a function successfully fits, calculate the error using the chi-squared formula in the code.
  4. Compute the derivative of the fitted function to analyze velocity trends.

--------------------------------------------------------

Contact:
If you have any questions, feel free to reach out:
Email: oaljebri@gmu.edu