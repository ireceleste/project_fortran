# project_fortran
Project for Scientific Programming II Course

## Description

This program implements a maximum-likelihood fitting procedure with a gradient descent minimisation algorithm. It is possible to fit data provided by the user or to data generated accordingly to the chosen fit model. 

## Preparation

Compile the code using the provided Makefile:

```
make
```

## Usage
To run the gradient descent fitter, having set all the input parameters in the file card.dat, run the following command:
```
bin/main card.dat
```
The usage instructions and the input card format requirements can be found by running:
```
bin/main --help
```
or 
```
bin/main -h
```
At the end of the help mode, the program can optionally provide a template for the input card. 
More details on the available fit models and on the parameters can be found by running:
```
bin/main --model
```
or 
```
bin/main -m
``` 
Once the fit has been performed, if the results have been saved to an output file, it is possible to plot the data with the best-fit model overlaid and the pulls:
```
plot.py <output_file>
```