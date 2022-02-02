# RUScal
Software to analyze resonant ultrasound spectroscopy measurements.

********** READ ME *********************************
Release 2022
Last update: 12/8/2021

RUScal is open-source research software for numerical calculations of the eigenvalues of a free-body mechanical resonator
based on Rayleigh-Ritz approximation. This software was developed in part at Oak Ridge National Laboratory under 
the Laboratory Directed Research and Development Program of Oak Ridge National Laboratory, managed by UT-Battelle, 
LLC, for the U. S. Department of Energy.

For publications or presentations that made use of this software, please cite the following:
[insert when accepted]

********** Contacts ********************************
Author: Raphael Hermann, hermannrp@ornl.gov
Additional contact: James Torres, torresdrex@gmail.com 

********** Getting Started *************************
There is no auto-installer included, so the software is compiled manually. 
The source code may be downloaded directly from the GitHub database.
URL: https://github.com/JT-ORNL/RUScal.git

Also included in the repository is a GUI that runs using the latest version of Python 3.x.
Additional Python libraries required:
openpyxl
numpy
pandas
pillow

---------- Windows ---------------------------------
Download ruscal.exe and move to the desired directory. To execute this program system-wide, add the path of the
directory containing ruscal.exe to the Path of System Variables.

---------- Mac and Linux ---------------------------
Ensure the following files are within the same directory:
rus.c
matrix.c
matrix.h
rusin.dat
MakeFile

Compile code and create the application by running MakeFile.

********** Instructions ****************************
The default input file is named rusin.dat, one of which has been provided to modify with appropriate parameters.

---------- Using the GUI ----------------------------
1. Enter all sample parameters and frequencies OR click "Read input file." The latter will read the rusio in the
   active directory.
2. Click "Write rusin.dat" to create or overwrite the rusin file in the active directory.
3. Click "Start ruscal" to execute RUScal code.

Use the follow procedure to make a new rusin file:

1. Enter sample/experiment information in the Sample Name entry box. Use a ':' to add a tracking parameter, eg.,
   temperature - place the value after the colon.
2. Choose a sample shape.
3. Choose a Bravais class. This will correspond to the number of independent elastic constants left activated on 
   the right side of the screen.
4. Choose a polynomial order. Higher polynomial leads to increased accuracy but longer computation time. 10 to 12 is
   recommended.
5. Choose a calculation mode:
   5a. 1 = Calculates the integer number of resonant frequencies entered in the box below the radio button.
   5b. 0 = Triggers the code to fit the free parameters that minimize the difference between measured and calculated
           resonant frequencies.
   5c. -1 = Monte Carlo search for elastic constants. Change the values of Bounds for each elastic constant to the
            desired maximum deviation in percent of the initial value. E.g., for 20% deviation, set Bounds to 20.
   5d. -3 = Calculates the chi-squared error of each elastic constant within +-0.05% of given value in steps of 0.01%.
   5e. -2 = Calculates resonant frequencies as a function of dimension or Euler angle systematically. From the adjacent 
            dropdown menu, select 1 dimension or Euler angle to vary. The number of calculated frequencies is equal to 
            the number of measured frequencies given in the Frequencies table at the bottom of the screen. 
   5f. -4 = Calculates chi-squared error for each elastic constant similar to Mode=-3 but with user-defined bounds.
            Set the Bounds value of each elastic constant to a number equal to 1/1000 deviation. E.g., entering a 
            value of 40 is equivalent to varying the elastic constant by 40/1000 or 4%. Only three constants can be
            varied at a time. Entering nonzero Bounds for more than 3 will only calculate for the first three.
6. Toggle the calculation of derivatives, df/dcxx. Checked = calculate.
7. Select number of mirror planes. Use default of zero if unsure.
8. Enter values for each Euler angle, if applicable (for single crystals).
9. Check boxes for each Euler angle: checked = fixed, unchecked = allow to be fit.
10. Enter dimensions in centimeters. Selection of a shape will update the comment under each entry box. Potato shape 
    unlocks separate entry boxes.
11. Check box for dimensions: checked = fix dimensions, unchecked = allow all to fit (volume is kept constant).
12. Enter sample mass in grams. Re-selecting the shape again will calculate the density.
13. Load resonant frequencies. 
    13a. Choose the file type to upload.
    13b. Click "Choose File," which opens the current directory. Select file containing 3 column list of data:
         column 1 = measured frequencies in MHz, col. 2 = calculated frequencies in MHz, col. 3 = weights between 0 and 1.
         If any measured and/or calculated frequencies are unknown, then enter zeros in the file.
    13c. Click "Load File." Bottom-left table will populate with data. 
    13d. Edit a row of data: 1) click on row, 2) click "Select," 3) edit values, and 4) click "Save/Update."
14. Enter elastic constants at the top right. Units are in Mbar (1 Mbar = 100 GPa).
15. Check boxes for each elastic constant: checked = fixed, unchecked = allow to fit. 
    Bounds are used only for the Monte Carlo (-1) and Grid Search (-4) modes.
16. Click "Write rusin.dat" to finish. The rusin file will be created or overwritten in the active directory.


