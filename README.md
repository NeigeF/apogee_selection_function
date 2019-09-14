# apogee_selection_function
fits file for apogee selection function (including distance dependent selection due to extinction)


Initially computed for and used in Frankel+19 'The inside-out growth of the Milky Way disk' only with apogee-1 disk fields. Expanded below to some fields from apogee-2.

Provides a simple way and examples to deal with selection effects with most APOGEE disk fields (but not all)

## Data
The file is available in this repository at
https://github.com/NeigeF/apogee_selection_function/blob/master/sf_apogee_dr14_disk.fits

## Tutorial
An example tutorial showing:
- how to load the file
- how to bring an apogee dataset down to the available selection function
- how to include the selection function in a density model
- and fit this model to apogee data
- how to sample stars in apogee fields (sort of efficiently)

is available as a jupyter notebook in this repository at
https://github.com/NeigeF/apogee_selection_function/blob/master/selection_function_ap1-2_publ.ipynb

if the notebook does not render here, please visit:
https://nbviewer.jupyter.org/github/NeigeF/apogee_selection_function/blob/master/selection_function_ap1-2_publ.ipynb





