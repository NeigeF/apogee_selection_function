# Simplified APOGEE selection function
We provide a pre-computed table for the APOGEE selection function (including distance dependent selection due to extinction) in the form of a fits file.

This table was initially computed by Frankel, Sanders, Rix, Ting & Ness (2019) to characterise the inside-out growth of the Milky Way disk, using only APOGEE-1 disk fields. Since this, we have extended the calculation to some fields from APOGEE-2.

It provides a simple way to deal with selection effects for most APOGEE disk fields (but not all). We provide details and examples below and in a jupyter notebook.

## Data
The file is available in this repository at  <br />
https://github.com/NeigeF/apogee_selection_function/blob/master/sf_apogee_dr14_disk.fits

## Tutorial

#### Content
An example tutorial is available as a jupyter notebook in this repository at  <br />
https://github.com/NeigeF/apogee_selection_function/blob/master/selection_function_ap1-2_publ.ipynb.

The example covers
- how to load the file,
- how to use an APOGEE dataset down with the available selection function,
- how to include the selection function in a density model,
- how to fit this model to APOGEE red clump data,
- how to sample mock stars from a density model in APOGEE fields (relatively efficiently),
- details of the literature this work is based on.

If the notebook does not render here, please visit:  <br />
https://nbviewer.jupyter.org/github/NeigeF/apogee_selection_function/blob/master/selection_function_ap1-2_publ.ipynb

#### Authors
Neige Frankel & Jason Sanders.

Please cite Frankel, Sanders, Rix, Ting & Ness (2019) if you find this code useful in your research.


