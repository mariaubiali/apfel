[![Circle CI](https://circleci.com/gh/scarrazza/apfel/tree/master.svg?style=svg)](https://circleci.com/gh/scarrazza/apfel/tree/master)

![alt text](https://github.com/scarrazza/apfel/raw/master/resources/logoapfel.png "Logo APFEL")

# APFEL: A PDF Evolution Library

Visit: http://apfel.hepforge.org and http://apfel.mi.infn.it/
 
APFEL is able to perform DGLAP evolution up to NNLO in QCD and to LO
in QED, both with pole and MSbar masses. The coupled DGLAP QCD+QED
evolution equations are solved in x-space by means of higher order
interpolations and Runge-Kutta techniques, allowing to explore
different options for the treatment of subleading terms.

The APFEL library implements the following tools:

- APFEL Web: APFEL provides an online web-app which integrates several
HEP softwares providing a complete suite for PDF analysis
(http://apfel.mi.infn.it/).

- APFEL GUI: APFEL provides the user with a graphical interface which
performs PDF, luminosity and DIS observables plots in real time,
granting an easy access to all the APFEL's functionalities.

## Download

You can obtain APFEL directly from the github repository:

https://github.com/scarrazza/apfel/releases

For the last development version you can clone the master code:

```Shell
git clone https://github.com/scarrazza/apfel.git
```

For the latest tag:

```Shell
git tag -l
git checkout tags/tag_name
```

## Installation 

Checkout the code and compile the code using the
following procedure:

```Shell
cd apfel
./configure --prefix=/where/install/apfel #(optional)
make && make install
```

By the default, if prefix is not set the program is installed in
/usr/local. If you define a custom prefix, remember to export
apfel/lib to the LD_LIBRARY_PATH. APFEL GUI requires ROOT (> 5), qmake
(> 4.5) and APFEL (> 1.0.1). The installation steps are:

```Shell
cd apfel/apfelGUI
qmake
make && make install
apfelgui
```
If you observe problems when running the GUI please set your LC_ALL=C.

## References

- V. Bertone, S. Carrazza, J. Rojo, *Doped Parton Distributions*, [arXiv:1509.04022](http://arxiv.org/abs/1509.04022).
- V. Bertone, S. Carrazza, D. Pagani, M. Zaro, *On the Impact of Lepton PDFs*, [arXiv:1508.07002](http://arxiv.org/abs/1508.07002).
- V. Bertone, S. Carrazza, E. Nocera, *Reference results for time-like evolution up to O(α_s^3)*, [arXiv:1501.00494](http://arxiv.org/abs/1501.00494).
- S. Carrazza, A. Ferrara, D. Palazzo, J. Rojo, *APFEL Web: a web-based application for the graphical visualization of parton distribution functions*, [arXiv:1410.5456](http://arxiv.org/abs/1410.5456).
- V. Bertone, S. Carrazza, J. Rojo, *APFEL: A PDF Evolution Library with QED corrections*, [arXiv:1310.1394](http://arxiv.org/abs/arXiv:1310.1394).

## Contact Information

Maintainers: Valerio Bertone, Stefano Carrazza

Homepage: http://apfel.hepforge.org/
