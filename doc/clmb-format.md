Coulomb file format
===================

This file format indicates the DMTP (Mulliken, ChelpG, DMA, CAMM etc.)
distributions. The format is free. It has three major fields:
  * distributed moments
  * atomic symbols and coordinates
  * origins of distributed moments (optional)

If the latter field is absent it means that origins=structure. Note, that
everything which is beyond these field has no meaning so you can write
as many lines as you wish at the beginning to describe your distribution
of moments precisely (molecule etc.)
The examplary files can be found [here](https://github.com/globulion/clmb/tree/master/doc/examples). 
As you can see, the first field contains:
  * `Zeroth-order property` - charges
  * `First-order property`  - dipoles
  * `Second-order property` - quadrupoles
  * `Third-order property`  - octupoles

At the end of this field always there is an information
about the type of distributed quadrupoles and octupoles
which can be either in **promitive** or
**traceless** form, respectively. To understand the meaning of these
forms see the [documentation](https://github.com/globulion/clmb/blob/master/doc/coulomb.pdf)
of *Coulomb.py* package. If the number if distributed centers is
equal to one (molecular multipole moments) then additionally 
the origin position (in **Angstroms**) is printed just next to the information about form.

### Reading the coulomb format file

If you have installed LIBBBG you can use `ParseDMA` function in Python
using `coulomb` as a specifier of the format:
```
camm = ParseDMA('file.clmb','coulomb')
```
