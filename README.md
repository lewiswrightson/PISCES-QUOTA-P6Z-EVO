# PISCES QUOTA P6Z EVO

PISCES QUOTA P6Z EVO builds upon PISCES QUOTA P6Z (Wrightson et al. 2022) to include thermal adaptation of _Crocosphaera_ for the standard P5Z PISCES QUOTA framework in NEMO4.2.3.

This model implements thermal performance curves for _Crocosphaera_ growth from long term laboratory cultures investigating thermal adaptation (Duan et al. _in prep_).

This version implements a simple evolutionary clock approach which enables thermal adaptation of the diazotroph phytoplankton functional type (PFT).

Two clocks are implemented that enables adaptation to both warmer (<32<sup>o</sup>C) and colder (<32<sup>o</sup>C) temperatures. 

The first clock counts up when the temperature is above 32<sup>o</sup>C and switches the diazotroph from a 28<sup>o</sup>C adapted strain to a 32<sup>o</sup>C adapted strain once the diazotroph has been exposed to temperatures warmer than 32oC for the prescribed adaptation time scale.

Following the thermal adaptation of the diazotroph a second clock activates thatcounts up when the temperature is below 32<sup>o</sup>C and switches the diazotroph back to a 28<sup>o</sup>C adapted strain when the temperature has remained lower that 32<sup>o</sup>C for the prescribed adaptation time scale.

The clocks are tracers which are advected in the model to track the temperatures the diazotrophs are exposed to.

Wrightson, L., Yang, N., Mahaffey, C., Hutchins, D. A., & Tagliabue, A. (2022). Integrating the impact of global change on the niche and physiology of marine nitrogen-fixing cyanobacteria. Global Change Biology, 28(23), 7078-7093. https://doi.org/10.1111/gcb.16399 

### Prerequisites / Getting Started:

MY_SRC (Model subroutines)

EXPREF directory (namelists, file_def, field_def and other xml files)

This version of this model requires flags to be set in **namelist_pisces_ref.bait.p6z.evo** to enable thermal adaptation:

Flag **ln_p6z** to select explicit diazotrophy

Flag **ln_tricho** = false (_Crocosphaera_ selected)

Flag **ln_evolve** activates thermal adaptation

Flag **ln_tiue** = false (no temperature dependence of the nitrogen fixation Fe use efficiency)

Flag **ln_tpue** - false (no temperature dependence of the nitrogen fixation P use efficiency)

namelist parameter **evotime** sets the adaptation timescale: 14 = 2 weeks (short) and 730 = 2 years (long)

### Installing: 

NEMO 4.2.3:
git clone --branch 4.2.3 https://forge.nemo-ocean.eu/nemo/nemo.git nemo_4.2.3

The configuration is compiled using ./makenemo -r P6Z_EVO

Model restart file: https://zenodo.org/records/20395315

### Authors:

Lewis Wrightson and Alessandro Tagliabue 

### Acknowledgements:

This work was supported by the Natural Environment Research Council grant NE/X014908/1 awarded to Alessandro Tagliabue and U.S. National Science Foundation grants OCE 2149837 and OCE 2336534 to David A. Hutchins and Fei-Xue Fu.
