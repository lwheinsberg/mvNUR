mvBIMBAM Installation Hack
================
Lacey W. Heinsberg



# Copyright information

Copyright 2024, University of Pittsburgh. All Rights Reserved. License:
[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# mvBIMBAM

## Mac

The motivation behind modifying the makefile was that during
installation (of mvBIMBAM “first release”), test users were getting the
following error:

    cd src && /Applications/Xcode.app/Contents/Developer/usr/bin/make
    g++ -static-libgcc -DIMPUTATION   -O3 control.o fpmath.o indiv.o diploid.o haploid.o model.o param.o  fp.o -lm libgsl.a libgslcblas.a  -o bimbam
    clang: error: unsupported option '-static-libgcc'
    make[1]: *** [fp] Error 1
    make: *** [all] Error 2

which suggests there is an issue with the build process. To get around
this, we “hacked” a solution by modifying the Makefile (located within
the `src` folder). See a summary below, or a “track changes” comparison
in `makefile_compare_AppleIntel.docx` or
`makefile_compare_AppleSilicon_M1.docx` (located here in the
`InstallationAppendix` folder) for details on how the makefile was
modified.

Note that the makefile changes are different for an Apple Silicon (e.g.,
2021 M1) vs. an Apple Intel. If you are not sure which Mac you have, you
can click on the apple symbol in the upper left corner of your machine.
If you have an Apple silicon, the Chip will be listed as Apple M1 or M2.
If you have the Apple Intel, you will see the processor listed as Intel
Core i5, i7, or similar.

Note that this modified makefile and redistribution of mvBIMBAM software
is redistributed as permitted under the creator’s original GPL-3
license.

### APPLE INTEL INSTRUCTIONS

**CHANGE A** Add GSL Compiler and Linker Flags

    CFLAGS += -I/usr/local/Cellar/gsl/2.7.1/include
    LDFLAGS += -L/usr/local/Cellar/gsl/2.7.1/lib

These flags include the GSL headers and dynamically specify the location
of the GSL libraries on a machine.

**CHANGE B** Update LIBS Variable. Specifically, modify:

    LIBS += -lm libgsl.a libgslcblas.a

to:

    LIBS += -lm -lgsl -lgslcblas

This directly links the GSL libraries for math, GSL core, and GSL CBLAS.

**CHANGE C** Update Compilation Command: Remove the `-static-libgcc`
flag. Specifically, modify this line:

    fp: fp.o $(OBJS); $(CC) -static-libgcc $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam

to

    fp: fp.o $(OBJS); $(CC) $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam

The `-static-libgcc` flag is related to linking the GCC (GNU Compiler
Collection) runtime libraries statically. However, this flag is not
supported by the Clang compiler, which is commonly used on macOS
(including Apple Silicon M1 Macs). By removing it, the build rule
becomes compatible with Clang and more common macOS build setups.

A copy of the original Makefile, modified Makefile, and a Word document
comparing the two can be found in the `Install_Problems` folder.

Note that this problem has been opened as an issue on the mvBIMBAM
GitHub page <https://github.com/heejungshim/mvBIMBAM/issues/3>. Please
check there for updates from the program architects.

### APPLE SILICON (M1)

**CHANGE A** Add GSL Compiler and Linker Flags

    CFLAGS += `gsl-config --cflags`
    LIBS += `gsl-config --libs`

These flags include the GSL headers and dynamically specify the location
of the GSL libraries on your machine.

**CHANGE B** Update LIBS Variable: Specifically, modify:

    LIBS += -lm libgsl.a libgslcblas.a

to:

    LIBS += -lm -lgsl -lgslcblas

This directly links the GSL libraries for math, GSL core, and GSL CBLAS.

**CHANGE C** Update Compilation Command: The compilation command for the
“fp” target should be updated to remove the `-static-libgcc` flag and
add a linker. Specifically, modify:

    fp: fp.o $(OBJS); $(CC) -static-libgcc $(CFLAGS) $(OBJS) fp.o $(LIBS) -o bimbam

to

    fp: fp.o $(OBJS); $(CC) $(CFLAGS) $(OBJS) fp.o $(LDFLAGS) $(LIBS) -o bimbam

The `-static-libgcc` flag is related to linking the GCC (GNU Compiler
Collection) runtime libraries statically. However, this flag is not
supported by the Clang compiler, which is commonly used on macOS
(including Apple Silicon M1 Macs). By removing it, the build rule
becomes compatible with Clang and more common macOS build setups.

Note also that Apple Silicon also requires the extra `$(LDFLAGS)` (where
Intel does not) which is used to specify linker-related flags and
options that may be specific to the architecture of M1 (not needed by
Intel).

A copy of the original Makefile, modified Makefile, and a Word document
comparing the two can be found in the `Install_Problems` folder.

Note that this problem has been opened as an issue on the mvBIMBAM
GitHub page <https://github.com/heejungshim/mvBIMBAM/issues/3>. Please
check there for updates from the program architects.

# Session information

``` r
sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.2.1 (2022-06-23)
    ##  os       macOS 14.3.1
    ##  system   aarch64, darwin20
    ##  ui       X11
    ##  language (EN)
    ##  collate  en_US.UTF-8
    ##  ctype    en_US.UTF-8
    ##  tz       America/Toronto
    ##  date     2024-07-22
    ##  pandoc   3.1.1 @ /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/ (via rmarkdown)
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date (UTC) lib source
    ##  cli           3.6.2   2023-12-11 [1] CRAN (R 4.2.3)
    ##  digest        0.6.35  2024-03-11 [1] CRAN (R 4.2.3)
    ##  evaluate      0.23    2023-11-01 [1] CRAN (R 4.2.0)
    ##  fastmap       1.2.0   2024-05-15 [1] CRAN (R 4.2.3)
    ##  htmltools     0.5.8.1 2024-04-04 [1] CRAN (R 4.2.3)
    ##  knitr       * 1.46    2024-04-06 [1] CRAN (R 4.2.1)
    ##  rlang         1.1.3   2024-01-10 [1] CRAN (R 4.2.3)
    ##  rmarkdown     2.27    2024-05-17 [1] CRAN (R 4.2.3)
    ##  rstudioapi    0.16.0  2024-03-24 [1] CRAN (R 4.2.3)
    ##  sessioninfo   1.2.2   2021-12-06 [1] CRAN (R 4.2.0)
    ##  xfun          0.44    2024-05-15 [1] CRAN (R 4.2.3)
    ##  yaml          2.3.8   2023-12-11 [1] CRAN (R 4.2.3)
    ## 
    ##  [1] /Users/law145/Library/R/arm64/4.2/library
    ##  [2] /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
