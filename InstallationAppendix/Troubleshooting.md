mvBIMBAM Installation Troubleshooting (Check for GSL)
================
Lacey W. Heinsberg



# Copyright information

Copyright 2024, University of Pittsburgh. All Rights Reserved. License:
[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# Troubleshooting mvBIMBAM installation

Note that you need to have GSL (GNU Scientific Library) installed on
your machine before you can install mvBIMBAM. GSL is an open-source
software library that provides a wide range of mathematical and
statistical functions for scientific and numerical computing. If you are
having trouble installing mvBIMBAM, lack of GSL seems to be the most
frequent cause.

Using the terminal, check to see if you have GSL installed on your
machine by asking for the version number:

    gsl-config --version

For example, on my machine, this command returns “2.7.1”.

If GSL is not installed, the simplest way to install it is via homebrew.
Homebrew is a popular package manager for macOS and Linux. It provides a
convenient way to install, manage, and update software packages and
libraries on your computer. Homebrew simplifies the process of setting
up development tools and applications.

First, install homebrew using the following terminal command:

    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"

Confirm that it has been installed, and check for the path:

    brew --version
    brew --prefix

For example, on my machine, the first command returns “Homebrew 4.0.23”
and the second returns “/opt/homebrew”. Note that this tutorial was
developed on an Apple silicon (M1) machine. If you have an Apple Intel,
the default homebrew install path will be “usr/local/Cellar” or
“usr/local”.

Now let’s move on to installing GSL:

    brew install gsl

Confirm install and path:

    gsl-config --version 
    gsl-config --prefix

Voila! (Hopefully :))

Please move back to the installation instructions in the README file and
try again!

# Contact information

If you have any questions or comments, please feel free to contact me!

Lacey W. Heinsberg, PhD, RN: <law145@pitt.edu>
