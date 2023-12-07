# Jormungandr
A program to classify nondiagonal triples on nice diagrams to construct indefinite Einstein solvmanifolds that are not of pseudo-Iwasawa type. Basd on the following paper:

Diego Conti, Federico A. Rossi and Romeo Segnan Dalmasso, A construction of Einstein solvmanifolds not based on nilsolitons [arxiv.org/abs/2312.03125](https://arxiv.org/abs/2312.03125)

The tables are

## Requirements

*Jormungandr* requires [cmake](https://cmake.org/) and [Wedge](https://github.com/diego-conti/wedge)

To compile, run

	mkdir build
	cd build
	cmake ..
	cmake --build .

You may need to set the WEDGE_PATH environment variable to point to your installation of Wedge, e.g.

	export WEDGE_PATH=/home/user/wedge

## Usage

The tables contained in the paper are obtained by invoking the program as 

    ./jormungandr --dimension n --surjective --free-parameters "=0"

If ```--surjective``` is not indicated, *Jormungandr* does not exclude cases where the root matrix is nonsurjective, and does not attempt to compute metrics.

Other command line options are available; run ```jormungandr ``` without parameters for a list.

The resulting output is in the subdirectory ```  tables ``` 
