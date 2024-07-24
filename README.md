Energy of Quantum Dots (QDs) - Hamiltonian Eigenvalue Calculation

This repository contains Python code for calculating the eigenvalues of a Hamiltonian matrix for quantum dots (QDs). The code constructs an NxN Hamiltonian matrix based on the excitations provided in CSV files and diagonalizes the matrix to obtain the eigenvalues. The results are saved in a .txt file based on the selected material and volume.

Contents
EigenvaluesCalculator.py: The main script for creating and diagonalizing the Hamiltonian matrix.
tablas1.csv, tablas2.csv, tablas3.csv, tablas4.csv, tablas5.csv: CSV files containing the possible excitations for holes and electrons.
README.md: This file.
Getting Started
Prerequisites
Ensure you have Python 3.x installed on your system. The required Python packages are listed below:

math
numpy
scipy
csv
time
multiprocessing
You can install the required packages using pip:

bash

  pip install numpy scipy

Usage

Clone the repository:

bash

git clone https://github.com/CarlosABQ34/Energy-of-Qds.git
cd Energy-of-Qds
Run the EigenvaluesCalculator.py script:

bash

python EigenvaluesCalculator.py

Explanation of the Code

The script reads the excitation data from the tablas*.csv files. Each file contains possible excitations for holes and electrons.
An NxN Hamiltonian matrix is constructed, where N is the number of elements in the CSV files.
The script allows for up to 4 levels of excitation for each dimension (tridimensional).
The Hamiltonian matrix is diagonalized to find its eigenvalues.
The results are saved in a .txt file, labeled by the selected material and specific volume.
Customizing
To reproduce the results or run the code with different parameters:

Modify the EigenvaluesCalculator.py script to select the desired material and volume.

Example
An example usage can be included here once the specific details are provided.

Contributing
We welcome contributions to this project. If you have suggestions for improvements or new features, please open an issue or submit a pull request.

License
This project is licensed under the MIT License. See the LICENSE file for details.

Contact
For any questions or further information, please contact:

Carlosbohorquezq@gmail.com
GitHub: CarlosABQ34
