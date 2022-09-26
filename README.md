Paper data

Data and code to show screening process and results

ML_screening.ipynb - jupyter notebook with our code. The notebook has been used to demonstrate the process of mass screening NCAA by the model in the article. The notebook has been run in 5 steps from start to finish, only the second step needs to be run with another tool or script, the other steps can be run here. It has been tested and can be run successfully.

molslibrary - A folder containing files before and after screening of small molecule libraries and csv files after integration of small molecule libraries.

model - A folder containing the building BT model covered in the article, used here for model prediction. Also included are amino acid descriptor file, data file and small molecule smiles format file for model training.

data - A folder containing the chemopy descriptor data file for the intergrated small molecule library, and the csv file with the protein features added for model prediction. Chemopy descriptors was calculated by PyBioMed(https://github.com/gadsbyfly/PyBioMed).

pybiomed.yaml - conda environment for calculating chemopy descriptors.


System Requirements
The notebook is runs in Windows.

Hardware requirements
Only a standard computer is required. 

Software requirements
Anaconda can be download by https://www.anaconda.com/products/individual

Related Websites
ChemDes[1] (http://www.scbdd.com/chemopy_desc/index/) 
PyBioMed[2] (https://github.com/gadsbyfly/PyBioMed)

[1]. Jie Dong, Dong-Sheng Cao, Hong-Yu Miao, Shao Liu, Bai-Chuan Deng, Yong-Huan Yun, Ning-Ning Wang, Ai-Ping Lu, Wen-Bin Zeng, Alex Chen. ChemDes: an integrated web-based platform for molecular descriptor and fingerprint computation. Journal of Cheminformatics 2015, 7:60
[2]. Dong J, Yao Z J, Zhang L, et al. PyBioMed: a python library for various molecular representations of chemicals, proteins and DNAs and their interactions[J]. Journal of cheminformatics, 2018, 10(1): 16.

Package  requirements
If you expect to run steps other than step 2, simply install the following packages:
Python 3.8.13
numpy 1.22.3
pandas 1.4.3 
sklearn 0.24.1 
rdkit 2022.03.4 
joblib 1.1.0 

If you expect to calculate chemopy molecular descriptors using your own computer, you could create a conda environment by 'pybiomae.yaml'
```
conda env create -f pybiomed
```
 

License
this project is covered under the MIT licence