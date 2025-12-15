# Multiome analysis
This repository used MOFA2 and quadbio to analyze multiome data: scRNA-seq + scATAC-seq from disease sample compare with control sample to find the transcription factor that cause the disease. You can also try SCENIC+.

How to install SCENIC+  
Python version to avoid error: <=3.11.8,>=3.8
```
git clone https://github.com/aertslab/scenicplus
module load miniconda
conda create --name scenicplus python=3.11.8
conda activate scenicplus
git checkout development
pip install .
```
