# Carcinogenicity Prediction using Metabokiller
 <br>
<div align="center">
<img src="Images/Asset 1.png"></div>
<br>

## Introduction

Metabokiller offers a novel, machine learning-based approach that accurately recognizes carcinogens by quantitatively assessing their chemical composition as well as potential to induce proliferation, oxidative stress, genomic instability, alterations in epigenetic signatures, and activation of anti-apoptotic pathways, and therefore, obviates the absolute need for bonafide (non)carcinogens for training model. Concomitant with the carcinogenicity prediction, it also reveals the contribution of the aforementioned biochemical processes in carcinogenicity, thereby making the proposed approach highly interpretable. <br/><br/>


The only strong dependency for this resource is [**RDKit**](https://www.rdkit.org/) which can be installed in a local [Conda](https://conda.io/) environment.

```
$ conda create -c conda-forge -n my-rdkit-env rdkit
$ conda activate my-rdkit-env
```

**Major dependencies**
1. [Signaturizer(v1.1.11)](https://gitlabsbnb.irbbarcelona.org/packages/signaturizer)
2. [LIME](https://github.com/marcotcr/lime)

The installation procedure takes less than 5 minutes.
```
$ pip install signaturizer
$ pip install lime
```

**Minor dependencies**
1. os
2. [scikit-learn v1.0.2](https://scikit-learn.org/stable/whats_new/v1.0.html)
3. [pandas](https://pandas.pydata.org/)
4. [numpy](https://numpy.org)
5. [tqdm](https://tqdm.github.io)
6. [joblib](https://pypi.org/project/joblib/)
7. [matplotlib](https://pypi.org/project/matplotlib/)
8. io 
9. [importlib](https://pypi.org/project/importlib/)


## How to use Metabokiller?


### Installation using pip 
```
$ pip install Metabokiller
```

#### Examples

To get predictions for individual carcinogenic properties:<br/>
```
>>> from Metabokiller import mk_predictor as mk
```
Prepare a list of canonical SMILES (Openbabel generated) strings
```
>>> smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Run predictions on any of the carcinogenic property of interest (e.g. epigenetic modifications)
```
>>> mk.Epigenetics(smiles)
```
Save the result as Pandas dataframe
```
result = mk.Epigenetics(smiles)
```

##### Metabokiller supported carcinogen-specific biochemical properties:

1. Epigenetic Alterations 
```
>>> mk.Epigenetics()
```

2. Oxidative stress 
```
>>> mk.Oxidative()
```

3. Electrophilic Property 
```
>>> mk.Electrophile()
```

4. Genomic Instability 
```
>>> mk.GInstability()
```

5. Pro-proliferative response 
```
>>> mk.Proliferation()
```

6. Anti-apoptotic response 
```
>>> mk.Apoptosis()
```


##### To get predictions for all available carcinogenic properties along with their explainability:
```
>>> from Metabokiller import EnsembleMK
```

Prepare a list of canonical SMILES (Openbabel generated) strings
```
>>> smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Run predictions for all available carcinogenic properties
```
>>> EnsembleMK.predict(smiles)
```
Save the result as Pandas dataframe
```
>>> result = EnsembleMK.predict(smiles)
```

##### LIME
	The biochemical property-focused Metabokiller, by the virtue of its construction, offers interpretability by implementing Local Interpretable Model-agnostic Explanations (LIME). An algorithm that provides interpretability with respect to carcinogen-specific biochemical properties for each SMILES provided.


##### To activate interpretability using LIME:

```
>>> result,explanation = EnsembleMK.predict(['ClCC=C', 'C=CCOC(=O)CC(C)C'],explainability=True)
```


```
# getting output from the explainability object
>>> from matplotlib.backends.backend_pdf import PdfPages
>>> from matplotlib import pyplot as plt

>>> pdf = PdfPages("Ensmble-Result.pdf")
>>> for fig in explanation:
...	fig.savefig(pdf, format='pdf')
>>> pdf.close()
```
<!-- comment -->
