# Metabokiller
Metabokiller: Artificial Intelligence uncovers carcinogenic human metabolites<br/><br/>

## Workflow 

<img src="Images/GH_Cover.png"> 

## Introduction

Metabokiller is a carcinogen-independent ensemble model for carcinogenicity prediction. It is a novel ensemble classifier that accurately recognizes carcinogens by quantitatively assessing their chemical composition as well as potential to induce proliferation, oxidative stress, genotoxicity, alterations in epigenetic signatures, and activation of anti-apoptotic pathways.<br/><br/>

## How to use Metabokiller?


### Installation using pip 
```
pip install MetaboKiller

```

- Example

To get predictions for individual carcinogenic properties:<br/>
```
from MetaboKiller import mk_predictor as mk
```
Prepare a list of canonical SMILES (Openbabel generated) strings
```
smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Run predictions on any of the carcinogenic property of interest (e.g. epigenetic modifications)
```
mk.Epigenetics(smiles)
```
Save the result as Pandas dataframe
```
result = mk.Epigenetics(smiles)
```

- List of carcinogenic properties available in  **mk** 
```
mk.Epigenetics()
mk.Oxidative()
mk.GInstability()
mk.Electrophile()
mk.Proliferation()
mk.Apoptosis()
```

- To get predictions for all available carcinogenic properties along with their explainability:
```
from MetaboKiller import EnsembleMK
```

Prepare a list of canonical SMILES (Openbabel generated) strings
```
smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Run predictions for all available carcinogenic properties
```
EnsembleMK.predict(smiles)
```
Save the result as Pandas dataframe
```
result = EnsembleMK.predict(smiles)
```
To get explainability of the results on the basis of individual carcinogenic properties (for each SMILES)
```
result,explaination = EnsembleMK.predict(sa,explainability=True)
```


```
# getting output from the explainability object
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt

pdf = PdfPages("Ensmble-Result.pdf")
for fig in explaination:
	fig.savefig(pdf, format='pdf')
pdf.close()
```
<!-- comment -->
