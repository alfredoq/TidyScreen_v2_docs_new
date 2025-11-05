---
title: 3. Receptors modeling and validation
---

**Step 1**: Retrieve receptors crystallographic structures

Using bash after project creation

```python
# Retrieve a structure corresponding to CZP
$ mkdir -p /PATH/TO/PROJECT/docking/raw_data/2OZ2
$ wget -O /PATH/TO/PROJECT/docking/raw_data/2OZ2/2OZ2.pdb https://files.rcsb.org/download/2OZ2.pdb
```

Alternatively you can download 2OZ2.pdb from [here](/downloads/LA_Workshop/2OZ2.pdb) (*right-click and download*)


**Step 2**: Processing of raw crystallographic structures

Processing of CZP crystal structure pdb code: 2OZ2

```python 
>>> from tidyscreen import tidyscreen as ts
>>> from tidyscreen.moldock import moldock as md

### Activate the project
>>> la_workshop = ts.ActivateProject("la_workshop_2025")

### Instantiate a MolDock object to process the receptor models
>>> la_workshop_moldock = md.MolDock(la_workshop)

# Process CZP crystal structure
structure_file = "/PATH/TO/PROJECT/docking/raw_data/2OZ2/2OZ2.pdb"
>>> la_workshop_moldock.process_raw_pdb(structure_file, x_coord=6.7,  y_coord=-3.7, z_coord=-22.0, x_points=40, y_points=40, z_points=40) # Note that coordinates and size of the grid box are provided


# Outputs
Please provide a brief description to store with the receptor model: # Provide the requested information

# Outputs
Chains found in the pdb file: {'C', 'A'} # Two chains were detected in the crystal
Chains detected: {'C', 'A'} # ID of the detected chains
Provide the chain identifier to mantain: ('all' to keep all chains) # Indicate the chain to retain. In this tutorial we will use chain 'A'

# Outputs
Non-standard residues detected:
The following non-standard residues were found in the pdb file: ['D1R', 'SO4']. 
Do you want to mantain ONE of them as a REFERENCE pdb file? (y/n): # answer 'y' if a reference ligand is to be kept in a separate file. 

# Outputs

Type the 3-letter code of the residue you want to save as REFERENCE FILE:  # provide 'D1R' to retain K777 as reference ligands

# Outputs
The following non-standard residues were found in the pdb file: ['D1R', 'SO4']. 
Do you want to mantain ONE of them in the processed receptor pdb file? (y/n): # Answer 'n'(as many times as non-standard residues exists)

# Outputs
# Informs the creation of:
# /PATH/TO/PROJECT/docking/raw_data/2OZ2/receptor.mol2
# /PATH/TO/PROJECT/docking/raw_data/2OZ2/receptor.pdbqt
```

**Step 3**: Compute Autodock4 grid files

```bash
$ cd /PATH/TO/PROJECT/docking/raw_data/2OZ2/
$ autogrid4 -p receptor.gpf -l receptor.glg

# The corresponding grid maps will be created in the receptor folder
```

You should end up with a receptor folder a included in asdsa [this](/downloads/LA_Workshop/2OZ2.zip) folder (*right-click to download compressed folder*)


**Step 4**: Import the prepared receptor into TidyScreen

```python 
>>> la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/2OZ2/")

# Outputs
Provide a brief description of the receptor model: # Input the required info

# Outputs
Succesfully stored receptor model
Successfully stored the receptor model in the database.
```

**Step 5**: For the purposes of this tutorial, we will be using a diverse set of receptor models:

* Refined CZP structure [donwload zipped folder](/downloads/LA_Workshop/CZP_refined.zip)
* hCatL X-ray structure [donwload zipped folder](/downloads/LA_Workshop/8hfv.zip)
* Refined hCatL structure [donwload zipped folder](/downloads/LA_Workshop/hCatL_refined.zip)

After downloading and extracting the folders, these alternative receptor models can be imported:

```python 
>>> la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/8hfv")
...
>>> la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/CZP_refined")
...
>>> la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/hCatLs_refined")
```


Activate ChemSpace to input ligands

```python
from tidyscreen.chemspace import chemspace as cs
la_workshop_cs = cs.ChemSpace(la_workshop)
```

Create K777 from the SMILES
```python 
csv_file = "/PATH/TO/PROJECT/chemspace/raw_data/K777.csv"
la_workshop_cs.input_csv(csv_file)
```

Generate mols for K777
```python
la_workshop_cs.generate_mols_in_table("K777", conf_rank=25, timeout=300, pdbqt_method="meeko", charge_method="gas")
```

Start creating docking conditions, Create a custom set of docking parameters
```python 
la_workshop_moldock.create_docking_params_set()
```

Extra: Use of sqlite_web to create a second set of docking parameters with n_runs = 100

Create a docking assay for K777 against the CZP crystal structure receptor
```python 
la_workshop_moldock.dock_table("K777",id_receptor_model=1,id_docking_params=2) # docking assay 1
```

Activate Docking Analysis dimension
```python
from tidyscreen.docking_analysis import docking_analysis as docking_analysis
la_workshop_docking_analysis = docking_analysis.DockingAnalysis(la_workshop)
```

Analyze the docking results for assay 1 (docking of K777 against CZP crystal receptor)
```python
la_workshop_docking_analysis.process_docking_assay(assay_id=1,max_poses=10, extract_poses=1) # Notice that we extract only 10 poses per molecule and save them a pdb file
```

Using SQL queries, inspect the results of the K777 docking against CZP RX raw receptor
```bash
sqlite3 -header -table docking/docking_assays/assay_1/assay_1.db "SELECT LigName, sub_pose, cluster_size, docking_score FROM Results"'

################
### Outputs ####
################
# +-----------------------------+--------------------------------+--------------+---------------+
# |           LigName           |            sub_pose            | cluster_size | docking_score |
# +-----------------------------+--------------------------------+--------------+---------------+
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 10           | -4.94         | # This Docked pose does not reproduce RX
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 5            | -2.24         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 3            | -1.64         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 1            | -1.54         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 3            | -1.38         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 2            | -1.33         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 7            | -1.06         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 1            | -0.86         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 6            | -0.8          |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 3            | -0.76         |
# +-----------------------------+--------------------------------+--------------+---------------+
```

Create a docking assay for K777 against the CZP refined structure receptor
```python
la_workshop_moldock.dock_table("K777",id_receptor_model=3,id_docking_params=2) # docking assay 2
```

Analyze the docking results for assay 2 (docking of K777 against CZP refined receptor)
```
la_workshop_docking_analysis.process_docking_assay(assay_id=2,max_poses=10, extract_poses=1) # Notice that we extract only 10 poses per molecule and save them a pdb file
```

Using SQL queries, inspect the results of the K777 docking against CZP RX raw receptor
```bash
sqlite3 -header -table docking/docking_assays/assay_2/assay_2.db "SELECT LigName, sub_pose, cluster_size, docking_score FROM Results"

################
### Outputs ####
################
# +-----------------------------+--------------------------------+--------------+---------------+
# |           LigName           |            sub_pose            | cluster_size | docking_score |
# +-----------------------------+--------------------------------+--------------+---------------+
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 29           | -5.42         | # This Docked pose reproduces the reference
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 7            | -5.34         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 2            | -4.64         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 1            | -4.61         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 5            | -4.61         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 2            | -4.55         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 7            | -4.34         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 1            | -4.31         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 4            | -4.3          |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 2            | -4.28         |
# +-----------------------------+--------------------------------+--------------+---------------+
```

Create a docking assay for K777 against the hCatL RX Raw receptor
```python
la_workshop_moldock.dock_table("K777",id_receptor_model=2,id_docking_params=2) # docking assay 3
```

Analyze the docking results for assay 3 (docking of K777 against CZP refined receptor)
```python
la_workshop_docking_analysis.process_docking_assay(assay_id=3,max_poses=10, extract_poses=1) # Notice that we extract only 10 poses per molecule and save them a pdb file
```

Using SQL queries, inspect the results of the K777 docking against CZP RX raw receptor
```bash
sqlite3 -header -table docking/docking_assays/assay_3/assay_3.db "SELECT LigName, sub_pose, cluster_size, docking_score FROM Results"

################
### Outputs ####
################
# +-----------------------------+--------------------------------+--------------+---------------+
# |           LigName           |            sub_pose            | cluster_size | docking_score |
# +-----------------------------+--------------------------------+--------------+---------------+
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 6            | -4.08         | # This Docked pose does not reproduce RX
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 7            | -3.93         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 3            | -2.97         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 4            | -2.62         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 5            | -2.0          |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 8            | -1.98         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 4            | -1.52         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 1            | -1.33         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 4            | -1.25         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 3            | -0.75         |
# +-----------------------------+--------------------------------+--------------+---------------+

```

Create a docking assay for K777 against the hCatL refined structure receptor
```python
la_workshop_moldock.dock_table("K777",id_receptor_model=4,id_docking_params=2) # docking assay 4
```

Analyze the docking results for assay 4 (docking of K777 against hCatL refined receptor)
```python 
la_workshop_docking_analysis.process_docking_assay(assay_id=4,max_poses=10, extract_poses=1) # Notice that we extract only 10 poses per molecule and save them a pdb file
```

Using SQL queries, inspect the results of the K777 docking against hCatL refined receptor
```bash
sqlite3 -header -table docking/docking_assays/assay_4/assay_4.db "SELECT LigName, sub_pose, cluster_size, docking_score FROM Results"

################
### Outputs ####
################
# -----------------------------+--------------------------------+--------------+---------------+
# |           LigName           |            sub_pose            | cluster_size | docking_score |
# +-----------------------------+--------------------------------+--------------+---------------+
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 22           | -5.17         | # This Docked pose reproduces the reference
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 7            | -5.01         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 1            | -4.92         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 4            | -4.67         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 7            | -4.58         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 8            | -4.32         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 1            | -4.25         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 2            | -4.12         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 4            | -3.98         |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 2            | -3.97         |
# +-----------------------------+--------------------------------+--------------+---------------+

```

Conclusion: Based on the previous results, Receptor Models 3 (CZP) and 4 (hCatL) refined receptors are selected for further docking studies


Compute the MMPBSA fingerprints for assays 2 and 4 
```python
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(assay_id=2,mmgbsa=1,prolif=0) # Input receptor field indexes when requested
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(assay_id=4,mmgbsa=1,prolif=0) # Input receptor field indexes when requested
```

Using SQL queries, inspect the results of the K777 MMGBSA total energy to CZP refined receptor
```bash
sqlite3 -header -table docking/docking_assays/assay_2/assay_2.db "SELECT sub_pose, delta_g_total FROM mmgbsa_fingerprints"

################
### Outputs ####
################
# +--------------------------------+---------------+
# |            sub_pose            | delta_g_total |
# +--------------------------------+---------------+
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | -50.2336      | # The reference pose is far more stable than the rest. Ref value: -50.2336
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | -38.5798      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | -42.6062      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | -39.7699      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | -34.289       |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | -39.1617      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | -34.7015      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | -31.1668      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | -32.5242      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | -35.3487      |
# +--------------------------------+---------------+
```

Using SQL queries, inspect the results of the K777 MMGBSA total energy to hCatL refined receptor
```bash
sqlite3 -header -table docking/docking_assays/assay_4/assay_4.db "SELECT sub_pose, delta_g_total FROM mmgbsa_fingerprints"

################
### Outputs ####
################
# +--------------------------------+---------------+
# |            sub_pose            | delta_g_total |
# +--------------------------------+---------------+
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | -43.6419      | # The reference pose is not more stable than the rest. Ref value: -43.6419
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | -36.1583      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | -40.368       |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | -36.2159      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | -31.9582      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | -34.0084      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | -43.9791      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | -30.7953      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | -31.2136      |
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | -22.9392      |
# +--------------------------------+---------------+
```

