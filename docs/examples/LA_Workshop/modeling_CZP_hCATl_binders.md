---
title: 3. Receptor modeling and validation
---

Before moving on to large-scale virtual screening, it is essential to validate the docking conditions on a trusted reference system. In this section, we will **prepare and import the receptor structures** used throughout the workshop, and perform a **validation docking experiment** using *K777* as reference.

---

## 1️⃣ Retrieve and process raw crystallographic structures

Experimental PDB files from the [Protein Data Bank (PDB)](https://www.rcsb.org/) are the starting point for most SBDD projects. However, these files usually cannot be used directly for molecular docking. Experimental PDB structures often contain multiple artifacts or missing data that can interfere with docking calculations, such as missing hydrogens, alternate conformations, crystallographic waters, co-crystallized ligands, and/or unresolved residues.  

Therefore, the first step is to **process and clean** the raw crystal structure to generate a ready-to-dock receptor model in the correct format.

Let's process the CZP crystal structure as a first example

We will use the CZP crystal structure with PDB code **2OZ2**, which is bound to **K777**. To retrieve the structure, execute the following commands in a bash terminal:

```bash
# Retrieve a structure corresponding to CZP
$ mkdir -p /PATH/TO/PROJECT/docking/raw_data/2OZ2
$ wget -O /PATH/TO/PROJECT/docking/raw_data/2OZ2/2OZ2.pdb https://files.rcsb.org/download/2OZ2.pdb
```

Alternatively you can download 2OZ2.pdb from [here](/downloads/2OZ2.pdb) (*right-click and download*).

The `process_raw_pdb()` function from **TidyScreen** automates all of the steps of receptor preparation for AutoDock in one integrated pipeline.   

#### 💬 Interactive prompts

By executing it, you will also be prompted to decide:

- Short description of the receptor model;  
`Please provide a brief description to store with the receptor model:`

- Which **chain(s)** to retain (e.g., `'A'` for this workshop, `'all'` to keep all chains);  
`Chains found in the pdb file: {'C', 'A'}`  
`Chains detected: {'C', 'A'}`  
`Provide the chain identifier to maintain ('all' to keep all chains):`

- If you wish to **preserve non-standard residues** - for example, co-crystallized ligands like *K777* (`D1R` residue code), and if you want to extract such ligands as **reference files** for later analysis or redocking validation  
`Non-standard residues detected:`
`The following non-standard residues were found in the pdb file: ['D1R', 'SO4'].`    
`Do you want to mantain ONE of them as a REFERENCE pdb file? (y/n):`
    - Type 'y' to keep K777.  

    `Type the 3-letter code of the residue you want to save as REFERENCE FILE:`  

    - Type 'D1R' to retain K777 as reference ligand.    

    - Afterwards, type 'n' for others like SO4, water, or cosolvent molecules (as many times as non-standard residues exists).

In particular, Molecular docking with **AutoDock4** requires receptor and ligand files in the `.pdbqt` format. This format extends the standard `.pdb` by including **partial atomic charges** and **atom types** defined by the AutoDock force field.  

Moreover, **AutoDock4** explores possible ligand poses inside a predefined **grid box**. This box must be centered in the **region of pharmacological interest** (in our case, the *catalytic site* of both enzymes).  

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/docking_box.png" alt="Description of image" width="400"/>
  <figcaption> Visualization of the AutoDock4 grid box covering the active site of CZP, defining the 3D region explored during ligand docking simulations. </figcaption>
  </p>
</figure>
---

To define these grids, we need to specify the dimensions of the box:

- The **central coordinates** (`x`, `y`, `z`)  
- The **number of grid points** in each direction (`x_points`, `y_points`, `z_points`)  

These values determine the volume within which the ligand can move and bind.

```python title="workshop.py"
from tidyscreen import tidyscreen as ts
from tidyscreen.moldock import moldock as md

# Activate the project
la_workshop = ts.ActivateProject("la_workshop_2025")

# Instantiate the MolDock class to manage receptor preparation
la_workshop_moldock = md.MolDock(la_workshop)

# Process the CZP crystal structure
structure_file = "/PATH/TO/PROJECT/docking/raw_data/2OZ2/2OZ2.pdb"
la_workshop_moldock.process_raw_pdb(
    structure_file, 
    x_coord=6.7, y_coord=-3.7, z_coord=-22.0,  # Center of the catalytic pocket
    x_points=40, y_points=40, z_points=40      # Grid box dimensions
)
```

Upon successful processing, the following files will be automatically created:

`/PATH/TO/PROJECT/docking/raw_data/2OZ2/D1R.pdb` ###Fredy, esto es así?  
`/PATH/TO/PROJECT/docking/raw_data/2OZ2/receptor.mol2`  
`/PATH/TO/PROJECT/docking/raw_data/2OZ2/receptor.pdbqt`  

In addition to the crystallographic binding pose of K777 (`D1R.pdb`), you get the cleaned, ready-to-dock receptor models compatible with AutoDock4. The `.mol2` file retains atom typing and bonding information, while the `.pdbqt` file will be used in all docking simulations.

:::tip[💡 Tip]
Always inspect the generated receptor visually (for instance using **VMD**, **PyMOL** or **UCSF Chimera**) to confirm that the structure is well-prepared.
:::


## 2️⃣ Compute AutoDock4 grid files

Once the receptor has been properly cleaned and converted into **AutoDock4-compatible formats**,  
the next step is to generate the **grid maps** with `autogrid4`. These define how the docking algorithm evaluates ligand-receptor interactions within the specified grid box.

AutoDock4 does not calculate binding interactions on the fly; instead, it precomputes energy maps for each atom type within the defined 3D grid. This drastically speeds up the docking process, allowing the program to focus on exploring ligand conformations and orientations. In addition to atom-specific energy maps, electrostatic and desolvation grids are also generated to account for polar interactions and solvation effects.

To compute the corresponding grid maps, navigate to your receptor directory and execute:

```bash
$ cd /PATH/TO/PROJECT/docking/raw_data/2OZ2/
$ autogrid4 -p receptor.gpf -l receptor.glg
```

After execution, several `.map` files will be generated, one per atom type (e.g., C, N, O, S),
along with the log file `receptor.glg`, all stored in the same receptor folder.

:::tip[💡 Tip]
If you modify the grid box size or center, remember to regenerate these maps before running new docking assays. Otherwise, the docking search may be performed in an incorrect spatial region.
:::

<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/docking_grids.png" alt="AutoDock4 grid maps and docking box on Cruzipain catalytic site" width="400"/>
  <figcaption>
  Visualization of the AutoDock4 grid maps. Colored mesh surfaces represent specific energy maps (OA=white, HD=iceblue, electrostatics=orange), precomputed within the 3D grid box used during docking simulations.
  </figcaption>
  </p>
</figure>


You should now have a fully prepared receptor folder similar to the one provided [here](/downloads/2OZ2.zip) (Right-click and save if you wish to explore the precomputed files). ###Fredy, se me descargan archivos vacíos.


## 3️⃣ Import precomputed receptor models into TidyScreen

In some cases, receptor preparation might have been performed outside of TidyScreen to generate all required files (`.pdbqt`, `.gpf`, `.map`, `.glg`, etc.). Instead of repeating those steps, you can **directly import the complete receptor folder** into your active TidyScreen project. This integrates the receptor into the project’s internal database and makes it immediately available for docking assays.

Use the following command to register the precomputed receptor in your active project:

```python title="workshop.py"
la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/2OZ2/")
```
During execution, you’ll be prompted to provide a brief description for the receptor model; for example: *“CZP raw crystal structure with precomputed grid maps”*.

*Outputs*: `Provide a brief description of the receptor model: # Input the required info`  
`Succesfully stored receptor model`  
`Successfully stored the receptor model in the database.`  

To accelerate the workflow, we will directly import a set of three receptor models that have been preprocessed and validated for docking. This will allow us to perform **comparative virtual screening** between the *therapeutic target* (CZP) and the *off-target homolog* (human Cathepsin L, hCatL).

These models include:
* **hCatL X-ray structure**, crystal structure co-crystallized with K777 (PDB ID: 8HFV) ([download zipped folder](/downloads/8hfv.zip)).
* **Refined CZP structure**, obtained after molecular dynamics (MM-MD) symulations ([download zipped folder](/downloads/CZP_refined.zip)).
* **Refined hCatL structure**, optimized through MM-MD simulation to improve binding site definition ([download zipped folder](/downloads/hCatL_refined.zip)).

After downloading and extracting these folders inside `/PATH/TO/PROJECT/docking/raw_data/`, import each receptor into TidyScreen:

```python title="workshop.py"
la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/8hfv")
la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/CZP_refined")
la_workshop_moldock.input_receptor("/PATH/TO/PROJECT/docking/raw_data/hCatL_refined")
```
###Fredy, los tres archivos vacíos.

Each receptor will be automatically indexed and stored in the project database, ready to be selected during docking setup.

## 5️⃣ Ligand preparation: K777  

Before running large-scale docking campaigns, it is essential to **validate the protocol** using reference ligands whose binding mode are already known. In this case, we will use **K777**, the well-characterized *vinyl sulfone* inhibitor co-crystallized with both CZP and hCatL. Reproducing its binding pose will help ensure that the docking parameters and receptor grids are correctly configured.  

To do that, in the `$PATH/chemspace/raw_data/` create a file called `K777.csv` containing the SMILES of the ligand.

```bash
cd $PATH/chemspace/raw_data/
echo "O=S(C1=CC=CC=C1)(/C=C/[C@H](CCC2=CC=CC=C2)NC([C@H](CC3=CC=CC=C3)NC(N4CCN(C)CC4)=O)=O)=O,K777,2" > K777.csv
```

This file includes the SMILES, ligand name, and an arbitrary flag column (here set to 2, the IC50 against CZP).

Then, activate the chemspace module within `workshop.py` and import the ligand into the `chemspace.db`.

```python title="workshop.py"
from tidyscreen.chemspace import chemspace as cs
la_workshop_cs = cs.ChemSpace(la_workshop)

csv_file = "/PATH/TO/PROJECT/chemspace/raw_data/K777.csv"
la_workshop_cs.input_csv(csv_file)
```

After preparing the receptor models, the next step is to generate all necessary *ligand files*. AutoDock4 requires ligands in `.pdbqt` format (which includes partial charges and atom types). Charges are assigned using the `prepare_ligand4.py` utility (default: Gasteiger), but this can be customized through the `charge_method` parameter. In this workflow, we rely on [`Meeko`](https://meeko.readthedocs.io/en/release-doc/), an RDKit-based toolkit that modernizes ligand preparation for AutoDock.

```python title="workshop.py"
la_workshop_cs.generate_mols_in_table(
    "K777", 
    conf_rank=25, 
    timeout=300, 
    pdbqt_method="meeko", 
    charge_method="gas"
)
```

This step generates multiple conformers (up to 25), optimizes their geometry, assigns charges, and prepares ready-to-dock input files within the `chemspace.db`.

At this stage, all the necessary input files have been prepared, and we are ready to proceed with the docking study.

## 6️⃣ Setting and launching molecular docking studies

Before running the docking experiment, it is necessary to define the **docking parameters**.  
The function `create_docking_params_set()` generates a database (`docking/params/docking_params.db`) containing the **default AutoDock conditions**.  

```python title="workshop.py"
la_workshop_moldock.create_docking_params_set()
```

Once the default set has been created, you may duplicate and customize the conditions to better suit your study. 
For example, you can increase the number of docking runs from 20 to 100 (*--nrun=100*) or enable AutoDock to identify and report ligand–receptor interactions for each pose (*--contact_analysis=1*).

If you prefer a graphical interface for editing these settings, you can open the parameters database using `sqlite_web`; a lightweight, browser-based SQLite viewer. This allows you to inspect, duplicate, and modify entries interactively, making it easier to fine-tune docking conditions without using SQL commands.

``` bash
sqlite_web /PATH/TO/PROJECT/docking/params/docking_params.db
```

Now we can launch the docking experiment with the `dock_table()` function. 

You must define 3 variables:
- *table_name*: name of the table containing the ligand(s) you want to use for the docking study, as stated in the `chemspace/processed_data/chemspace.db` database.
- *id_receptor_model*: ID of the target model you want to use for the docking study, as stated in the `docking/receptors/receptor.db` database.
- *id_docking_params*: ID of the docking conditions you want to apply for the study, as stated in the `docking/params/docking_params.db` database.

Once again you will be asked to provide a brief **description** of the docking assay. Try to include as many relevant details as possible.

:::note[📝 IDs]
Make sure that the IDs (*id_receptor_model* and *id_docking_params*) correspond to the entries you created earlier.
If in doubt, you can query the databases to list available IDs before running the docking.
:::


```python title="workshop.py"
la_workshop_moldock.dock_table(table_name="K777", id_receptor_model=1, id_docking_params=2)
```

After running `dock_table()`, *TidyScreen* will:

- Create a summary entry of the assay (variables + description) in `docking/docking_registers/docking_registers.db`.
- Create a new assay folder under `docking/docking_assays/assay_<ASSAY_ID>` containing:
  - *receptor/* and *ligands/* subfolders with all input files,
  - a ready-to-run execution script: `docking_execution.sh`.


You can launch the assay in your terminal with:

```jsx
./docking_execution.sh
```
:::note[📝 AutoDock]
TidyScreen prepares all required input files and executes the docking through AutoDock-GPU, the GPU-accelerated implementation of AutoDock4. This allows a substantial speed-up compared to the CPU version, while maintaining compatibility with the same scoring function and parameter set.
:::

## 7️⃣ Docking results analysis


#### Autodock output

In this first example, analysing the `.dlg` file and the docked poses is relatively straightforward, since the study involved only a single ligand. However, when working with larger ligand sets, manual inspection quickly becomes impractical.  For this reason, *TidyScreen* provides a `DockingAnalysis` module and a `process_docking_assay()` function to automate the post-processing of docking results.


```python title="workshop.py"
#Activate the module
from tidyscreen.docking_analysis import docking_analysis as docking_analysis
la_workshop_docking_analysis = docking_analysis.DockingAnalysis(la_workshop)
```

The `process_docking_assay()` function parses the AutoDock outputs of a given assay, leveraging the [Ringtail](https://ringtail.readthedocs.io/en/latest/) library to efficiently handle the results. It ranks the poses by docking score (most negative = best) and stores the selected results for downstream analysis.

You should define the following variables:

- *assay_id* (int, required): ID of the docking assay to analyze (as registered by dock_table()).
- *max_poses* (int, default = 10): Maximum number of poses per ligand to store in the results database, ordered from best (lowest binding energy) to worse.
- *extract_poses* (int [0,1], default = 0): If 1, writes PDB files for each selected pose. If 0, keeps only the database records and existing PDBQT/DLG outputs.


```python title="workshop.py"
#Set and run the docking analysis
la_workshop_docking_analysis.process_docking_assay(assay_id=1, max_poses=10, extract_poses=1)
```

As a result, a new database `assay_<ID>.db`  will be generated inside the corresponding assay docking folder.

This file contains several tables that organize the results of the analysis:

* *interaction_indices*: dictionary of interaction types (e.g., hydrogen bonds H and van der Waals V) that can be identified between ligands and receptor.

* *interactions*: records which interactions were detected for each docked pose of every ligand.

* *ligands*: stores ligand-specific data and a summary of their detected interactions.

* *receptor*: summary information of the target used in the assay.

* *Results*: pose-level data extracted from the .dlg file, including binding energies, cluster ranks, and other docking details.

If *extract_poses=1*, a new subfolder *docked_1_per_cluster/* will be created inside the assay folder, containing the PDB files of the top docked poses (one representative per cluster) for each ligand, making it easier to visualize and compare binding modes.

#### Interpreting docking outcomes

Let us start by analyzing the results obtained. The following query lists all poses ranked by their docking score:

```bash
sqlite3 -header -table docking/docking_assays/assay_1/assay_1.db \
"SELECT LigName, sub_pose, cluster_size, docking_score FROM Results;"
```

*Output*

| **LigName**                      | **sub_pose**                        | **cluster_size** | **docking_score** |
|----------------------------------|-------------------------------------|------------------|-------------------|
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 10 | -4.94 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 5  | -2.24 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 3  | -1.64 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 1  | -1.54 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 3  | -1.38 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 2  | -1.33 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 7  | -1.06 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 1  | -0.86 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 6  | -0.80 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 3  | -0.76 |


Despite producing a reasonable number of clusters, the top-ranked pose (-4.94 kcal.mol<sup>-1</sup>) does not reproduce the experimental binding mode observed in the crystal. When visualized, the ligand often appears misoriented with respect to the crystallographic pose, with its warhead shifted away from the catalytic cysteine.

In this case, the receptor derived directly from the crystal was used “as is,” without prior optimization, resulting in a non-productive binding mode for K777. 

So that, next we will repeat this same redocking using the refined CZP receptor, which was previously equilibrated via molecular dynamics, and confirm that the predicted pose now matches the crystal orientation of K777.


```python title="workshop.py"
la_workshop_moldock.dock_table(table_name="K777", id_receptor_model=1, id_docking_params=2)

la_workshop_docking_analysis.process_docking_assay(assay_id=2, max_poses=10, extract_poses=1)
```

Query the database to inspect the new docking results:

```bash
sqlite3 -header -table docking/docking_assays/assay_2/assay_2.db \
"SELECT LigName, sub_pose, cluster_size, docking_score FROM Results;"
```

*Output*

| **LigName**                      | **sub_pose**                        | **cluster_size** | **docking_score** |
|----------------------------------|-------------------------------------|------------------|-------------------|
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 29 | -5.42 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 7  | -5.34 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 2  | -4.64 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 1  | -4.61 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 5  | -4.61 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 2  | -4.55 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 7  | -4.34 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 1  | -4.31 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 4  | -4.30 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 2  | -4.28 |


Unlike the previous case, the top-ranked pose (-5.42 kcal.mol<sup>-1</sup>) now reproduces the experimental binding mode of K777 with excellent geometric accuracy. This improved fit reflects the enhanced flexibility and energetic realism of the refined receptor, which allows the pocket to adapt to the inhibitor and minimizes steric strain. This result validates both the receptor preparation pipeline and the docking parameters, confirming that the workflow can be safely applied to screen new ligands.

#### Cross-docking validation on human Cathepsin L (hCatL)

To further validate the docking workflow and evaluate **selectivity**, we repeated the same procedure against hCatL. 

* Docking using the raw X-ray hCatL structure

```python title="workshop.py"
la_workshop_moldock.dock_table("K777", id_receptor_model=2, id_docking_params=2)
la_workshop_docking_analysis.process_docking_assay(assay_id=3, max_poses=10, extract_poses=1)
```
Inspect the results:

```bash
sqlite3 -header -table docking/docking_assays/assay_3/assay_3.db \
"SELECT LigName, sub_pose, cluster_size, docking_score FROM Results;"
```

*Output*

| **LigName**                      | **sub_pose**                        | **cluster_size** | **docking_score** |
|----------------------------------|-------------------------------------|------------------|-------------------|
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 6  | -4.08 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 7  | -3.93 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 3  | -2.97 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 4  | -2.62 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 5  | -2.00 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 8  | -1.98 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 4  | -1.52 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 1  | -1.33 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 4  | -1.25 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 3  | -0.75 |


In this case, the top-scoring poses (≈ -4.0 kcal.mol<sup>-1</sup>) failed to reproduce the experimental covalent binding orientation of K777. This mismatch likely may reflect the rigid conformation of the crystal structure, which was obtained with a covalently bound inhibitor, thus restricting access to the catalytic cleft for non-covalent docking.

* Docking using the refined hCatL structure 

We then repeated the docking using the refined hCatL receptor (ID 4), previously equilibrated by molecular dynamics.

```python title="workshop.py"
la_workshop_moldock.dock_table("K777", id_receptor_model=4, id_docking_params=2)
la_workshop_docking_analysis.process_docking_assay(assay_id=4, max_poses=10, extract_poses=1)
```
Inspect the results:

```bash
sqlite3 -header -table docking/docking_assays/assay_4/assay_4.db \
"SELECT LigName, sub_pose, cluster_size, docking_score FROM Results;"
```
*Output*

| **LigName**                      | **sub_pose**                        | **cluster_size** | **docking_score** |
|----------------------------------|-------------------------------------|------------------|-------------------|
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1  | 22 | -5.17 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_2  | 7  | -5.01 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_3  | 1  | -4.92 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_4  | 4  | -4.67 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_5  | 7  | -4.58 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_6  | 8  | -4.32 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_7  | 1  | -4.25 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_8  | 2  | -4.12 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_9  | 4  | -3.98 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_10 | 2  | -3.97 |

The best pose (-5.17 kcal.mol<sup>-1</sup>) successfully reproduced the experimental covalent orientation. The improvement observed after receptor refinement highlights how structural relaxation enhances the accuracy of docking simulations. The refined hCatL pocket allowed better accommodation of K777, maintaining the proper orientation of the vinyl sulfone warhead and hydrogen-bonding pattern characteristic of the experimental complex.

Together with the results obtained for CZP, these findings confirm that receptor models 3 (CZP) and 4 (hCatL) offer the best compromise between structural realism and computational tractability, and will be used in the subsequent large-scale virtual screening.

## 8️⃣ Post-docking refinement with MMGBSA calculations 

While docking scores provide a quick estimate of binding affinity, they rely on simplified scoring functions that may overlook dynamic-related effects.  
To obtain more physically grounded interaction energies, *TidyScreen* integrates an [**MMGBSA module**](https://pubmed.ncbi.nlm.nih.gov/25835573/) (*Molecular Mechanics Generalized Born Surface Area*) that re-evaluates the best poses considering implicit solvent effects and molecular mechanics terms.

We can compute MMGBSA energies for any docking assay using the function `compute_fingerprints_for_whole_assay()`:

```python title="workshop.py"
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(assay_id=2, mmgbsa=1, prolif=0)
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(assay_id=4, mmgbsa=1, prolif=0)
```

During execution, you’ll be prompted to specify which receptor fields to include in the energy calculation. Once completed, TidyScreen automatically stores the per-pose MMGBSA energy terms in the corresponding assay databases.

Let's inspect the MMGBSA total binding free energies for assay 2.

```bash
sqlite3 -header -table docking/docking_assays/assay_2/assay_2.db \
"SELECT sub_pose, delta_g_total FROM mmgbsa_fingerprints;"
```
*Output*

| **sub_pose**                     | **delta_g_total** |
|----------------------------------|-----------------------------|
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_1   | -50.2336 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_2   | -38.5798 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_3   | -42.6062 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_4   | -39.7699 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_5   | -34.2890 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_6   | -39.1617 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_7   | -34.7015 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_8   | -31.1668 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_9   | -32.5242 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_10  | -35.3487 |


The reference pose (-50.23 kcal.mol<sup>-1</sup>) is markedly more stable than the others, confirming a correct reproduction of the experimental conformation.

Now, for assay 4 (refined hCatL)...

```bash
sqlite3 -header -table docking/docking_assays/assay_4/assay_4.db \
"SELECT sub_pose, delta_g_total FROM mmgbsa_fingerprints;"
```
*Output*

| **sub_pose**                     | **delta_g_total** |
|----------------------------------|-----------------------------|
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_1   | -43.6419 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_2   | -36.1583 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_3   | -40.3680 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_4   | -36.2159 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_5   | -31.9582 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_6   | -34.0084 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_7   | -43.9791 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_8   | -30.7953 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_9   | -31.2136 |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N_10  | -22.9392 |

Unlike for CZP, the experimental reference pose (-43.64 kcal.mol<sup>-1</sup>) is not the most stable one according to MMGBSA.

These results illustrate a critical aspect for our project:
* For Cruzipain (CZP), the most negative docking score and the most stable MMGBSA energy coincide; both correspond to the bioactive conformation observed experimentally.
* For hCatL, however, MMGBSA ranking suggests slightly more favorable alternative poses that do not match the crystallographic orientation.

This discrepancy reflects that while MMGBSA provides a more detailed energetic profile, it is still dependent on the initial docking geometry and receptor conformation. In this context, prioritizing docking scores may be a more reliable first criterion for hit identification, at least within this small validation set. However, the relative reliability of docking versus MMGBSA can vary depending on the target flexibility, ligand series, and parameterization, etc.

Importantly, the MMGBSA values obtained here can serve as useful reference benchmarks for future comparisons, particularly when analyzing larger or chemically diverse ligand libraries.

A more statistically robust conclusion requires performing this workflow on a broader, chemically diverse dataset, which will be explored in the **next tutorial** on expanded training/test sets and pose validation pipelines.


:::note[Key takeaways]

At the end of this tutorial, you should be able to:

✅ **Prepare and clean receptor structures**, from raw crystallographic data to AutoDock-ready models.  
✅ **Define and visualize docking grids**, focusing the search on the catalytic site or region of pharmacological interest.  
✅ **Import receptor folders** directly into TidyScreen and manage multiple models (raw vs. refined).  
✅ **Generate ligand input files**.  
✅ **Set and customize docking parameters** using the integrated SQLite databases.  
✅ **Launch docking runs**.  
✅ **Analyze docking results** efficiently using the DockingAnalysis module and SQL queries.  
✅ **Compute MMGBSA energies** to refine binding affinity predictions and compare pose stabilities.  
✅ **Interpret results**, recognizing when each metric is most reliable.  
:::
