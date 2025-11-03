---
title: Chemical space design and vHTS workflow
---

# Chemical Space Design and vHTS Workflow

In this final section, we will extend the previously validated docking methodology to generate and evaluate **new triazole-based Targeted Covalent Inhibitors (TCIs)** of Cruzipain (CZP).

Following the successful validation of the docking workflow using a training set of triazole-containing derivatives featuring a methyleneamine linker as bioisosteric replacements of peptide bonds and the well-known efficacy of the vinylsulfone-containing inhibitor **K777**, our goal is to design a new series of molecules that **combine these three structural elements**.

Given the absence of reported analogues combining these motifs, we will perform a massive virtual library generation and screening campaign aimed at identifying candidates with **high affinity for CZP** and **minimal affinity for hCatL**, while ensuring **synthetic feasibility** for subsequent *wet-lab* validation.

---

## Overview of the Protocol

The overall workflow can be divided into five major steps:

1. **Extraction of building blocks**  
   Retrieve potential reactants from the *eMolecules* chemical database.

2. **Chemical space reduction**  
   Filter the combinatorial chemical space based on:
   - *Drug-likeness* criteria 
   - *Synthetic feasibility*
   - *Affordability*, using models provided through **Ersilia Hub**.

3. **In silico synthesis of derivatives**  
   Combine the selected reactants through feasible *in silico* synthetic routes to generate the final compounds, replicating realistic reactions that can be performed in the wet lab.

4. **Toxicity filtering**  
   Further reduce the library by eliminating potentially cytotoxic molecules, applying an **Ersilia Hub cytotoxicity prediction model**.

5. **Docking and analysis**  
   Perform docking-based *vHTS* on the safest and most promising analogues (approximately **5000 compounds**) and analyze their predicted affinity and selectivity profiles toward CZP and hCatL.

---

🎯 This will allow prioritization of candidates for *in vitro* and *in vivo* validation, completing the *dry*-to-*wet lab* pipeline proposed in this tutorial.

---

## 1. **Extraction of building blocks**  


#### 🧭 Synthetic plan

We will combine three independently tractable substituents (**R<sup>1</sup>**, **R<sup>2</sup>**, **R<sup>3</sup>**) on a common scaffold through a concise, wet-lab-ready sequence:

*i.)* **R<sup>1</sup>** will be introduced from commercial **amino acids**, which will be derivatized to the corresponding **azides** by means of *diazotransfer* reactions.

*ii.)* **R<sup>2</sup>** and **R<sup>3</sup>** will be derived from **aldehyde** and **methyleneamine** derivatives, respectively, which will be combined via *A3 coupling* to obtain **propargylamine** intermediates.

*iii.)* The azides and alkynes obtained from (i) and (ii) will be joined through *CuAAC* to yield **1,4-disubstituted 1,2,3-triazole** derivatives.

*iv.)* Finally, to introduce the **phenyl-vinylsulfone warhead (WH)**, a *Horner-Wadsworth-Emmons (HWE) olefination* will be replicated *in silico*.

---

Thus, we will generate lists of potential **reactants** - namely **amino acids**, **aldehydes**, and **methyleneamines** - that can be purchased and used as starting points for the *in silico* synthesis.  
To this purpose, we will filter compounds from the **eMolecules** database. This will ensure that all designed structures are **synthetically feasible and accessible** for *wet-lab* validation.


:::note[📝 About eMolecules]

**[eMolecules](https://www.emolecules.com/)** is a large supplier-aggregated catalog of **purchasable reagents** and **make-on-demand** building blocks. It was founded in 2005 with a vision to reduce drug discovery timelines through improved efficiencies in the compound search and acquisition process.
:::

:::warning[Important]

Please check if the *chemspace module* is active within `workshop.py`.

```python title="workshop.py"
la_workshop_cs = chemspace.ChemSpace(la_workshop)
```
:::

### 1.1. 🗂️ Import eMolecules database

Import the chemical database file from **eMolecules** (usually stored in `chemspace/raw_data/`) into the local **TidyScreen** project database (`chemspace/processed_data/chemspace.db`).

This step converts the original `.csv` file into a structured SQLite database that can be efficiently queried by TidyScreen for compound selection and filtering.

```python title="workshop.py"
emolecules_database_file = "/PATH/TO/FILE/emolecules.csv" #Set your own path!
la_workshop_cs.input_csv(emolecules_database_file)
```

### 1.2. 🧩 Chemical scaffolds filtering

To subset the types of compounds we need, **TidyScreen** uses chemical filters defined by **SMARTS patterns**.  
By default, a list of these filters is already included and can be inspected using the following function:


```python title="workshop.py"
# List available SMARTS filters to obtain reactants
la_workshop_cs.list_available_smarts_filters()
```
You should get the following output: 


```Available SMARTS filters:```  
```Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]```  
```...```  
```Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]```  
```Filter_id: 11, Filter_Name: PrimAmines, SMARTS: [NX3;H2;!$(NC=O)]```  
```Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]```  
```...```  
```Filter_id: 53, Filter_Name: sulfonamide, SMARTS: [SX4](=O)(=O)N```  

One of the classes of reactants we need corresponds to α-amino acids (`Filter_id = 1`).
While it is possible to directly select all eMolecules entries containing at least one such group, this approach would be inefficient and chemically irrational. The output would likely include compounds with
multiple amino acid groups (unsuitable for diazotransfer), additional reactive functional groups interfering with subsequent steps, rare isotopes or heteroatoms uncommon in drug-like molecules, chemically unstable or undesirable moieties, etc.

To refine the chemical space, TidyScreen allows creating sequential filtering workflows combining multiple SMARTS rules.

The function `create_smarts_filters_workflow()` accepts a dictionary in the form:

`{ filter_id : maximum_allowed_occurrences }`

Each key-value pair defines how many times a particular substructure is permitted in a molecule (e.g., 1:1 means exactly one α-amino acid group; 2:0 means exclude all boron-containing compounds).

In this example, we will filter molecules that:

- Contain exactly one α-amino acid group (1:1)
- Contain exactly one primary amine (11:1)
- Contain exactly one carboxylic acid group (26:1)
- Exclude molecules containing any of the following: Boron (2:0), Silicon (3:0), Azides (4:0), Terminal Alkynes (5:0), Deuterium (6:0), Tritium (7:0), <sup>13</sup>C (8:0), <sup>15</sup>N (9:0), Selenium (17:0), Amides (21:0), Thiols (25:0), Esters (35:0), Sulfonamides (53:0).


```python title="workshop.py"
# Create a custom filtering workflow for α-amino acids
la_workshop_cs.create_smarts_filters_workflow({1:1,11:1,26:1,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,17:0,21:0,25:0,35:0,53:0})
```

When prompted, you will need to provide a short **description of the filter** (i.e.: *"Filtering of amino acids compatible with CuAAC-centered reaction scheme"*).

`Provide a description for the SMARTS filters workflow: `

Once confirmed, the workflow is stored in the internal database.  
You can verify that it was correctly created using the `list_available_smarts_filters_workflows()` function:

```python title="workshop.py"
# List all available SMARTS filtering workflows
la_workshop_cs.list_available_smarts_filters_workflows()
```

*Output*: `Workflow_id: 1, Filter_Specs: {"1": 1, "11": 1, "26": 1, "2": 0, "3": 0, "4": 0, "5": 0, "6": 0, "7": 0, "8": 0, "9": 0, "17": 0, "21": 0, "25": 0, "35": 0, "53": 0}, Description: Filtering of amino acids compatible with CuAAC-centered reaction scheme`


Now we can apply the corresponding workflow to filter the desired building blocks using the function `subset_table_by_smarts_workflow()`. 
You need to specify two parameters: *The table name* (as stored in the chemspace.db) and the *Workflow ID* to be applied.

```python title="workshop.py"
# Filter α-amino acids
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",1)
```

When prompted, provide a short **description of the subset** (i.e.: *"Amino acids for diazotransfer reaction"*).

`Enter a description for the subset:`

A new table containing only the filtered compounds will be created in the `chemspace.db`. Each table follows the naming convention `emolecules_subset_X`, where "X" corresponds to the subset number (incremental per run), not to the Workflow ID. This ensures that your filtered set is stored as an independent subset, ready to be used for subsequent in silico synthesis steps.

`Table 'emolecules_subset_1' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'`  
`Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '1'`


Next, we will subset the two remaining classes of starting materials - **aldehydes** and **methyleneamines** - following the same procedure used for α-amino acids.

The following workflow will retain only molecules containing **one aldehyde group** (`Filter_id = 10`) while excluding other reactive or undesired functional groups.

```python title="workshop.py"
# Create a custom filtering workflow for aldehydes
la_workshop_cs.create_smarts_filters_workflow({10:1,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,11:0,12:0,17:0,21:0,25:0,26:0,35:0,53:0})

# Filter aldehydes
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",2)
```

`Table 'emolecules_subset_2' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'`  
`Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '2'`

The `emolecules_subset_2` table will now contain only aldehyde derivatives suitable for A3 coupling with methyleneamines.


For **methyleneamines**, we do not have a pre-defined SMARTS filter.  
However, TidyScreen provides a specific function for adding custom SMARTS-based filters: `add_smarts_filter()`.  
You simply need to specify the SMARTS pattern and a descriptive filter name.


```python title="workshop.py"
# Add a custom SMARTS filter for primary methylene amines
la_workshop_cs.add_smarts_filter("[NX3;H2][CX4;H2]","Primary_Amines_custom")
```

*Output*: `Successfully added SMARTS filter: '[NX3;H2][CX4;H2]'`

To confirm that the new SMARTS filter was added successfully, you can list all available filters again:

```python title="workshop.py"
# List available SMARTS filters to obtain reactants
la_workshop_cs.list_available_smarts_filters()
```

```Available SMARTS filters:```  
```Filter_id: 1, Filter_Name: Aminoacids, SMARTS: [NX3H2,NX4+H3][CX4H]([*])[CX3H0](=[OX1])[OX2H,OX1-]```  
```...```  
```Filter_id: 10, Filter_Name: Aldehydes, SMARTS: [CX3H1](=O)[#6]```  
```Filter_id: 11, Filter_Name: PrimAmines, SMARTS: [NX3;H2;!$(NC=O)]```  
```Filter_id: 12, Filter_Name: SecAmines, SMARTS: [NX3;H1;!$(NC=O)]```  
```...```  
```Filter_id: 53, Filter_Name: sulfonamide, SMARTS: [SX4](=O)(=O)N```  
```Filter_id: 54, Filter_Name: Primary_Amines_custom, SMARTS: [NX3;H2][CX4;H2]```  

Now that the new filter is available, we can create and apply a corresponding workflow to subset the desired primary methyleneamine derivatives:

```python title="workshop.py"
# Create a custom filtering workflow for primary methylene amines
la_workshop_cs.create_smarts_filters_workflow({54:1,10:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,12:0,17:0,21:0,25:0,35:0,53:0})

# Filter primary methylene amines
la_workshop_cs.subset_table_by_smarts_workflow("emolecules",3)
```

`Table 'emolecules_subset_3' created in: '/PATH/TO/PROJECT/chemspace/processed_data/chemspace.db'`  
`Succesfully subseted table: 'emolecules' by SMARTS filters workflow with ID: '3'`

At this point, your chemspace database contains three curated subsets ready for in silico synthesis.

### 1.3. 📸 Inspection of a sample of filtered building blocks

Before moving forward, it is useful to visually inspect a random sample of compounds from each curated subset to ensure that the filtering process performed as expected.

TidyScreen provides the function `depict_ligand_table()`, which automatically generates molecular depictions for a given table in the **chemspace database** and saves them in the `processed_data/misc` directory.


```python title="workshop.py"
# Depict a random sample of amino acids
la_workshop_cs.depict_ligand_table("emolecules_subset_1", limit=25, random=True) # Outputs to '/PATH/TO/PROJECT/chemspace/processed_data/misc'
# Depict a random sample of aldehydes
la_workshop_cs.depict_ligand_table("emolecules_subset_2", limit=25, random=True) # Outputs to '/PATH/TO/PROJECT/chemspace/processed_data/misc'
# Depict a random sample of primary methyleneamines
la_workshop_cs.depict_ligand_table("emolecules_subset_3", limit=25, random=True) # Outputs to '/PATH/TO/PROJECT/chemspace/processed_data/misc'
```

## 2. Building blocks prioritization

In order to narrow down the huge number of potential reactants, we will compute a set of standard **drug-like properties** and apply **Ersilia models** for each curated subset and then **prioritize** entries by applying simple, reproducible property filters.


### 2.1. ⚗️ Subsetting by drug-like properties (RDKit)

TidyScreen relies on the open-source cheminformatics toolkit **[RDKit](https://www.rdkit.org/)** to calculate commonly used molecular descriptors directly from SMILES strings.  
The following properties are computed by default running the `compute_properties()` function:

- **MolWt** → Molecular weight  
- **MolLogP** → Octanol-water partition coefficient (logP)  
- **NumHDonors** → Number of hydrogen bond donors  
- **NumHAcceptors** → Number of hydrogen bond acceptors  
- **NumRotatableBonds** → Number of rotatable bonds  
- **TPSA** → Topological polar surface area

These descriptors are standard in medicinal chemistry and provide a simple first layer of prioritization to retain drug-like, synthetically tractable molecules.

```python title="workshop.py"
# α-Amino acids
la_workshop_cs.compute_properties("emolecules_subset_1") 
# Aldehydes
la_workshop_cs.compute_properties("emolecules_subset_2") 
# Methyleneamines
la_workshop_cs.compute_properties("emolecules_subset_3") 
```

:::warning[!!!]

Please note that when running the `compute_properties()` function you'll get the following warning: 

`The table named 'emolecules_subset_X' already exists in the database. I will replace it, are you ok with that? (y/n):`.

When prompted, type *`y`* to update the tables with the newly computed properties.
:::

💡 Now that all tables contain the computed molecular descriptors, you can explore their value distributions to get a sense of the chemical diversity within each subset.  
This helps identify whether your current chemical space fits the intended *drug-like* range before applying filters.

With all molecular descriptors available, we can constrain the chemical space to **compounds within defined physicochemical ranges**. This ensures that the selected building blocks remain within desirable regions of drug-like chemical space, minimizing the inclusion of overly flexible, lipophilic/hydrophilic, or bulky structures that may hinder synthesis or bioavailability.

To perform this filtering, TidyScreen provides the function `subset_table_by_properties()`, which allows you to select compounds that meet specific descriptor criteria.  

In this example, we keep entries with:

- 200 ≤ MolWt ≤ 500
- 1.5 ≤ MolLogP ≤ 3.0
- NumRotatableBonds ≤ 2

To retain compact, moderately lipophilic building blocks with low flexibility.


```python title="workshop.py"
drug_like_filters=["MolWt>=200","MolWt<=500","MolLogP>=1.5","MolLogP<=3","NumRotatableBonds<=2"]
# α-Amino acids
la_workshop_cs.subset_table_by_properties("emolecules_subset_1",drug_like_filters)
# Aldehydes
la_workshop_cs.subset_table_by_properties("emolecules_subset_2",drug_like_filters)
# Methyleneamines
la_workshop_cs.subset_table_by_properties("emolecules_subset_3",drug_like_filters)
```

:::note[Output]
When using the `subset_table_by_properties()` function, you will be prompted to enter a **short description** of the new subset. Being as detailed as possible will help you avoid misunderstandings and ensure clear tracking and reproducibility.
:::

After running the filtering step, TidyScreen automatically creates **new tables** inside your project database (`chemspace.db`). Each of them contains only the compounds that **match all defined property conditions**, effectively reducing the chemical space to its most promising and drug-like entries.  

You can identify these new tables by their names (e.g., `emolecules_subset_1_subset_4`, `emolecules_subset_2_subset_5`, etc.), each corresponding to the filtered version of the original subsets.

### 2.2. 💸 Price-based prioritization with Ersilia

To further improve **synthetic feasibility**, we will estimate reagent prices using an *Ersilia Open Source Initiative* model and then subset our building blocks by **price ranges**.  
The model `eos7a45` (CopriNet-style predictor) adds a new column to your table (i.e., `eos7a45_coprinet`) with an estimated price score.

```python title="workshop.py"
# α-Amino acids
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_1_subset_4", "eos7a45")

# Aldehydes
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_2_subset_5", "eos7a45")

# Methyleneamines
la_workshop_cs.apply_ersilia_model_on_table("emolecules_subset_3_subset_6", "eos7a45")
```

To inspect the range of predicted prices (min-max) per subset and set realistic thresholds, we can use  `sqlite3` in a *bash terminal*.

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db \
"SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_1_subset_4 WHERE eos7a45_coprinet IS NOT NULL;"

sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db \
"SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_2_subset_5 WHERE eos7a45_coprinet IS NOT NULL;"

sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db \
"SELECT MIN(eos7a45_coprinet), MAX(eos7a45_coprinet) FROM emolecules_subset_3_subset_6 WHERE eos7a45_coprinet IS NOT NULL;"
```

You should get the following output: 

```2.66|5.60 #Amino acids```  
```1.24|7.24 #Aldehydes```  
```2.19|7.09 #Methyleneamines```  

💬 **How to interpret these values**

The `eos7a45_coprinet` field represents a **predicted cost index ($/mmol)** rather than an absolute monetary price.
Lower scores indicate cheaper, easier-to-source reagents, while higher scores correspond to rarer or more complex compounds.

These values help define cutoff thresholds that balance cost and chemical diversity, guiding which compounds will be retained for downstream in silico synthesis.

🧾 **Create tables of “cheap, drug-like” building blocks**


After computing predicted prices with the Ersilia model, we can further refine our chemical space by retaining only those compounds that are both drug-like and synthetically affordable.  
This additional filter helps ensure that the selected reagents are not only *chemically sound* but also *economically realistic* for laboratory synthesis - a key factor when planning large-scale combinatorial campaigns.

In this step, we will use *bash commands* to interact directly with the project’s SQLite database and create new tables containing only the **cost-effective reagents**, using `sqlite3` to query the project database (`chemspace.db`).

Each line below:
1. Selects only molecules within the specified price range (`BETWEEN X AND Y`), and  
2. Stores the result as a new table labeled `*_druglike_cheap`, ready for the next stages of *in silico* library construction.

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aminoacids_druglike_cheap AS SELECT * FROM emolecules_subset_3_subset_6 WHERE eos7a45_coprinet BETWEEN 1 AND 3.5 AND eos7a45_coprinet IS NOT NULL;" # Returns 83 amino acid building blocks

sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE aldehydes_druglike_cheap AS SELECT * FROM emolecules_subset_1_subset_4 WHERE eos7a45_coprinet BETWEEN 1 AND 2.3 AND eos7a45_coprinet IS NOT NULL;" # Returns 135 aldehyde building blocks

sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "CREATE TABLE amines_druglike_cheap AS SELECT * FROM emolecules_subset_2_subset_5 WHERE eos7a45_coprinet BETWEEN 1 AND 2.6 AND eos7a45_coprinet IS NOT NULL;" # Returns 138 methyleneamine building blocks
```

You now have three price-filtered, drug-like reagent pools — `aldehydes_druglike_cheap`, `amines_druglike_cheap`, and `aminoacids_druglike_cheap`.  
These curated sets combine **structural suitability**, **favorable physicochemical properties**, and **economic accessibility**, ensuring a realistic starting point for *in silico* combinatorial synthesis in the next module.

## 3. Combinatorial virtual synthesis

In this section, we will virtually replicate the multistep synthetic route that combines the three types of prioritized building blocks - **amino acids, aldehydes, and methylene amines** - to generate the desired **triazole-based potential TCIs**.

#### Reaction scheme overview
---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/Synthetic_Scheme_Whorkshop.png" alt="Description of image" width="900"/>
  <figcaption>General synthetic workflow leading to triazole-based inhibitors. </figcaption>
  </p>
</figure>
---


The synthetic plan consists of six consecutive transformations:
1. **Acylation** to prepare an ester intermediate from α-amino acids.  
2. **Diazotransfer** to generate α-azido esters.  
3. **A3 coupling** to afford *propargylamines*.  
4. **CuAAC reaction** (*click chemistry*) to form 1,4-disubstituted 1,2,3-triazoles.  
5. **DIBAL reduction** to partially reduce carbonyl intermediates.  
6. **Horner–Wadsworth–Emmons olefination** to append the *phenyl-vinylsulfone warhead*.

#### SMARTS-based combinatorial synthesis

TidyScreen uses **SMARTS** to define chemical transformations in a machine-readable way.  Each reaction rule is a *pattern-based instruction* describing how specific functional groups in reactants (defined before `>>`) are converted into new substructures in products (defined after `>>`).

This approach enables **combinatorial virtual synthesis**, allowing you to simulate multi-step reaction sequences across large datasets of reagents to produce *synthetically plausible* compound libraries.

### 3.1. 🧠 Defining SMARTS reactions

We begin by checking whether any reactions are already defined. 
List available SMARTS reactions (empty at first):

```python title="workshop.py"
la_workshop_cs.list_available_smarts_reactions()
```

*Output*: `SMARTS reactions table does not exist yet. Add reactions to the database first`.

Let’s now populate the database with the six transformations corresponding to our synthetic route:

Add the SMARTS reaction for diazotransfer
```python title="workshop.py"
la_workshop_cs.add_smarts_reaction("[CX3:1](=[O:2])[OX2H,OX1-:3]>>[C:1](=[O:2])[O:3][C]", "Acylation")

la_workshop_cs.add_smarts_reaction("[NX3;H2:1][CX4:2][CX3,H0:3]>>[N-]=[N+]=[NX2;H0:1][CX4:2][CX3,H0:3]", "Diazotransfer")

la_workshop_cs.add_smarts_reaction("[N:1].[CX3H1:2](=[O:3])>>[NX4+:1][C@H:2][C:4]#[C:5]", "A3 coupling")

la_workshop_cs.add_smarts_reaction("[NX1-:1]=[NX2+:2]=[NX2:3].[CX2H1:4]#[CX2H0:5]>>[NX2+0:1]1=[NX2+0:2][N:3]-[C:4]=[C:5]1", "CuAAC")

la_workshop_cs.add_smarts_reaction("[CX3:1](=[O:2])[OX2,OX1-:3]>>[CX3H1:1](=[O:2])", "DIBAL reduction")

la_workshop_cs.add_smarts_reaction("[CX4:1][CX3H1:2](=[O:3])>>[CX4:1]\[CX3H1:2]=[CX3H1]\[S](=O)(=O)c1ccccc1", "Horner_Wadsworth_Emmons Olefination")
```

We can verify that all reactions were correctly registered:

```python title="workshop.py"
la_workshop_cs.list_available_smarts_reactions()
```

*Output*:  
`Available SMARTS reactions:`  
`Reaction_id: 1, Name: Acylation`  
`Reaction_id: 2, Name: Diazotransfer`  
`Reaction_id: 3, Name: A3 coupling`  
`Reaction_id: 4, Name: CuAAC`  
`Reaction_id: 5, Name: DIBAL reduction`  
`Reaction_id: 6, Name: Horner_Wadsworth_Emmons Olefination`  

### 3.2. 🧪 Creating and applying a reaction workflow

TidyScreen allows the definition of multi-step reaction workflows, enabling sequential application of transformations to progressively build a virtual library.

Here, we will chain the six reactions defined above into a single workflow that mimics our synthetic plan:

```python title="workshop.py"
la_workshop_cs.add_smarts_reaction_workflow([1,2,3,4,5,6])
```

List the available workflows to confirm creation:

```python title="workshop.py"
la_workshop_cs.list_available_reactions_workflows()
```

Now, execute the workflow, specifying which reagent sets participate in each stage. The order and arrows (`->:-1`, etc.) indicate how intermediate products feed into the next transformation.

```python title="workshop.py"
la_workshop_cs.apply_reaction_workflow(1,           # Workflow ID 
    [   ["aminoacids_druglike_cheap"],              # Step 1: Acylation
        ["->:-1"],                                  # Step 2: Diazotransfer
        ["amines_druglike_cheap", "aldehydes_druglike_cheap"],  # Step 3: A3 coupling
        ["->:-2", "->:-1"],                         # Step 4: CuAAC
        ["->:-1"],                                  # Step 5: DIBAL reduction
        ["->:-1"]                                   # Step 6: HWE olefination
    ])
```

This process may take several minutes depending on the number of reactants and system resources.
All generated intermediates and final products are automatically stored in the project’s chemspace database under the name `reaction_set_1`.

To get a sense of the diversity and plausibility of your newly generated compounds, you can visualize a random subset of the virtual triazole derivatives:

```python title="workshop.py"
la_workshop_cs.depict_ligand_table("reaction_set_1", limit=25, random=True)
```

The output will generate a grid of 2D structures and save it in `/PATH/TO/PROJECT/chemspace/processed_data/misc`

This figure provides a quick visual overview of the chemical diversity, substitution patterns, and structural coherence of the in silico library; a valuable sanity check before proceeding to molecular docking.

✅ You have now constructed a multi-step combinatorial virtual synthesis pipeline, producing a library of triazole-based analogues ready for the next step: virtual screening and docking against CZP and hCatL. Wait... ready? 🤔

## 4. Candidates prioritization

To reduce late-stage attrition and keep *wet-lab* efforts focused, we will filter the virtual products by **predicted human cytotoxicity** using again an **Ersilia** model. Then, we will keep candidates with low predicted cytotoxicity for the docking stage.

Here we use the EOS model `eos21q7`, which provides a probability-like score stored as `eos21q7_dili_probability` (*higher values* = *greater risk*).

```python title="workshop.py"
la_workshop_cs.apply_ersilia_model_on_table("reaction_set_1","eos21q7")
```

Query the project database from bash using `sqlite3` to explore the distribution and set realistic thresholds:

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos21q7_dili_probability), MAX(eos21q7_dili_probability) FROM reaction_set_1 WHERE eos21q7_dili_probability IS NOT NULL;
```

*Outputs*: `0.28|0.73`

Here we select the 5000 compounds with the lowest predicted cytotoxicity and store them in a new table. Only those compounds will be submitted to molecular docking.

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db \
"CREATE TABLE candidates_for_docking AS
 SELECT *
 FROM reaction_set_1
 WHERE eos21q7_dili_probability IS NOT NULL
 ORDER BY eos21q7_dili_probability ASC
 LIMIT 5000;"
```

🙌 You now have a **price-aware, drug-like, synthetically feasible, and low-cytotoxicity** candidate set; effectively defining a **rational and efficient chemical space** for *in silico* exploration. These compounds are ready for the **virtual screening stage** against **CZP**, followed by **counter-screening** versus **hCatL** to assess selectivity.


## 5. vHTS

### 5.1. ⚙️ Molecular docking of prioritized candidates

With the chemical space now refined and filtered, we are ready to perform *molecular docking* of the 5K safest and most promising compounds against **CZP**; and later, for selectivity assessment, against **hCatL**.

First, prepare all ligands in the `candidates_for_docking` table.  

```python title="workshop.py"
la_workshop_cs.generate_mols_in_table("candidates_for_docking")
```

Instantiate a MolDock object to enable access to docking-related functions:

```python title="workshop.py"
la_workshop_moldock = md.MolDock(la_workshop)
```

Set the docking campaign using your prepared receptor and validated parameter sets. In this example, we use receptor model 3 (the refined CZP structure) and docking parameters 2.

```python title="workshop.py"
la_workshop_moldock.dock_table("candidates_for_docking",
                               id_receptor_model=3,
                               id_docking_params=2)
```

Once the assay is created, TidyScreen automatically generates multiple executable scripts, each containing 1000 docking runs (in total, 5000 compounds = 5 scripts).

You can launch them directly from the terminal:

```bash
for i in $(seq 1 11); do ./docking_execution_${i}.sh ; done
```

:::warning[Important]
Docking execution requires access to a local GPU board or a compatible GPU computing environment. If no GPU resources are available, these scripts cannot be executed, and docking should instead be run on a properly configured workstation or HPC cluster.
:::

### 5.2. 📊 Docking results analysis


### 5.3. 🎯 Selection of candidates for synthesis


