---
title: 5. Chemical space design and vHTS workflow
---

# Chemical Space Design and vHTS Workflow

<div style={{ textAlign: "justify" }}>

In this final section, we will extend the previously validated docking methodology to generate and evaluate **new triazole-based Targeted Covalent Inhibitors (TCIs)** of Cruzipain (CZP).

Following the successful validation of the docking workflow using a training set of triazole-containing derivatives featuring a methyleneamine linker as bioisosteric replacements of peptide bonds and the well-known efficacy of the vinylsulfone-containing inhibitor **K777**, our goal is to design a new series of molecules that **combine these three structural elements**.

Given the absence of reported analogues combining these motifs, we will perform a massive virtual library generation and screening campaign aimed at identifying candidates with **high affinity for CZP** and **minimal affinity for hCatL**, while ensuring **synthetic feasibility** for subsequent *wet-lab* validation.

---

## Overview of the Protocol

The overall workflow can be divided into five major steps:

**5.1- Extraction of building blocks**  
   Retrieve potential reactants from the *eMolecules* chemical database.

**5.2- Chemical space reduction**  
   Filter the combinatorial chemical space based on:
   - *Drug-likeness* criteria 
   - *Synthetic feasibility*
   - *Affordability*, using models provided through **Ersilia Hub**.

**5.3- In silico synthesis of derivatives**  
   Combine the selected reactants through feasible *in silico* synthetic routes to generate the final compounds, replicating realistic reactions that can be performed in the wet lab.

**5.4- Toxicity filtering**  
   Further reduce the library by eliminating potentially hepatotoxic molecules, applying an **Ersilia Hub DILI prediction model**.

**5.5- Docking and analysis**  
   Perform docking-based *vHTS* on the safest and most promising analogues (approximately **5000 compounds**) and analyze their predicted affinity and selectivity profiles toward CZP and hCatL.

---

🎯 This will allow prioritization of candidates for *in vitro* and *in vivo* validation, completing the *dry*-to-*wet lab* pipeline proposed in this tutorial.

---

## 5.1- **Extraction of building blocks**  


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

### 5.1.1- 🗂️ Import eMolecules database

Import the chemical database file from **eMolecules** (usually stored in `chemspace/raw_data/`) into the local **TidyScreen** project database (`chemspace/processed_data/chemspace.db`).

This step converts the original `.csv` file into a structured SQLite database that can be efficiently queried by TidyScreen for compound selection and filtering.

```python title="workshop.py"
emolecules_database_file = "/PATH/TO/FILE/emolecules.csv" #Set your own path!
la_workshop_cs.input_csv(emolecules_database_file)
```

### 5.1.2- 🧩 Chemical scaffolds filtering

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

### 5.1.3- 📸 Inspection of a sample of filtered building blocks

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

## 5.2- Building blocks prioritization

In order to narrow down the huge number of potential reactants, we will compute a set of standard **drug-like properties** and apply **Ersilia models** for each curated subset and then **prioritize** entries by applying simple, reproducible property filters.


### 5.2.1- ⚗️ Subsetting by drug-like properties (RDKit)

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

### 5.2.2- 💸 Price-based prioritization with Ersilia

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

## 5.3- Combinatorial virtual synthesis

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

### 5.3.1- 🧠 Defining SMARTS reactions

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

### 5.3.2- 🧪 Creating and applying a reaction workflow

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

## 5.4- Candidates prioritization

To reduce late-stage attrition and keep *wet-lab* efforts focused, we will filter the virtual products by **predicted human Drug-Induced Liver Injury (DILI)** using again an **Ersilia** model. Then, we will keep candidates with low predicted totoxicity for the docking stage.

Here we use the EOS model `eos21q7`, which provides a probability-like score stored as `eos21q7_dili_probability` (*higher values* = *greater risk*).

```python title="workshop.py"
la_workshop_cs.apply_ersilia_model_on_table("reaction_set_1","eos21q7")
```

Query the project database from bash using `sqlite3` to explore the distribution and set realistic thresholds:

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db "SELECT MIN(eos21q7_dili_probability), MAX(eos21q7_dili_probability) FROM reaction_set_1 WHERE eos21q7_dili_probability IS NOT NULL;
```

*Outputs*: `0.28|0.73`

Here we select the 5000 compounds with the lowest predicted toxicity and store them in a new table. Only those compounds will be submitted to molecular docking.

```bash
sqlite3 $PATH/TO/PROJECT/chemspace/processed_data/chemspace.db \
"CREATE TABLE candidates_for_docking AS
 SELECT *
 FROM reaction_set_1
 WHERE eos21q7_dili_probability IS NOT NULL
 ORDER BY eos21q7_dili_probability ASC
 LIMIT 5000;"
```

According to the same EOS DILI prediction model, the reference inhibitor **K777** exhibits a DILI probability of `0.56`, indicating a relatively high risk of hepatotoxicity.

In contrast, the 5000 compounds retained for docking fall within a much safer range (`0.28-0.38`). This confirms that, beyond synthetic feasibility and drug-likeness, our chemical space is biased toward lower predicted human toxicity, effectively prioritizing compounds with improved safety margins compared to known reference inhibitors.

🙌 You now have a **price-aware, drug-like, synthetically feasible, and low-toxicity** candidate set; effectively defining a **rational and efficient chemical space** for *in silico* exploration. These compounds are ready for the **virtual screening stage** against **CZP**, followed by **counter-screening** versus **hCatL** to assess selectivity.


## 5.5- vHTS

### 5.5.1- ⚙️ Molecular docking of prioritized candidates

With the chemical space now refined and filtered, we are ready to perform *molecular docking* of the 5K safest and most promising compounds against **CZP**; and later, for selectivity assessment, against **hCatL**.

First, prepare all ligands in the `candidates_for_docking` table.  

```python title="workshop.py"
la_workshop_cs.generate_mols_in_table("candidates_for_docking", conf_rank=25, timeout=300, pdbqt_method="meeko", charge_method="gas")
```

Instantiate a MolDock object to enable access to docking-related functions (If not previusly done in the workflow script):

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

### 5.5.2- 📊 Docking results analysis

After completing the large-scale docking run (`assay_6`), we can post-process the results and select candidates with the lowest binding energies.
As discussed in the previous tutorials (when working with K777 and the training set), the docking score proved to be the most consistent criterion for identifying bioactive poses. Therefore, we will prioritize it as our primary selection metric here.

```python title="workshop.py"
# Analyze the run and extract only the best pose per ligand (lowest-energy)
la_workshop_docking_analysis.process_docking_assay(
    assay_id=6,
    max_poses=1,
    extract_poses=1
)
```

Use SQL to subset the **Top 100 compounds by docking score**:

```bash
sqlite3 -header -table docking/docking_assays/assay_6/assay_6.db \
"SELECT LigName, sub_pose, docking_score
 FROM Results
 ORDER BY docking_score ASC
 LIMIT 100;"
```

<details>
<summary> <b>Output based on docking score</b></summary>
<div style={{ fontSize: '90%' }}>
| **LigName**                     | **sub_pose**                    | **docking_score** |
|---------------------------------|---------------------------------|------------------:|
| CYWLRULGWUKHND-QDRCWIAMSA-O     | CYWLRULGWUKHND-QDRCWIAMSA-O_1   | -8.15             |
| AABZCBAGAYRJDZ-ROLKXXGCSA-N     | AABZCBAGAYRJDZ-ROLKXXGCSA-N_1   | -8.00             |
| MBTIUSBPXHCMMC-GHDUTYINSA-O     | MBTIUSBPXHCMMC-GHDUTYINSA-O_1   | -7.99             |
| ULUDSBKRRIDGPN-LQTKSSSVSA-O     | ULUDSBKRRIDGPN-LQTKSSSVSA-O_1   | -7.98             |
| XPIFPVCHVHNQHJ-RILXKHBHSA-O     | XPIFPVCHVHNQHJ-RILXKHBHSA-O_1   | -7.94             |
| YKLJSRDWNOQPQO-OQBDGEOBSA-N     | YKLJSRDWNOQPQO-OQBDGEOBSA-N_1   | -7.94             |
| VBFTWLUZOOZBPB-RJCWRSEOSA-O     | VBFTWLUZOOZBPB-RJCWRSEOSA-O_1   | -7.88             |
| RYFVONWDPDPHHB-ZGCOHECXSA-O     | RYFVONWDPDPHHB-ZGCOHECXSA-O_1   | -7.81             |
| NFHBDVKWICVFDQ-APDVVUBJSA-O     | NFHBDVKWICVFDQ-APDVVUBJSA-O_1   | -7.81             |
| WERVHOCRXRNYBH-REISJJRISA-O     | WERVHOCRXRNYBH-REISJJRISA-O_1   | -7.74             |
| PFZHGYDOQTZDPQ-HUPUPJIPSA-N     | PFZHGYDOQTZDPQ-HUPUPJIPSA-N_1   | -7.72             |
| MMKPYJXQEXUYDH-OIBPPNFKSA-O     | MMKPYJXQEXUYDH-OIBPPNFKSA-O_1   | -7.71             |
| LLSCRUHCEVGKSQ-XKHREPBBSA-O     | LLSCRUHCEVGKSQ-XKHREPBBSA-O_1   | -7.70             |
| DOCNEKCJXXPFEZ-NVOBCCQGSA-O     | DOCNEKCJXXPFEZ-NVOBCCQGSA-O_1   | -7.68             |
| DYRCPBZRUZAXIX-JKPVDXOUSA-O     | DYRCPBZRUZAXIX-JKPVDXOUSA-O_1   | -7.67             |
| OMQWCNCSDAOYQS-OJPUXPQHSA-O     | OMQWCNCSDAOYQS-OJPUXPQHSA-O_1   | -7.65             |
| GEGJOGATSOGXEP-CCVDMMCGSA-O     | GEGJOGATSOGXEP-CCVDMMCGSA-O_1   | -7.65             |
| CPECBCLWZCCFIR-DWLBMTHHSA-O     | CPECBCLWZCCFIR-DWLBMTHHSA-O_1   | -7.61             |
| JSATUSIYCVGCOF-GXCCJSRQSA-O     | JSATUSIYCVGCOF-GXCCJSRQSA-O_1   | -7.61             |
| MSWGHHCHSHNSRS-JMRBFUEQSA-O     | MSWGHHCHSHNSRS-JMRBFUEQSA-O_1   | -7.60             |
| UZLBPXKFEJTBJW-QAYYRARUSA-O     | UZLBPXKFEJTBJW-QAYYRARUSA-O_1   | -7.60             |
| YQBRHEDHYQJPGY-NVOBCCQGSA-O     | YQBRHEDHYQJPGY-NVOBCCQGSA-O_1   | -7.60             |
| KGFVBUCFAKEIIE-HXCDMFHQSA-O     | KGFVBUCFAKEIIE-HXCDMFHQSA-O_1   | -7.59             |
| NQZUNDAQXMGKSL-XAWYSJHESA-O     | NQZUNDAQXMGKSL-XAWYSJHESA-O_1   | -7.59             |
| QMUPQKJUACHEQD-OUCRSZIVSA-N     | QMUPQKJUACHEQD-OUCRSZIVSA-N_1   | -7.58             |
| XBSZIXSVKWJTSD-OGPZQOLISA-O     | XBSZIXSVKWJTSD-OGPZQOLISA-O_1   | -7.58             |
| ZLAUXDJLABCRDU-OURVLYCJSA-O     | ZLAUXDJLABCRDU-OURVLYCJSA-O_1   | -7.57             |
| UGDNQIHMIAUMOF-LEEQTAOHSA-O     | UGDNQIHMIAUMOF-LEEQTAOHSA-O_1   | -7.57             |
| SRXXWAWPFVBRME-HDKPYIOWSA-O     | SRXXWAWPFVBRME-HDKPYIOWSA-O_1   | -7.57             |
| YWHRELYUHVGYPT-GLLYUYBOSA-O     | YWHRELYUHVGYPT-GLLYUYBOSA-O_1   | -7.55             |
|...|...|...|
| YGZRUQQIZMBSDD-XXEOEALZSA-O     | YGZRUQQIZMBSDD-XXEOEALZSA-O_1   | -7.32             |
| DDCLCWAKUGQATH-KXZHXVJYSA-O     | DDCLCWAKUGQATH-KXZHXVJYSA-O_1   | -7.32             |
| OXUZYRNVCOIXKL-ZGJBHGAESA-O     | OXUZYRNVCOIXKL-ZGJBHGAESA-O_1   | -7.31             |
| PRKCCYOLKCDSCQ-LFCSLWFPSA-O     | PRKCCYOLKCDSCQ-LFCSLWFPSA-O_1   | -7.30             |
| RBFWAFPREHYGNU-ZEMQAGKDSA-O     | RBFWAFPREHYGNU-ZEMQAGKDSA-O_1   | -7.30             |
| DYYZDMIHWINFIF-PXDOLEPXSA-O     | DYYZDMIHWINFIF-PXDOLEPXSA-O_1   | -7.30             |
| CCNSRCIALOOJHQ-IYNFVXDOSA-O     | CCNSRCIALOOJHQ-IYNFVXDOSA-O_1   | -7.29             |
| FPMTZBSRRXATEP-KQFRVKRDSA-O     | FPMTZBSRRXATEP-KQFRVKRDSA-O_1   | -7.29             |
| IJTFVGGDMZVVCY-NEJGVSAFSA-O     | IJTFVGGDMZVVCY-NEJGVSAFSA-O_1   | -7.29             |
| JXDRXFUSYAZPDG-SGGLPICOSA-O     | JXDRXFUSYAZPDG-SGGLPICOSA-O_1   | -7.28             |
| ABYYJZBHRQSKOH-SFBQEAGJSA-O     | ABYYJZBHRQSKOH-SFBQEAGJSA-O_1   | -7.28             |
| YVYRLDPJKNDTKF-FRNNXQESSA-O     | YVYRLDPJKNDTKF-FRNNXQESSA-O_1   | -7.28             |
| TZNGNDMWLRPWHC-AKVYSXHKSA-O     | TZNGNDMWLRPWHC-AKVYSXHKSA-O_1   | -7.28             |
| MBRMELWIASNFDQ-BTXCHTAFSA-O     | MBRMELWIASNFDQ-BTXCHTAFSA-O_1   | -7.27             |
| OIMCTJQEECQBCO-ZGJDELOYSA-O     | OIMCTJQEECQBCO-ZGJDELOYSA-O_1   | -7.27             |
</div>
</details>


Upon inspecting the output, you will notice that the docking scores are considerably more negative than those obtained for K777 bound to CZP (≈ -4.08 kcal.mol<sup>-1</sup>), suggesting that more stable encounter complexes may have been identified in this large-scale screening.

To reinforce the validity of these findings, we will now redock the top 100 potential inhibitors to verify reproducibility and ensure that the observed ranking trends remain consistent across independent docking runs.

Create a table named `top_100_CZP_binders` to store the selected compounds for further analysis:

```bash
sqlite3 chemspace/processed_data/chemspace.db \
"ATTACH DATABASE 'docking/docking_assays/assay_6/assay_6.db' AS db2;
 CREATE TABLE top_100_CZP_binders AS
 SELECT * FROM candidates_for_docking
 WHERE inchi_key IN (
     SELECT LigName FROM db2.Results
     ORDER BY docking_score ASC
     LIMIT 100
 );
 DETACH DATABASE db2;"
```

Create a new docking assay for `top_100_CZP_binders` using the refined CZP receptor and the validated docking parameters to confirm reproducibility:

```python title="workshop.py"
la_workshop_moldock.dock_table(
    "top_100_CZP_binders",
    id_receptor_model=3,
    id_docking_params=2
)
```
Analyze the docking results for `assay_7`: 

```python title="workshop.py"
# Extract only the lowest-energy pose per ligand
la_workshop_docking_analysis.process_docking_assay(
    assay_id=7,
    max_poses=1,
    extract_poses=1
)  
```

Finally, compute the MMGBSA energy fingerprints for the newly generated complexes to evaluate their thermodynamic stability:

```python title="workshop.py"
# Provide receptor field indices when prompted
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(
    assay_id=7,
    mmgbsa=1,
    prolif=0
)  
```

Although the detailed results are not displayed here, the redocking step confirmed excellent reproducibility across independent assays. The calculated MMGBSA values were highly consistent with the previous run, showing only minor deviations (~0.2 kcal.mol<sup>-1</sup>) for a few ligands, well within the expected precision range of the method.

---- 

Having validated the workflow and confirmed the robustness of our docking parameters for CZP, the next logical step is to apply the same protocol against the refined hCatL receptor. This comparative analysis will allow us to **identify ligands exhibiting preferential binding to CZP over hCatL**, thus highlighting selectivity determinants relevant for the design of next-generation antichagasic agents.


Create a new docking assay for `top_100_CZP_binders` using the refined hCatL receptor and the validated docking parameters to confirm reproducibility:

```python title="workshop.py"
la_workshop_moldock.dock_table(
    "top_100_CZP_binders",
    id_receptor_model=4,
    id_docking_params=2
)
```
Analyze the docking results for `assay_8`: 

```python title="workshop.py"
# Extract only the lowest-energy pose per ligand
la_workshop_docking_analysis.process_docking_assay(
    assay_id=8,
    max_poses=1,
    extract_poses=1
)  
```

Finally, compute the MMGBSA energy fingerprints for the newly generated complexes to evaluate their thermodynamic stability:

```python title="workshop.py"
# Provide receptor field indices when prompted
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(
    assay_id=8,
    mmgbsa=1,
    prolif=0
)  
```

Now that both docking and MMGBSA analyses have been completed for CZP and hCatL, we can directly compare their results to identify ligands that preferentially bind to CZP.

### 5.5.3- 🎯 Selection of candidates for synthesis

The following SQL query joins both assay databases to retrieve, for each compound, the docking scores and MMGBSA total energies obtained against both receptors.
By computing the difference between the values (`diff_dock_score` and `diff_d_g_tot`), we can easily pinpoint candidates with improved binding affinity for CZP relative to hCatL.

:::warning[Important]
Please make sure to **replace `/PATH/TO/PROJECT/`** with your **own local or cluster directory path** before running the command below. If not updated, the SQL query will fail because the database cannot be located.
:::

``` bash
sqlite3 -header -table docking/docking_assays/assay_7/assay_7.db \
"ATTACH DATABASE '/PATH/TO/PROJECT/docking/docking_assays/assay_8/assay_8.db' AS db2;
 SELECT t1.sub_pose,
        t1.docking_score AS CZP_ds,
        t2.delta_g_total AS CZP_d_g_tot,
        db2.Results.docking_score AS hCatL_ds,
        db2.mmgbsa_fingerprints.delta_g_total AS hCatL_d_g_tot,
        ROUND((t1.docking_score - db2.Results.docking_score),2) AS diff_dock_score,
        ROUND((t2.delta_g_total - db2.mmgbsa_fingerprints.delta_g_total),2) AS diff_d_g_tot
 FROM Results t1
 JOIN mmgbsa_fingerprints t2 ON t1.LigName = t2.LigName
 JOIN db2.Results ON t1.LigName = db2.Results.LigName
 JOIN db2.mmgbsa_fingerprints ON t1.LigName = db2.mmgbsa_fingerprints.LigName
 ORDER BY diff_dock_score ASC
 LIMIT 10;
 DETACH DATABASE db2;"
``` 

----

<details> 
<summary> <b>Output ranked by docking score difference</b></summary>
| **sub_pose**                    | **CZP_ds** | **CZP_d_g_tot** | **hCatL_ds** | **hCatL_d_g_tot** | **diff_dock_score** | **diff_d_g_tot** |
|---------------------------------|------------:|----------------:|--------------:|------------------:|--------------------:|----------------|
| CUFNDGFOXHNKKU-MGUYBSQASA-O_1  | -7.24       | -51.6476        | -5.85         | -39.9543          | -1.39               | -11.69          |
| CYWLRULGWUKHND-QDRCWIAMSA-O_1  | -8.09       | -34.9057        | -6.73         | -29.0774          | -1.36               | -5.83           |
| TZCGFRFEIJJRGN-VKFROBJISA-N_1  | -7.44       | -49.0613        | -6.20         | -39.6637          | -1.24               | -9.40           |
| GEGJOGATSOGXEP-CCVDMMCGSA-O_1  | -7.66       | -39.9591        | -6.48         | -38.7983          | -1.18               | -1.16           |
| IHRQBMNUDKYYTK-HUAVMZDFSA-N_1  | -7.28       | -51.4805        | -6.24         | -45.7701          | -1.04               | -5.71           |
| DDCLCWAKUGQATH-KXZHXVJYSA-O_1  | -7.49       | -42.1068        | -6.47         | -41.6397          | -1.02               | -0.47           |
| KHUSNQDEVMQYSO-GEKJTBGDSA-O_1  | -7.51       | -35.6408        | -6.59         | -36.5528          | -0.92               | 0.91            |
| AABZCBAGAYRJDZ-ROLKXXGCSA-N_1  | -8.19       | -38.4541        | -7.29         | -44.8960          | -0.90               | 6.44            |
| MSWGHHCHSHNSRS-JMRBFUEQSA-O_1  | -7.60       | -48.8381        | -6.70         | -47.1861          | -0.90               | -1.65           |
| JXDRXFUSYAZPDG-SGGLPICOSA-O_1  | -7.29       | -49.5987        | -6.45         | -52.0911          | -0.84               | 2.49            |
</details>

Upon inspection of the results, we will focus on the three top-scoring ligands according to their docking scores, as they exhibit the most promising binding affinities toward CZP.

<div style={{ fontSize: '90%' }}>
| **sub_pose**                    | **CZP_ds** | **CZP_d_g_tot** | **hCatL_ds** | **hCatL_d_g_tot** | **diff_dock_score** | **diff_d_g_tot** |
|---------------------------------|------------:|----------------:|--------------:|------------------:|--------------------:|----------------|
| CUFNDGFOXHNKKU-MGUYBSQASA-O_1  | -7.24       | -51.6476        | -5.85         | -39.9543          | -1.39               | -11.69          |
| CYWLRULGWUKHND-QDRCWIAMSA-O_1  | -8.09       | -34.9057        | -6.73         | -29.0774          | -1.36               | -5.83           |
| TZCGFRFEIJJRGN-VKFROBJISA-N_1  | -7.44       | -49.0613        | -6.20         | -39.6637          | -1.24               | -9.40           |
| K777  | -5.42       | -50.2336        | -5.17         | -43.6419          | -0.25               | -6.59           |
</div>

1. **CUFNDGFOXHNKKU-MGUYBSQASA-O** 

In the first case, the compound **CUFNDGFOXHNKKU-MGUYBSQASA-O** stands out due to its large affinity gap between CZP and hCatL, both in docking score and MMGBSA; the latter exceeding even the 6.6 kcal.mol<sup>-1</sup> difference observed for K777.

However, upon visual inspection, we identified a critical issue: this compound was incorrectly generated during the in silico synthesis. 

#################### FIGURE 2D ####################

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/CUFNDGFOXHNKKU-MGUYBSQASA-O.png" alt="Description of image" width="400"/>
  <figcaption>Chemical structure of candidate 'CUFNDGFOXHNKKU-MGUYBSQASA-O'. </figcaption>
  </p>
</figure>
---

Specifically, the nitrile group of the amine derivative (R<sup>3</sup>) was spuriously used as the reacting atom during the A3 coupling step. This misreaction originated from an overly permissive SMARTS pattern of the N in the reaction definition:

`la_workshop_cs.add_smarts_reaction("[N:1].[CX3H1:2](=[O:3])>>[NX4+:1][C@H:2][C:4]#[C:5]", "A3 coupling")`

This compound was intentionally retained in the dataset for educational purposes, to emphasize the importance of thoroughly validating synthetic rules and performing sanity checks prior to large-scale docking campaigns, avoiding unnecessary computational costs and misinterpretations.

Having said that, most of the remaining ligands were correctly synthesized, allowing for a reliable comparative analysis of their binding behaviors.

2. **CYWLRULGWUKHND-QDRCWIAMSA-O** 

This compound exhibits a significant difference in docking score between CZP and hCatL, favoring CZP binding; in contrast to K777, for which the difference was reversed.
However, the ΔG (MMGBSA) difference is smaller, and the absolute binding energy is also less favorable.

Nevertheless, visual inspection suggests this is likely the ideal outcome:
* For CZP, the predicted pose corresponds to a bioactive conformation, where the warhead is properly oriented toward the catalytic cysteine, and the triazole ring effectively mimics the amide bond through bioisosteric replacement.
* Conversely, in hCatL, the predicted pose is not geometrically compatible with covalent inhibition, indicating that this compound could display promising selectivity toward CZP and thus represents a synthetically feasible candidate for further exploration.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/CYWLRULGWUKHND-QDRCWIAMSA-O.png" alt="Description of image" width="400"/>
  <figcaption>Chemical structure of candidate 'CYWLRULGWUKHND-QDRCWIAMSA-O'. </figcaption>
  </p>
</figure>
---


<div style={{display: 'flex', justifyContent: 'space-around', alignItems: 'center', gap: '20px'}}>
  <figure style={{textAlign: 'center', flex: 1}}>
    <img src="/TidyScreen_v2_docs_new/img/CYWLRULGWUKHND-QDRCWIAMSA-O-CZP.png" alt="CZP fingerprint" width="285"/>
    <figcaption><b>CYWLRULGWUKHND-QDRCWIAMSA-O_1 bound to CZP</b></figcaption>
  </figure>
  <figure style={{textAlign: 'center', flex: 1}}>
    <img src="/TidyScreen_v2_docs_new/img/CYWLRULGWUKHND-QDRCWIAMSA-O-CatL.png" alt="hCatL fingerprint" width="450"/>
    <figcaption><b>CYWLRULGWUKHND-QDRCWIAMSA-O_1 bound to hCatL</b></figcaption>
  </figure>
</div>



3. **TZCGFRFEIJJRGN-VKFROBJISA-N**

In this case, all indicators are favorable: both the docking score and MMGBSA ΔG values are improved for CZP relative to hCatL, and the absolute ΔG values are comparable to those obtained for K777.
Visual inspection further supports this observation: the compound adopts bioactive-like poses in both targets, properly aligned for potential covalent inhibition.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/TZCGFRFEIJJRGN-VKFROBJISA-N.png" alt="Description of image" width="300"/>
  <figcaption>Chemical structure of candidate 'TZCGFRFEIJJRGN-VKFROBJISA-N'. </figcaption>
  </p>
</figure>
---


<div style={{display: 'flex', justifyContent: 'space-around', alignItems: 'center', gap: '20px'}}>
  <figure style={{textAlign: 'center', flex: 1}}>
    <img src="/TidyScreen_v2_docs_new/img/TZCGFRFEIJJRGN-VKFROBJISA-N-CZP.png" alt="CZP fingerprint" width="400"/>
    <figcaption><b>TZCGFRFEIJJRGN-VKFROBJISA-N_1 bound to CZP</b></figcaption>
  </figure>
  <figure style={{textAlign: 'center', flex: 1}}>
    <img src="/TidyScreen_v2_docs_new/img/TZCGFRFEIJJRGN-VKFROBJISA-N-CatL.png" alt="hCatL fingerprint" width="400"/>
    <figcaption><b>TZCGFRFEIJJRGN-VKFROBJISA-N_1 bound to hCatL</b></figcaption>
  </figure>
</div>


While this may suggest moderate selectivity between CZP and hCatL, it also indicates a favorable binding profile that could translate into reduced off-target toxicity; particularly when compared to K777, as preliminarily estimated by Ersilia’s toxicity predictor. 

## 🧩 5.6- Summary of the Chemical Space Exploration

Taken together, this small-scale benchmark exemplifies the diversity of outcomes that can emerge from rational in silico design workflows.
Among the top-ranked ligands, we identified three representative cases:

- a **synthetic artefact** (*CUFNDGFOXHNKKU-MGUYBSQASA-O*), which highlights the need for **rigorous validation** of virtual reaction schemes and filters before large-scale screening;

- a **balanced binder** (*TZCGFRFEIJJRGN-VKFROBJISA-N*), displaying **comparable affinities** for both targets but potentially improved safety relative to K777; and

- a **selective inhibitor candidate** (*CYWLRULGWUKHND-QDRCWIAMSA-O*), favoring CZP over hCatL and reproducing a bioactive-like orientation.

Altogether, these findings illustrate how in silico chemical space exploration can efficiently combine synthetic accessibility, predicted safety, and mechanistic insight to refine candidate selection before synthesis. 


<details> 
<summary> <b>💭 WHAT IF...? </b></summary>

So far, we have relied on manual inspection and energetic metrics (Docking score, ΔG<sub>MMGBSA</sub>) to evaluate candidate poses.
But what if we could go one step further, quantifying the per-residue interactions and learning from them automatically?

Each docking result already contains per-residue interaction fingerprints, which can be extracted directly from the `mmgbsa_fingerprints` table as shown below:

```bash
sqlite3 docking/docking_assays/assay_2/assay_2.db \
"SELECT writefile('fp_K777-CZP.csv', mmgbsa_csv_file) 
 FROM mmgbsa_fingerprints 
 WHERE sub_pose = 'RHJLQMVZXQKJKB-FPHSVDBKSA-N_1';"  # K777 bound to CZP

sqlite3 docking/docking_assays/assay_4/assay_4.db \
"SELECT writefile('fp_K777-hCatL.csv', mmgbsa_csv_file) 
 FROM mmgbsa_fingerprints 
 WHERE sub_pose = 'RHJLQMVZXQKJKB-FPHSVDBKSA-N_1';"  # K777 bound to hCatL

sqlite3 docking/docking_assays/assay_7/assay_7.db \
"SELECT writefile('fp_CYWLRULGWUKHND-QDRCWIAMSA-O_1-CZP.csv', mmgbsa_csv_file) 
 FROM mmgbsa_fingerprints 
 WHERE sub_pose = 'CYWLRULGWUKHND-QDRCWIAMSA-O_1';"  # Candidate bound to CZP

sqlite3 docking/docking_assays/assay_8/assay_8.db \
"SELECT writefile('fp_CYWLRULGWUKHND-QDRCWIAMSA-O_1-hCatL.csv', mmgbsa_csv_file) 
 FROM mmgbsa_fingerprints 
 WHERE sub_pose = 'CYWLRULGWUKHND-QDRCWIAMSA-O_1';"  # Candidate bound to hCatL

sqlite3 docking/docking_assays/assay_7/assay_7.db \
"SELECT writefile('fp_TZCGFRFEIJJRGN-VKFROBJISA-N_1-CZP.csv', mmgbsa_csv_file) 
 FROM mmgbsa_fingerprints 
 WHERE sub_pose = 'TZCGFRFEIJJRGN-VKFROBJISA-N_1';"  # Candidate bound to CZP

sqlite3 docking/docking_assays/assay_8/assay_8.db \
"SELECT writefile('fp_TZCGFRFEIJJRGN-VKFROBJISA-N_1-hCatL.csv', mmgbsa_csv_file) 
 FROM mmgbsa_fingerprints 
 WHERE sub_pose = 'TZCGFRFEIJJRGN-VKFROBJISA-N_1';"  # Candidate bound to hCatL
```
These CSV files encode per-residue contributions (vdW, electrostatics, solvation, etc.), enabling:
* automated filtering of poses that reproduce key binding interactions ("hotspots"),
* quantitative comparison across targets, and
* use as training data for machine learning classifiers, capable of distinguishing bioactive from non-bioactive binding modes.

In practice, this approach can drastically reduce the need for visual inspection while providing interpretable molecular descriptors for downstream AI-driven design pipelines.

🔍 **Comparative analysis of interaction fingerprints**

| CZP fingerprint | hCatL fingerprint |
|-----------------|-------------------|
| ![](/img/mmgbsa_group_total_by_residue_CZP.png) | ![](/img/mmgbsa_group_total_by_residue_hCatL.png) |

In this example, K777 exhibits a similar interaction pattern across both targets, as expected for a non-selective covalent inhibitor.

In contrast, the compound *CYWLRULGWUKHND-QDRCWIAMSA-O* shows a marked difference in interaction profiles:

* Gln19 (CZP) / Gln20 (hCatL), both residues belong to the oxyanion hole, essential for catalytic stabilization.
* HIP162 (CZP) / HIP164 (hCatL), the histidine partner of Cys25/Cys26 forming the catalytic dyad.

Had we applied an interaction-based filter requiring the satisfaction of these two key contacts, the analysis would have automatically identified *CYWLRULGWUKHND-QDRCWIAMSA-O* as a CZP-selective binder, while *TZCGFRFEIJJRGN-VKFROBJISA-N* would emerge as bioactive in both targets.


</details>



</div>