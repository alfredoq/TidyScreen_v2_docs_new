---
title: 4. Evaluation of a training set
---

Before performing large-scale virtual screening, it is essential to validate the docking workflow using a **well-characterized benchmark set**.  
This allows verifying that the selected docking parameters, scoring functions, and receptor models can **reproduce experimental trends** and correctly discriminate between active and less active compounds.

To date, no **triazole-based inhibitors bearing a vinylsulfone warhead** have been reported, except for one analog described by *Brak et al.* (2008).  
In their pioneering work and follow-up publications, the same group explored triazole-containing analogues with alternative reactive moieties, mainly **4-fluorophenoxymethyl ketones (4FPMKs)**, originally developed for **cathepsin** inhibition.

Their design strategy consisted in **replacing the scissile amide bond** of peptidic inhibitors with a **1,4-disubstituted 1,2,3-triazole** ring while systematically varying the **R<sup>1</sup>, R<sup>2</sup>, and R<sup>3</sup>** substituents to explore ligand–target complementarity. Subsequent optimization led to a family of **20 triazole-based 4FPMK derivatives**, collectively showing a wide **inhibitory range against Cruzipain (CZP)**, from **3 nM to 23 µM**.  Although their **SAR remains only partially explored**, these analogues demonstrated good trypanocidal activity, acceptable physicochemical profiles, and promising *in vivo* efficacy - despite safety concerns linked to CYP3A4 inhibition. One representative compound, **Ts-370**, was co-crystallized with CZP (PDB ID: **3IUT**), providing valuable structural insights into the binding mode of this series and serving as a solid reference for further *in silico* modeling efforts. 

For the purpose of this workshop, we will use this small but well-defined set as a **test set**. The goal is to evaluate whether our SBDD workflow can reproduce the **binding trends and relative affinities** observed experimentally, using the available crystallographic and biochemical data.

While best practice in computational drug design recommends **larger and chemically diverse training/test sets**, this reduced dataset provides a **realistic and time-efficient framework** for validating the docking setup and analyzing performance under controlled conditions during this workshop.

#### Test set

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/Training_Set_Workshop.png" alt="Description of image" width="700"/>
  <figcaption>  Chemical structures of the compounds considered as the <b>test set</b>. The reference inhibitor <b>K777</b>, its triazolized analogue (<b>Tz-K777</b>), and the unique triazole-based vinylsulfone derivative reported to date (<b>KB-38</b>), along with the series of 20 <b>4FPMK</b> derivatives originally synthesized by Jonathan A. Ellman (check Table below), are included. </figcaption>
  </p>
</figure>
---


| **ID** | **R<sup>1</sup>** | **R<sup>2</sup>** | **X** | **R<sup>3</sup>** | **IC50 (nM)** |
|:-------:|:-------|:-------|:------|:-------|:--------------:|
| **K777** | - | - | - | - | **2** |
| **Tz-K777** | - | - | - | - | **n.d.** |
| **KB-38** | - | - | - | - | **n.d.** |
| **Ts-3** | *n*-Bu | *c*-Pent | CH₂ | PzPm | **3** |
| **Ts-5** | *n*-Bu | *c*-Pent | CH₂ | BnThia | **5** |
| **Ts-25** | msp | *c*-Pent | CH₂ | BnThia | **25** |
| **Ts-27** | *n*-Bu | *c*-Pent | CH₂ | Quin | **27** |
| **Ts-37** | *n*-Bu | *i*-Pr | CH₂ | BnThia | **37** |
| **Ts-63** | mop | *i*-Pr | CH₂ | BnThia | **63** |
| **Ts-75** | *n*-Bu | *i*-Pr | CH₂ | PzPm | **75** |
| **Ts-124** | mop | *c*-Pent | CH₂ | BnThia | **124** |
| **Ts-125** | msp | *i*-Pr | CH₂ | BnThia | **125** |
| **Ts-340** | mop | *i*-Pr | CH₂ | Quin | **340** |
| **Ts-370** | *n*-Bu | *i*-Pr | CH₂ | Quin | **370** |
| **Ts-467** | msp | *i*-Pr | CH₂ | PzPm | **467** |
| **Ts-560** | msp | *i*-Pr | CH₂ | Quin | **560** |
| **Ts-790** | msp | *c*-Pent | CH₂ | BnFur | **790** |
| **Ts-880** | msp | *i*-Pr | CH₂ | BnFur | **880** |
| **Ts-1100** | *n*-Bu | *i*-Pr | CH₂ | Py | **1100** |
| **Ts-1350** | *n*-Bu | *i*-Pr | C=O | Quin | **1350** |
| **Ts-2820** | mop | *c*-Pent | CH₂ | PzPm | **2820** |
| **Ts-3590** | *n*-Bu | *i*-Pr | C=O | BnThia | **3590** |
| **Ts-23000** | *n*-Bu | *i*-Pr | C=O | Py | **23000** |

> **Note:** IC50 values correspond to inhibition of *Cruzipain* (CZP) and are expressed in nanomolar (nM). Compound IDs reflect their approximate potency order.


To streamline this stage, we provide a precompiled CSV file containing the **SMILES**, **compound IDs**, and **inhibitory potency labels** for all molecules included in the test set.

This file can be directly imported into your local **TidyScreen** project and integrated into the **chemspace database** (`chemspace.db`) for subsequent docking validation.

---

<details>
<summary>📂 <b>Download: TS_Workshop.csv</b></summary>

[➡️ Click here to download the file](/files/TS_Workshop.csv)

</details>

After downloading, simply run the following command inside your Python environment to import the dataset into the chemspace module of your project:

```python title="workshop.py"
# Please check if the chemspace object is active! --> la_workshop_cs = chemspace.ChemSpace(la_workshop)
la_workshop_cs.input_csv("/PATH/TO/TS_Workshop.csv")
```

All entries will be automatically stored in your `chemspace.db` file under `/chemspace/processed_data/`, ready for use in the validation stage.

Depict the training set molecules to discuss SAR features

```python
la_workshop_cs.depict_ligand_table("TS_Workshop")
```

Now that the training/test dataset has been imported into your project, we will run a controlled **docking validation**  to ensure that the parameters previously optimized for CZP (and later to be used in large-scale screening) can reproduce reliable and consistent binding trends for known inhibitors.

So that, prepare all ligands in the `TS_Workshop` table, and convert them into docking-compatible formats.

```python title="workshop.py"
la_workshop_cs.generate_mols_in_table("TS_Workshop", conf_rank=25, timeout=300, pdbqt_method="meeko", charge_method="gas") # Note the parameters to generate conformers
```

Instantiate a MolDock object to activate all docking-related functions within your current project, and set up the docking campaign using the previously validated CZP receptor model (ID = 3) and the optimized docking parameters (ID = 2).

```python title="workshop.py"
la_workshop_moldock = md.MolDock(la_workshop)

la_workshop_moldock.dock_table("TS_Workshop",
                               id_receptor_model=3,
                               id_docking_params=2)
```

Once the docking assay is created, TidyScreen automatically generates an execution script.

To execute the job, you may run from your terminal:

```bash
./docking_execution_1.sh
```

:::warning[Important]
Docking execution requires access to a local GPU board or a compatible GPU computing environment. If no GPU resources are available, these scripts cannot be executed, and docking should instead be run on a properly configured workstation or HPC cluster.
:::

###### Analysis

Analyze the docking results for assay 5

```python title="workshop.py"
la_workshop_docking_analysis.process_docking_assay(assay_id=5,max_poses=5, extract_poses=1) # Notice that we extract only 5 poses per molecule and save them a pdb file
```

Compute the MMPBSA fingerprints for assay 5 
```python title="workshop.py"
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(assay_id=5,mmgbsa=1,prolif=0) # Input receptor field 
```

Use SQL statement to extract the best pose per ligand based on MMGBSA total energy

```bash
sqlite3 -header -table docking/docking_assays/assay_5/assay_5.db "SELECT LigName, sub_pose, delta_g_total FROM (SELECT LigName, sub_pose, delta_g_total, ROW_NUMBER() OVER (PARTITION BY LigName ORDER BY delta_g_total ASC) as rn FROM mmgbsa_fingerprints) ranked WHERE rn = 1;"

################
### Outputs ####
################
# +-----------------------------+-------------------------------+---------------+
# |           LigName           |           sub_pose            | delta_g_total |
# +-----------------------------+-------------------------------+---------------+
# | BFNJAXKMPJHXHW-DCFHFQCYSA-O | BFNJAXKMPJHXHW-DCFHFQCYSA-O_4 | -51.3374      | # NON-Bioactive conformation
# | CLEDCFNOFQFOGX-KMRXNPHXSA-O | CLEDCFNOFQFOGX-KMRXNPHXSA-O_2 | -34.4763      | # Bioactive conformation
# | CROWHGLIUOMOOH-SIBVEZHUSA-O | CROWHGLIUOMOOH-SIBVEZHUSA-O_2 | -34.6089      | # Bioactive conformation
# | FSMFTHPDNIOSTO-FNZWTVRRSA-O | FSMFTHPDNIOSTO-FNZWTVRRSA-O_2 | -45.7016      | # NON-Bioactive conformation
# | FTWLHZZOIPWOOS-IADCTJSHSA-O | FTWLHZZOIPWOOS-IADCTJSHSA-O_5 | -53.6438      | # NON-Bioactive conformation
# | GFKPWKCDJBRXQK-DCFHFQCYSA-O | GFKPWKCDJBRXQK-DCFHFQCYSA-O_4 | -51.5091      | # NON-Bioactive conformation
# | GJWDUYBKAOXWIS-ZTOMLWHTSA-O | GJWDUYBKAOXWIS-ZTOMLWHTSA-O_1 | -50.1655      | # Bioactive conformation
# | LCNNMCFKBKHPHV-DCFHFQCYSA-O | LCNNMCFKBKHPHV-DCFHFQCYSA-O_4 | -36.9908      | # NON-Bioactive conformation
# | NRDAGRURPLJETJ-KMRXNPHXSA-O | NRDAGRURPLJETJ-KMRXNPHXSA-O_4 | -52.3813      | # NON-Bioactive conformation
# | NXGQNXSIIMHEHD-NTFOOJQESA-O | NXGQNXSIIMHEHD-NTFOOJQESA-O_3 | -47.5815      | # NON-Bioactive conformation
# | OQDQGSSUGQXEKZ-JHOBJCJYSA-O | OQDQGSSUGQXEKZ-JHOBJCJYSA-O_4 | -52.5189      | # NON-Bioactive conformation
# | PWQIDUKXKZPEMX-AUDOKZQXSA-O | PWQIDUKXKZPEMX-AUDOKZQXSA-O_4 | -38.6934      | # NON-Bioactive conformation
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1 | -50.2251      | # Bioactive conformation
# | RSSMYFQRFUEXTN-ZTOMLWHTSA-O | RSSMYFQRFUEXTN-ZTOMLWHTSA-O_4 | -41.0921      | # NON-Bioactive conformation
# | STUFFQIGUGGQAV-KMRXNPHXSA-O | STUFFQIGUGGQAV-KMRXNPHXSA-O_5 | -53.2883      | # NON-Bioactive conformation
# | VCCYIWMCXDRTJQ-MMTVBGGISA-N | VCCYIWMCXDRTJQ-MMTVBGGISA-N_2 | -43.8118      | # NON-Bioactive conformation
# | VHJLFUKKYWPYME-JHOBJCJYSA-N | VHJLFUKKYWPYME-JHOBJCJYSA-N_3 | -46.4965      | # NON-Bioactive conformation
# | WWSWXUCLNZWXKO-IADCTJSHSA-O | WWSWXUCLNZWXKO-IADCTJSHSA-O_1 | -38.12        | # Bioactive conformation
# | XBQOGPRXCCTYMI-IEWVHIKDSA-O | XBQOGPRXCCTYMI-IEWVHIKDSA-O_1 | -51.054       | # NON-Bioactive conformation
# | XNOSHWKNNVMTIA-ZTOMLWHTSA-O | XNOSHWKNNVMTIA-ZTOMLWHTSA-O_5 | -48.1715      | # NON-Bioactive conformation
# | YXKHZYNSINIHQT-NGQVCNFZSA-O | YXKHZYNSINIHQT-NGQVCNFZSA-O_5 | -49.2022      | # NON-Bioactive conformation
# | ZDDOCTLDSWCEMH-JHOBJCJYSA-O | ZDDOCTLDSWCEMH-JHOBJCJYSA-O_2 | -40.3028      | # Bioactive conformation
# | ZMTYWFIGVVVFSK-SIBVEZHUSA-N | ZMTYWFIGVVVFSK-SIBVEZHUSA-N_4 | -41.9176      | # NON-Bioactive conformation
# +-----------------------------+-------------------------------+---------------+
```

Use SQL statement to extract the best pose per ligand based on the docked score

```bash
sqlite3 -header -table docking/docking_assays/assay_5/assay_5.db "SELECT LigName, sub_pose, docking_score FROM (SELECT LigName, sub_pose, docking_score, ROW_NUMBER() OVER (PARTITION BY LigName ORDER BY docking_score ASC) as rn FROM Results) ranked WHERE rn = 1;"

################
### Outputs ####
################
# +-----------------------------+-------------------------------+---------------+
# |           LigName           |           sub_pose            | docking_score |
# +-----------------------------+-------------------------------+---------------+
# | BFNJAXKMPJHXHW-DCFHFQCYSA-O | BFNJAXKMPJHXHW-DCFHFQCYSA-O_1 | -5.23         | # Bioactive conformation
# | CLEDCFNOFQFOGX-KMRXNPHXSA-O | CLEDCFNOFQFOGX-KMRXNPHXSA-O_1 | -5.38         | # Bioactive conformation
# | CROWHGLIUOMOOH-SIBVEZHUSA-O | CROWHGLIUOMOOH-SIBVEZHUSA-O_1 | -5.74         | # NON-Bioactive conformation
# | FSMFTHPDNIOSTO-FNZWTVRRSA-O | FSMFTHPDNIOSTO-FNZWTVRRSA-O_1 | -4.31         | # Bioactive conformation
# | FTWLHZZOIPWOOS-IADCTJSHSA-O | FTWLHZZOIPWOOS-IADCTJSHSA-O_1 | -4.72         | # NON-Bioactive conformation
# | GFKPWKCDJBRXQK-DCFHFQCYSA-O | GFKPWKCDJBRXQK-DCFHFQCYSA-O_1 | -5.5          | # Bioactive conformation
# | GJWDUYBKAOXWIS-ZTOMLWHTSA-O | GJWDUYBKAOXWIS-ZTOMLWHTSA-O_1 | -6.29         | # Bioactive conformation
# | LCNNMCFKBKHPHV-DCFHFQCYSA-O | LCNNMCFKBKHPHV-DCFHFQCYSA-O_1 | -5.31         | # NON-Bioactive conformation
# | NRDAGRURPLJETJ-KMRXNPHXSA-O | NRDAGRURPLJETJ-KMRXNPHXSA-O_1 | -5.75         | # Bioactive conformation
# | NXGQNXSIIMHEHD-NTFOOJQESA-O | NXGQNXSIIMHEHD-NTFOOJQESA-O_1 | -5.66         | # Bioactive conformation
# | OQDQGSSUGQXEKZ-JHOBJCJYSA-O | OQDQGSSUGQXEKZ-JHOBJCJYSA-O_1 | -6.56         | # Bioactive conformation
# | PWQIDUKXKZPEMX-AUDOKZQXSA-O | PWQIDUKXKZPEMX-AUDOKZQXSA-O_1 | -7.36         | # Bioactive conformation
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1 | -5.43         | # Bioactive conformation
# | RSSMYFQRFUEXTN-ZTOMLWHTSA-O | RSSMYFQRFUEXTN-ZTOMLWHTSA-O_1 | -5.71         | # NON-Bioactive conformation
# | STUFFQIGUGGQAV-KMRXNPHXSA-O | STUFFQIGUGGQAV-KMRXNPHXSA-O_1 | -5.26         | # Bioactive conformation
# | VCCYIWMCXDRTJQ-MMTVBGGISA-N | VCCYIWMCXDRTJQ-MMTVBGGISA-N_1 | -4.33         | # NON-Bioactive conformation
# | VHJLFUKKYWPYME-JHOBJCJYSA-N | VHJLFUKKYWPYME-JHOBJCJYSA-N_1 | -4.3          | # Bioactive conformation
# | WWSWXUCLNZWXKO-IADCTJSHSA-O | WWSWXUCLNZWXKO-IADCTJSHSA-O_1 | -5.93         | # Bioactive conformation
# | XBQOGPRXCCTYMI-IEWVHIKDSA-O | XBQOGPRXCCTYMI-IEWVHIKDSA-O_1 | -4.88         | # NON-Bioactive conformation
# | XNOSHWKNNVMTIA-ZTOMLWHTSA-O | XNOSHWKNNVMTIA-ZTOMLWHTSA-O_1 | -5.67         | # NON-Bioactive conformation
# | YXKHZYNSINIHQT-NGQVCNFZSA-O | YXKHZYNSINIHQT-NGQVCNFZSA-O_1 | -6.0          | # Bioactive conformation
# | ZDDOCTLDSWCEMH-JHOBJCJYSA-O | ZDDOCTLDSWCEMH-JHOBJCJYSA-O_1 | -6.28         | # NON-Bioactive conformation
# | ZMTYWFIGVVVFSK-SIBVEZHUSA-N | ZMTYWFIGVVVFSK-SIBVEZHUSA-N_1 | -3.79         | # Bioactive conformation
# +-----------------------------+-------------------------------+---------------+
```

Use SQL statement to extract the best pose per ligand based on the cluster size

```bash
sqlite3 -header -table docking/docking_assays/assay_5/assay_5.db "SELECT LigName, sub_pose, cluster_size FROM (SELECT LigName, sub_pose, cluster_size, ROW_NUMBER() OVER (PARTITION BY LigName ORDER BY cluster_size DESC) as rn FROM Results) ranked WHERE rn = 1;"

################
### Outputs ####
################
# +-----------------------------+-------------------------------+--------------+
# |           LigName           |           sub_pose            | cluster_size |
# +-----------------------------+-------------------------------+--------------+
# | BFNJAXKMPJHXHW-DCFHFQCYSA-O | BFNJAXKMPJHXHW-DCFHFQCYSA-O_4 | 11           | # NON-Bioactive conformation
# | CLEDCFNOFQFOGX-KMRXNPHXSA-O | CLEDCFNOFQFOGX-KMRXNPHXSA-O_1 | 18           | # Bioactive conformation
# | CROWHGLIUOMOOH-SIBVEZHUSA-O | CROWHGLIUOMOOH-SIBVEZHUSA-O_2 | 13           | # Bioactive conformation
# | FSMFTHPDNIOSTO-FNZWTVRRSA-O | FSMFTHPDNIOSTO-FNZWTVRRSA-O_3 | 24           | # NON-Bioactive conformation
# | FTWLHZZOIPWOOS-IADCTJSHSA-O | FTWLHZZOIPWOOS-IADCTJSHSA-O_2 | 6            | # NON-Bioactive conformation
# | GFKPWKCDJBRXQK-DCFHFQCYSA-O | GFKPWKCDJBRXQK-DCFHFQCYSA-O_2 | 22           | # NON-Bioactive conformation
# | GJWDUYBKAOXWIS-ZTOMLWHTSA-O | GJWDUYBKAOXWIS-ZTOMLWHTSA-O_3 | 10           | # NON-Bioactive conformation
# | LCNNMCFKBKHPHV-DCFHFQCYSA-O | LCNNMCFKBKHPHV-DCFHFQCYSA-O_3 | 9            | # Bioactive conformation
# | NRDAGRURPLJETJ-KMRXNPHXSA-O | NRDAGRURPLJETJ-KMRXNPHXSA-O_2 | 18           | # NON-Bioactive conformation
# | NXGQNXSIIMHEHD-NTFOOJQESA-O | NXGQNXSIIMHEHD-NTFOOJQESA-O_2 | 21           | # NON-Bioactive conformation
# | OQDQGSSUGQXEKZ-JHOBJCJYSA-O | OQDQGSSUGQXEKZ-JHOBJCJYSA-O_1 | 13           | # Bioactive conformation 
# | PWQIDUKXKZPEMX-AUDOKZQXSA-O | PWQIDUKXKZPEMX-AUDOKZQXSA-O_1 | 26           | # Bioactive conformation
# | RHJLQMVZXQKJKB-FPHSVDBKSA-N | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1 | 28           | # Bioactive conformation
# | RSSMYFQRFUEXTN-ZTOMLWHTSA-O | RSSMYFQRFUEXTN-ZTOMLWHTSA-O_1 | 15           | # NON-Bioactive conformation
# | STUFFQIGUGGQAV-KMRXNPHXSA-O | STUFFQIGUGGQAV-KMRXNPHXSA-O_2 | 17           | # NON-Bioactive conformation
# | VCCYIWMCXDRTJQ-MMTVBGGISA-N | VCCYIWMCXDRTJQ-MMTVBGGISA-N_4 | 19           | # NON-Bioactive conformation
# | VHJLFUKKYWPYME-JHOBJCJYSA-N | VHJLFUKKYWPYME-JHOBJCJYSA-N_4 | 18           | # NON-Bioactive conformation
# | WWSWXUCLNZWXKO-IADCTJSHSA-O | WWSWXUCLNZWXKO-IADCTJSHSA-O_1 | 9            | # Bioactive conformation
# | XBQOGPRXCCTYMI-IEWVHIKDSA-O | XBQOGPRXCCTYMI-IEWVHIKDSA-O_4 | 11           | # NON-Bioactive conformation
# | XNOSHWKNNVMTIA-ZTOMLWHTSA-O | XNOSHWKNNVMTIA-ZTOMLWHTSA-O_2 | 16           | # Bioactive conformation
# | YXKHZYNSINIHQT-NGQVCNFZSA-O | YXKHZYNSINIHQT-NGQVCNFZSA-O_3 | 21           | # NON-Bioactive conformation
# | ZDDOCTLDSWCEMH-JHOBJCJYSA-O | ZDDOCTLDSWCEMH-JHOBJCJYSA-O_2 | 12           | # Bioactive conformation
# | ZMTYWFIGVVVFSK-SIBVEZHUSA-N | ZMTYWFIGVVVFSK-SIBVEZHUSA-N_5 | 8            | # NON-Bioactive conformation
# +-----------------------------+-------------------------------+--------------+
```

Conclusion:


###### hCatL