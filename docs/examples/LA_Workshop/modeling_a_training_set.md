---
title: 4. Evaluation of a training set
---

<div style={{ textAlign: "justify" }}>

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

Let's post-process the docking assay as learned in the previous tutorial, including **pose extraction for visualization** and applying **three alternative ranking criteria** for hit selection: **docking score**, **cluster size**, and **MMGBSA**.

```python title="workshop.py"
# Parse AutoDock outputs, keep top-5 poses per ligand and extract them as PDB
la_workshop_docking_analysis.process_docking_assay(
    assay_id=5,
    max_poses=5,
    extract_poses=1
)

# Computes ΔG estimates for each stored pose
la_workshop_docking_analysis.compute_fingerprints_for_whole_assay(
    assay_id=5,
    mmgbsa=1,
    prolif=0
)  # Provide receptor field indices if prompted

```

Now, for each ligand, we will extract the best pose according to three independent criteria, and compare which one most frequently corresponds to a bioactive conformation (as judged by visual inspection against the reference complex).

Use the following SQL statement to extract the best pose per ligand based on the **docking score**:

```bash
sqlite3 -header -table docking/docking_assays/assay_5/assay_5.db \
"SELECT LigName, sub_pose, docking_score
 FROM (
   SELECT LigName, sub_pose, docking_score,
          ROW_NUMBER() OVER (PARTITION BY LigName ORDER BY docking_score ASC) AS rn
   FROM Results
 ) ranked
 WHERE rn = 1;"
```

<details>
<summary> <b>Output based on docking score</b></summary>

| **LigName**                     | **sub_pose**                    | **docking_score** | **Conformation**           |
|---------------------------------|---------------------------------|------------------:|-----------------------------|
| BFNJAXKMPJHXHW-DCFHFQCYSA-O     | BFNJAXKMPJHXHW-DCFHFQCYSA-O_1   | -5.23             | Bioactive                  |
| CLEDCFNOFQFOGX-KMRXNPHXSA-O     | CLEDCFNOFQFOGX-KMRXNPHXSA-O_1   | -5.38             | Bioactive                  |
| CROWHGLIUOMOOH-SIBVEZHUSA-O     | CROWHGLIUOMOOH-SIBVEZHUSA-O_1   | -5.74             | Non-bioactive              |
| FSMFTHPDNIOSTO-FNZWTVRRSA-O     | FSMFTHPDNIOSTO-FNZWTVRRSA-O_1   | -4.31             | Bioactive                  |
| FTWLHZZOIPWOOS-IADCTJSHSA-O     | FTWLHZZOIPWOOS-IADCTJSHSA-O_1   | -4.72             | Non-bioactive              |
| GFKPWKCDJBRXQK-DCFHFQCYSA-O     | GFKPWKCDJBRXQK-DCFHFQCYSA-O_1   | -5.50             | Bioactive                  |
| GJWDUYBKAOXWIS-ZTOMLWHTSA-O     | GJWDUYBKAOXWIS-ZTOMLWHTSA-O_1   | -6.29             | Bioactive                  |
| LCNNMCFKBKHPHV-DCFHFQCYSA-O     | LCNNMCFKBKHPHV-DCFHFQCYSA-O_1   | -5.31             | Non-bioactive              |
| NRDAGRURPLJETJ-KMRXNPHXSA-O     | NRDAGRURPLJETJ-KMRXNPHXSA-O_1   | -5.75             | Bioactive                  |
| NXGQNXSIIMHEHD-NTFOOJQESA-O     | NXGQNXSIIMHEHD-NTFOOJQESA-O_1   | -5.66             | Bioactive                  |
| OQDQGSSUGQXEKZ-JHOBJCJYSA-O     | OQDQGSSUGQXEKZ-JHOBJCJYSA-O_1   | -6.56             | Bioactive                  |
| PWQIDUKXKZPEMX-AUDOKZQXSA-O     | PWQIDUKXKZPEMX-AUDOKZQXSA-O_1   | -7.36             | Bioactive                  |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N     | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1   | -5.43             | Bioactive                  |
| RSSMYFQRFUEXTN-ZTOMLWHTSA-O     | RSSMYFQRFUEXTN-ZTOMLWHTSA-O_1   | -5.71             | Non-bioactive              |
| STUFFQIGUGGQAV-KMRXNPHXSA-O     | STUFFQIGUGGQAV-KMRXNPHXSA-O_1   | -5.26             | Bioactive                  |
| VCCYIWMCXDRTJQ-MMTVBGGISA-N     | VCCYIWMCXDRTJQ-MMTVBGGISA-N_1   | -4.33             | Non-bioactive              |
| VHJLFUKKYWPYME-JHOBJCJYSA-N     | VHJLFUKKYWPYME-JHOBJCJYSA-N_1   | -4.30             | Bioactive                  |
| WWSWXUCLNZWXKO-IADCTJSHSA-O     | WWSWXUCLNZWXKO-IADCTJSHSA-O_1   | -5.93             | Bioactive                  |
| XBQOGPRXCCTYMI-IEWVHIKDSA-O     | XBQOGPRXCCTYMI-IEWVHIKDSA-O_1   | -4.88             | Non-bioactive              |
| XNOSHWKNNVMTIA-ZTOMLWHTSA-O     | XNOSHWKNNVMTIA-ZTOMLWHTSA-O_1   | -5.67             | Non-bioactive              |
| YXKHZYNSINIHQT-NGQVCNFZSA-O     | YXKHZYNSINIHQT-NGQVCNFZSA-O_1   | -6.00             | Bioactive                  |
| ZDDOCTLDSWCEMH-JHOBJCJYSA-O     | ZDDOCTLDSWCEMH-JHOBJCJYSA-O_1   | -6.28             | Non-bioactive              |
| ZMTYWFIGVVVFSK-SIBVEZHUSA-N     | ZMTYWFIGVVVFSK-SIBVEZHUSA-N_1   | -3.79             | Bioactive                  |

</details>


💡 Based on this benchmark, **docking scores successfully identified** the experimentally known bioactive conformations for roughly **65% of the ligands**.


Now, use the SQL statement to extract the best pose per ligand based on the cluster size (larger = more frequently sampled).

```bash
sqlite3 -header -table docking/docking_assays/assay_5/assay_5.db \
"SELECT LigName, sub_pose, cluster_size
 FROM (
   SELECT LigName, sub_pose, cluster_size,
          ROW_NUMBER() OVER (PARTITION BY LigName ORDER BY cluster_size DESC) AS rn
   FROM Results
 ) ranked
 WHERE rn = 1;"
```

<details>
<summary> <b>Output based on cluster size</b></summary>

| **LigName**                     | **sub_pose**                    | **cluster_size** | **Conformation**           |
|---------------------------------|---------------------------------|-----------------:|-----------------------------|
| BFNJAXKMPJHXHW-DCFHFQCYSA-O     | BFNJAXKMPJHXHW-DCFHFQCYSA-O_4   | 11               | Non-bioactive              |
| CLEDCFNOFQFOGX-KMRXNPHXSA-O     | CLEDCFNOFQFOGX-KMRXNPHXSA-O_1   | 18               | Bioactive                  |
| CROWHGLIUOMOOH-SIBVEZHUSA-O     | CROWHGLIUOMOOH-SIBVEZHUSA-O_2   | 13               | Bioactive                  |
| FSMFTHPDNIOSTO-FNZWTVRRSA-O     | FSMFTHPDNIOSTO-FNZWTVRRSA-O_3   | 24               | Non-bioactive              |
| FTWLHZZOIPWOOS-IADCTJSHSA-O     | FTWLHZZOIPWOOS-IADCTJSHSA-O_2   | 6                | Non-bioactive              |
| GFKPWKCDJBRXQK-DCFHFQCYSA-O     | GFKPWKCDJBRXQK-DCFHFQCYSA-O_2   | 22               | Non-bioactive              |
| GJWDUYBKAOXWIS-ZTOMLWHTSA-O     | GJWDUYBKAOXWIS-ZTOMLWHTSA-O_3   | 10               | Non-bioactive              |
| LCNNMCFKBKHPHV-DCFHFQCYSA-O     | LCNNMCFKBKHPHV-DCFHFQCYSA-O_3   | 9                | Bioactive                  |
| NRDAGRURPLJETJ-KMRXNPHXSA-O     | NRDAGRURPLJETJ-KMRXNPHXSA-O_2   | 18               | Non-bioactive              |
| NXGQNXSIIMHEHD-NTFOOJQESA-O     | NXGQNXSIIMHEHD-NTFOOJQESA-O_2   | 21               | Non-bioactive              |
| OQDQGSSUGQXEKZ-JHOBJCJYSA-O     | OQDQGSSUGQXEKZ-JHOBJCJYSA-O_1   | 13               | Bioactive                  |
| PWQIDUKXKZPEMX-AUDOKZQXSA-O     | PWQIDUKXKZPEMX-AUDOKZQXSA-O_1   | 26               | Bioactive                  |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N     | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1   | 28               | Bioactive                  |
| RSSMYFQRFUEXTN-ZTOMLWHTSA-O     | RSSMYFQRFUEXTN-ZTOMLWHTSA-O_1   | 15               | Non-bioactive              |
| STUFFQIGUGGQAV-KMRXNPHXSA-O     | STUFFQIGUGGQAV-KMRXNPHXSA-O_2   | 17               | Non-bioactive              |
| VCCYIWMCXDRTJQ-MMTVBGGISA-N     | VCCYIWMCXDRTJQ-MMTVBGGISA-N_4   | 19               | Non-bioactive              |
| VHJLFUKKYWPYME-JHOBJCJYSA-N     | VHJLFUKKYWPYME-JHOBJCJYSA-N_4   | 18               | Non-bioactive              |
| WWSWXUCLNZWXKO-IADCTJSHSA-O     | WWSWXUCLNZWXKO-IADCTJSHSA-O_1   | 9                | Bioactive                  |
| XBQOGPRXCCTYMI-IEWVHIKDSA-O     | XBQOGPRXCCTYMI-IEWVHIKDSA-O_4   | 11               | Non-bioactive              |
| XNOSHWKNNVMTIA-ZTOMLWHTSA-O     | XNOSHWKNNVMTIA-ZTOMLWHTSA-O_2   | 16               | Bioactive                  |
| YXKHZYNSINIHQT-NGQVCNFZSA-O     | YXKHZYNSINIHQT-NGQVCNFZSA-O_3   | 21               | Non-bioactive              |
| ZDDOCTLDSWCEMH-JHOBJCJYSA-O     | ZDDOCTLDSWCEMH-JHOBJCJYSA-O_2   | 12               | Bioactive                  |
| ZMTYWFIGVVVFSK-SIBVEZHUSA-N     | ZMTYWFIGVVVFSK-SIBVEZHUSA-N_5   | 8                | Non-bioactive              |
</details>

💡 When ranked by **cluster size**, **only about 39%** of the ligands were correctly associated with their experimentally known **bioactive conformations**.

Next, explore a third ranking criterion based on MMGBSA binding free energies.

```bash
sqlite3 -header -table docking/docking_assays/assay_5/assay_5.db \
"SELECT LigName, sub_pose, delta_g_total
 FROM (
   SELECT LigName, sub_pose, delta_g_total,
          ROW_NUMBER() OVER (PARTITION BY LigName ORDER BY delta_g_total ASC) AS rn
   FROM mmgbsa_fingerprints
 ) ranked
 WHERE rn = 1;"
```

<details>
<summary> <b>Output based on MMGBSA total energy</b></summary>
| **LigName**                     | **sub_pose**                    | **delta_g_total** | **Conformation**           |
|---------------------------------|---------------------------------|------------------:|-----------------------------|
| BFNJAXKMPJHXHW-DCFHFQCYSA-O     | BFNJAXKMPJHXHW-DCFHFQCYSA-O_4   | -51.3374          | Non-bioactive              |
| CLEDCFNOFQFOGX-KMRXNPHXSA-O     | CLEDCFNOFQFOGX-KMRXNPHXSA-O_2   | -34.4763          | Bioactive                  |
| CROWHGLIUOMOOH-SIBVEZHUSA-O     | CROWHGLIUOMOOH-SIBVEZHUSA-O_2   | -34.6089          | Bioactive                  |
| FSMFTHPDNIOSTO-FNZWTVRRSA-O     | FSMFTHPDNIOSTO-FNZWTVRRSA-O_2   | -45.7016          | Non-bioactive              |
| FTWLHZZOIPWOOS-IADCTJSHSA-O     | FTWLHZZOIPWOOS-IADCTJSHSA-O_5   | -53.6438          | Non-bioactive              |
| GFKPWKCDJBRXQK-DCFHFQCYSA-O     | GFKPWKCDJBRXQK-DCFHFQCYSA-O_4   | -51.5091          | Non-bioactive              |
| GJWDUYBKAOXWIS-ZTOMLWHTSA-O     | GJWDUYBKAOXWIS-ZTOMLWHTSA-O_1   | -50.1655          | Bioactive                  |
| LCNNMCFKBKHPHV-DCFHFQCYSA-O     | LCNNMCFKBKHPHV-DCFHFQCYSA-O_4   | -36.9908          | Non-bioactive              |
| NRDAGRURPLJETJ-KMRXNPHXSA-O     | NRDAGRURPLJETJ-KMRXNPHXSA-O_4   | -52.3813          | Non-bioactive              |
| NXGQNXSIIMHEHD-NTFOOJQESA-O     | NXGQNXSIIMHEHD-NTFOOJQESA-O_3   | -47.5815          | Non-bioactive              |
| OQDQGSSUGQXEKZ-JHOBJCJYSA-O     | OQDQGSSUGQXEKZ-JHOBJCJYSA-O_4   | -52.5189          | Non-bioactive              |
| PWQIDUKXKZPEMX-AUDOKZQXSA-O     | PWQIDUKXKZPEMX-AUDOKZQXSA-O_4   | -38.6934          | Non-bioactive              |
| RHJLQMVZXQKJKB-FPHSVDBKSA-N     | RHJLQMVZXQKJKB-FPHSVDBKSA-N_1   | -50.2251          | Bioactive                  |
| RSSMYFQRFUEXTN-ZTOMLWHTSA-O     | RSSMYFQRFUEXTN-ZTOMLWHTSA-O_4   | -41.0921          | Non-bioactive              |
| STUFFQIGUGGQAV-KMRXNPHXSA-O     | STUFFQIGUGGQAV-KMRXNPHXSA-O_5   | -53.2883          | Non-bioactive              |
| VCCYIWMCXDRTJQ-MMTVBGGISA-N     | VCCYIWMCXDRTJQ-MMTVBGGISA-N_2   | -43.8118          | Non-bioactive              |
| VHJLFUKKYWPYME-JHOBJCJYSA-N     | VHJLFUKKYWPYME-JHOBJCJYSA-N_3   | -46.4965          | Non-bioactive              |
| WWSWXUCLNZWXKO-IADCTJSHSA-O     | WWSWXUCLNZWXKO-IADCTJSHSA-O_1   | -38.1200          | Bioactive                  |
| XBQOGPRXCCTYMI-IEWVHIKDSA-O     | XBQOGPRXCCTYMI-IEWVHIKDSA-O_1   | -51.0540          | Non-bioactive              |
| XNOSHWKNNVMTIA-ZTOMLWHTSA-O     | XNOSHWKNNVMTIA-ZTOMLWHTSA-O_5   | -48.1715          | Non-bioactive              |
| YXKHZYNSINIHQT-NGQVCNFZSA-O     | YXKHZYNSINIHQT-NGQVCNFZSA-O_5   | -49.2022          | Non-bioactive              |
| ZDDOCTLDSWCEMH-JHOBJCJYSA-O     | ZDDOCTLDSWCEMH-JHOBJCJYSA-O_2   | -40.3028          | Bioactive                  |
| ZMTYWFIGVVVFSK-SIBVEZHUSA-N     | ZMTYWFIGVVVFSK-SIBVEZHUSA-N_4   | -41.9176          | Non-bioactive              |
</details>

💡 When ranked by **MMGBSA total energies**, **only 6 out of 23 ligands (~26%) matched their experimentally known bioactive conformations**; the lowest success rate among the three criteria tested.


:::note[📊 Comparative interpretation]

Across both validation stages (the single-ligand K777 benchmark and the broader training set analysis) the docking score consistently emerged as the most reliable ranking criterion for identifying bioactive conformations. While cluster size offered limited correlation with biological relevance (≈39% accuracy) and MMGBSA refinements further reduced predictive performance (≈26%), the docking energy values successfully reproduced known binding modes in approximately two-thirds of the cases.

Together, these results suggest that, within the current setup and receptor parametrization, docking scores provide the most robust first-level filter for pose selection and hit identification, whereas MMGBSA or cluster-based refinements may serve as complementary but not decisive tools for subsequent energetic or stability assessments.
:::

</div>
