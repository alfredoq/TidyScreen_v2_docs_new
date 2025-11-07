---
title: 1. Design of Targeted Covalent Inhibitors (TCIs) targeting Cruzipain
---
<div style={{ textAlign: "justify" }}>

## 1.1- Chagas Disease: A Neglected Global Challenge

Chagas disease (CD), also known as American trypanosomiasis, is a potentially life-threatening systemic infection caused by the protozoan parasite *Trypanosoma cruzi* (*T. cruzi*). It affects approximately 6.5 million people worldwide, with an annual incidence of about 30,000 new cases and up to 13,000 deaths. Despite its impact, CD remains one of the most neglected tropical diseases and the third most common parasitic infection globally, underscoring the urgent need for new public health and therapeutic strategies in endemic regions.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/LancetChagasMap.jpg" alt="Description of image" width="700"/>
  <figcaption>Estimated worldwide prevalence of Chagas disease. </figcaption>
  </p>
</figure>
---

Although significant efforts have been devoted to drug discovery for CD, **only two nonspecific agents** developed nearly five decades ago - **nifurtimox** (1940) and **benznidazole** (1974) - remain the sole treatment options. While both drugs are effective against the acute stage, their severe adverse effects and limited efficacy in chronic infection restrict their use, leading to poor adherence and the emergence of resistance.

## 1.2- Cruzipain (CZP): A Key Druggable Target

By the late 1970s, a family of cysteine proteases was identified as essential to the *T. cruzi* life cycle. Among them, **Cruzipain** (CZP), characterized by Cazzulo and colleagues, emerged as the predominant enzyme and a key druggable target. Beyond its role in parasite nutrition, CZP participates in metacyclogenesis and contributes to host cell invasion and virulence.

Structurally, CZP’s catalytic domain comprises 215 amino acids forming α-helices and antiparallel β-sheets arranged into two subdomains that define the active site. This site contains a **catalytic triad (Cys25, His162, Asn182)** responsible for hydrolysis and **four binding subsites (S1′, S1, S2, and S3)** that dictate substrate and inhibitor recognition. Extensive structure-activity relationship (SAR) studies have elucidated the distinct roles of these subsites in ligand binding and specificity.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/CZP_Subsites.png" alt="Description of image" width="500"/>
  <figcaption>3D structure of CZP (PDB ID: [3IUT](https://www.rcsb.org/structure/3IUT)) active site, depicting S1' (red), S1 (teal), S2 (blue), S3 (magenta) subsites, and the region corresponding to the oxyanion hole and catalytic triad (OAH+CT, olive). </figcaption>
  </p>
</figure>
---

A wide variety of compounds have been reported as CZP inhibitors, ranging from natural peptides to synthetic small molecules identified through drug repurposing and de novo design. These inhibitors act through either **non-covalent or covalent mechanisms**. 

- **Non-covalent inhibitors** interact reversibly within the active site, governed by an equilibrium constant (Ki), forming an energetically favourable enzyme-inhibitor encounter complex. 
- **Covalent inhibitors** incorporate a reactive warhead (WH) designed to form a chemical bond with the nucleophilic catalytic residue Cys25. The initial non-covalent recognition ensures proper WH orientation before covalent bond formation, characterized by a rate constant (k<sub>inact</sub>). When this bond formation is irreversible, inhibition is permanent; if reversible, both rate constants are comparable, leading to a dynamic equilibrium. 

Covalent inhibitors can achieve greater potency and longer duration of action but often raise safety concerns due to potential off-target reactivity. Such limitations can be mitigated by improving target selectivity and using WHs of controlled reactivity. This gave rise to **Targeted Covalent Inhibitors (TCIs)**, compounds bearing a reactive moiety specifically designed to form a covalent bond with a defined residue of the intended target, thereby ensuring potency with minimal off-target effects. Although their design requires sophisticated computational and structural approaches, several TCIs have successfully reached clinical use, validating this concept.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/TCI.png" alt="Description of image" width="700"/>
  <figcaption>Schematic representation of TCIs' inhibitory mechanism. </figcaption>
  </p>
</figure>
---


## 1.3- K777: Benchmark CZP Inhibitor

Among CZP inhibitors, **K777** stands out as the benchmark compound. It is a peptidomimetic CZP inhibitor (in vitro IC<sub>50</sub> = 2-10 nM) featuring a vinylsulfone WH whose mechanism has been confirmed both experimentally and computationally. K777 exhibits potent in vitro and in vivo activity across different *T. cruzi* strains and CZP subtypes, including those resistant to standard therapies. The crystal structure of CZP covalently bound to K777 (PDB ID: [2OZ2](https://www.rcsb.org/structure/2OZ2)) revealed critical interactions: its peptide backbone forms hydrogen bonds with Gly65 (S3) and Gly66 (S2); the sulfonyl group aligns with the oxyanion hole; and the methylpiperazine at R<sup>3</sup> establishes a stabilizing contact with Ser61 (S3). Modifications or removal of the homoPhe-like R<sup>1</sup> drastically reduce potency, and attempts to replace the phenyl group at R<sup>2</sup> failed to improve activity.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/K777IPandBioac.png" alt="Description of image" width="400"/>
  <figcaption>Intermolecular interactions within the CZP active site and biological data of K777. </figcaption>
  </p>
</figure>
---

Preclinical studies confirmed K777’s oral bioavailability and efficacy in rodents, dogs, and non-human primates, eliminating parasitemia in both acute and chronic infections. However, they were halted due to elevated hepatotoxicity markers - likely resulting from irreversible inhibition of CYP3A4 - and other adverse effects, even at low doses. Consequently, ongoing research aims to design analogues or alternative delivery strategies to preserve efficacy while minimizing toxicity.

## 1.4- CZP vs. Human Cathepsins: The Selectivity Challenge

Importantly, CZP shares significant structural similarity with the **human cathepsin family** of cysteine proteases, which participate in diverse physiological processes such as immune regulation, tissue remodeling, and hormone processing. Among them, cathepsin L (hCatL) exhibits the highest homology (>50%) with CZP. This resemblance poses a major selectivity challenge, as inhibition of hCatL can lead to undesired side effects. Notably, **K777 potently inhibits hCatL** (in vitro IC<sub>50</sub> = 0.2 nM), leading to its repurposing as **SLV213** for potential treatment of SARS-CoV-2, where hCatL mediates viral entry. Thus, minimizing hCatL inhibition remains a critical objective in the rational design of CZP-directed inhibitors.

---
<figure>
  <p align="center">
  <img src="/TidyScreen_v2_docs_new/img/CZPhCatL.png" alt="Description of image" width="400"/>
  <figcaption>Superimposed 3D-structures of CZP (cyan) and hCatL (red). </figcaption>
  </p>
</figure>
---


## 1.5- Perspectives and Outlook

The success and limitations of K777 collectively provide a compelling **proof of concept** for CZP as a therapeutic target.  
Future development of **selective, non-peptidic TCIs** should focus on:

- Optimizing **non-covalent recognition** across the four CZP subsites  
- Achieving precise control over **WH reactivity** toward Cys25  

This balanced design is essential to maximize **potency** and **selectivity**, minimize **off-target effects**, and overcome the **poor bioavailability** and **metabolic instability** typically associated with peptide-based inhibitors.  

</div>
