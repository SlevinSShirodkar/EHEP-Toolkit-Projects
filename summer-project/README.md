# C++ ROOT Analysis

This repository contains a custom C++ analysis framework built on the ROOT data analysis ecosystem. The project's primary objective was to independently reverse-engineer and deduce the underlying Feynman diagrams of four unlabelled Monte Carlo (MC) simulation samples by applying fundamental kinematics, particle object selection, and resonance plotting.

## Technical Skills Demonstrated
* **1)C++ & ROOT Framework:** Engineered analysis scripts to iterate over thousands of simulated events, filling custom histograms and applying complex physics logic[cite: 8, 12].
* **2)Particle Object Selection:** Programmed identification criteria for Muons, Electrons, Photons, Jets, and b-tagged Jets, applying strict transverse momentum ($p_T$), pseudorapidity ($\eta$), isolation ($\Delta R$), and ID cuts.
* **3)Kinematics & Dynamics:** Calculated and plotted multi-body invariant masses, missing transverse energy (MET), and transverse mass ($m_T$) to reconstruct unknown parent particles.
* **4)Truth-Matching:** Developed logic to match reconstructed leptons to their true generator-level counterparts (GenPart) utilizing $\Delta R$ minimization and tracing parent-child decay lineage (`GenMother`).

## Process Deduction & Analytical Logic
Successfully identified the physical processes governing four mystery MC datasets (Labelled A, B, C and D). By iteratively analyzing kinematic distributions (such as $M_{\mu\mu}$ and MET), applying targeted phase-space cuts (e.g., Z-boson mass windows and b-tagging requirements), and cross-referencing multijet and dilepton topologies, I deduced the underlying Standard Model processes purely from the simulated event data.

Few of the plots I generated using the mystery MC samples' information are uploaded to the repository. 
