[![powered by pyGIMLi](https://img.shields.io/badge/powered%20by-pyGIMLi-informational?style=flat&logo=python&logoColor=white)](https://www.pygimli.org)
# Minimum entropy constrained cooperative inversion of electrical resistivity, seismic and magnetic data <br>- Complementary Repository

By [Anton Ziegon](https://www.gim.rwth-aachen.de/team/alumni/), [Dr. Marc Boxberg](https://www.gim.rwth-aachen.de/team/marc-boxberg/) and [Prof. Dr. Florian Wagner](https://www.gim.rwth-aachen.de/team/florian-wagner/) of the [Teaching and Research Unit Geophysical Imaging and Monitoring (GIM)](https://www.gim.rwth-aachen.de/about/) of RWTH Aachen University.

![GIM Logo](https://www.gim.rwth-aachen.de/images/logos/gim_logo.svg)

---
## Abstract
Geophysical methods are widely used to gather information about the subsurface as they are non-intrusive and comparably cheap. However, the solution to the geophysical inverse problem is inherently non-unique, which introduces considerable uncertainties. As a partial remedy to this problem, independently acquired geophysical data sets can be jointly inverted to reduce ambiguities in the resulting multi-physical subsurface images. A novel cooperative inversion approach with joint minimum entropy constraints is used to create more consistent multi-physical images with sharper boundaries. Here, this approach is implemented in an open-source software and its applicability on electrical resistivity tomography (ERT), seismic refraction tomography (SRT) and magnetic data is investigated. A synthetic 2D ERT and SRT data study is used to demonstrate the approach and to investigate the influence of the governing parameters. The findings showcase the advantage of the joint minimum entropy (JME) stabilizer over separate, conventional smoothness-constrained inversions. The method is then used to analyze field data from Rockeskyller Kopf, Germany. 3D ERT and magnetic data is combined and results confirm the expected volcanic diatreme structure with improved details. The multi-physical images of both methods are consistent in some regions as similar boundaries are produced in the resulting models. Because of its sensitivity to hydrologic conditions in the subsurface, observations suggest that the ERT method senses different structures than the magnetic method. These structures in the ERT result do not seem to be enforced on the magnetic susceptibility distribution, showcasing the flexibility of the approach. Both investigations outline the importance of a suitable parameter and reference model selection for the performance of the approach and suggest careful parameter tests prior to the joint inversion. With proper settings, the JME inversion is a promising tool for geophysical imaging, however, this work also identifies some objectives for future studies and additional research to explore and optimize the method.

## Structure of Repostory
- [Code](./Code): Contains...
    - *JointEntropyInversion* class
    - Notebooks to produce figures and results (first number of filename indicates the chapter in the thesis that is affiliated with it)
    - Scripts with helper functions for plotting

-  [Data](./Data):
    - [Synthetic](./Data/Synthetic): Contains synthetic data used for the study
    - [Rockeskyll](./Data/Rockeskyll): Folder associated with the field data study at Rockeskyller Kopf, Germany. Contains...
        - Folder with histogram filtered ERT data of all 4 lines ([ERT_raw]('./Data/Rockeskyll/ERT_raw'))
        - Folder with acquisition geometry in 3D ([Geo](./Data/Rockeskyll/Geo))
        - Magnetic data (*Magnetic_data_corrected.csv*)
        - Output of Data preparation notebook (4.0)
        - 2D/3D Meshes (*.bms*)
        - Interpolation indices for extraction of 2D cross-sections out of 3D volume (*.npy*, see Notebook 4.4)
        - Folder with extra scripts that could help visualize additional results ([extra_scripts](./Data/Rockeskyll/extra_scripts))
        - Inversion results of Notebooks 4.1-4.3 ([Res_Conventional](./Data/Rockeskyll/Res_Conventional), [Res_ME1](./Data/Rockeskyll/Res_ME1),  [Res_ME2](./Data/Rockeskyll/Res_ME2),  [Res_JME1](./Data/Rockeskyll/Res_JME1),  [Res_JME2](./Data/Rockeskyll/Res_JME2))
        - Mesh containing JME inversion results in *.vtk*-format (for visualization in ParaView)
  - [Figures](./Figures):
      - [Ch-2](./Figures/Ch-2): Figures of chapter 2 (Theory & Methods)
      - [Ch-3](./Figures/Ch-3): Figures of chapter 3 (Synth. data study)
      - [Ch-4](./Figures/Ch-4): Figures of chapter 4 (Field data study: Rockeskyller Kopf) including ParaView screenshots
      - [Ch-5](./Figures/Ch-5): Figures of chapter 5 (Discussion)
      - Flowchart corresponding to the *Joint EntropyInversion* class
