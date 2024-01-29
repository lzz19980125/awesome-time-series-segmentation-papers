# Awesome Time Series Segmentation Papers

[![Awesome](https://awesome.re/badge.svg)](https://awesome.re) 

This repository contains a reading list of papers on **Time Series Segmentation**. **This repository is still being continuously improved.**

As a crucial time series preprocessing technique, **semantic segmentation** divides poorly understood time series into several discrete and homogeneous segments. This approach aims to uncover latent temporal evolution patterns, detect unexpected regularities and regimes, thereby rendering the analysis of massive time series data more manageable. 

Time series segmentation often intertwines with research in many domains. Firstly, the relationship between **time series segmentation, time series change point detection,** and some aspects of **time series anomaly/outlier detection** is somewhat ambiguous. Therefore, this repository includes a selection of papers from these areas. Secondly, time series segmentation can be regarded as a process of information compression in time series, hence papers in this field often incorporate concepts from **information theory** (e.g., using minimum description length to guide the design of unsupervised time series segmentation models). Additionally, the task of decomposing human actions into a series of plausible motion primitives can be addressed through methods for segmenting sensor time series. Consequently, papers related to motion capture from the fields of **computer vision** and **ubiquitous computing** are also included in this collection.

Generally, the subjects of unsupervised semantic segmentation can be categorized into:

* ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) univariate time series: ![](https://latex.codecogs.com/svg.image?\inline&space;1\times&space;T), where  ![](https://latex.codecogs.com/svg.image?\inline&space;T) is the length of the time series.
* ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) multivariate time series: ![](https://latex.codecogs.com/svg.image?\inline&space;M\times&space;T), where ![](https://latex.codecogs.com/svg.image?\inline&space;M) is the number of variables (channels).
* ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) tensor: ![](https://latex.codecogs.com/svg.image?\inline&space;N_{1}&space;\times&space;\cdots&space;\times&space;N_{k}&space;\times&space;M\times&space;T), where ![](https://latex.codecogs.com/svg.image?\inline&space;N_{1}&space;\times&space;\cdots&space;\times&space;N_{k}) denotes the dimensions other than time and variables.

In the field of time series research, unlike time series forecasting, anomaly detection, and classification/clustering, the number of papers on time series segmentation has been somewhat lukewarm in recent years (this observation may carry a degree of subjectivity from the author). Additionally, deep learning methods do not seem to dominate this area as they do in others. Some classic but solid algorithms remain highly competitive even today, with quite a few originating from the same research group. Therefore, in the following paper list, I will introduce them indexed by well-known researchers and research groups in this field.

## Some Additional Information

🚩 2024/1/27: **I have marked some recommended papers / datasets / implementations with 🌟 (Just my personal preference 😉).**

## Survey & Evaluation

NOTE: the ranking has no particular order.

|     TYPE      |               Venue               |             Paper Title and Paper Interpretation             |                             Code                             |
| :-----------: | :-------------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|    Dataset    |     *DARLI-AP@EDBT/ICDT '23*      | Time Series Segmentation Applied to a New Data Set for Mobile Sensing of Human Activities 🌟 | [MOSAD](https://github.com/ermshaua/mobile-sensing-human-activity-data-set) |
|    Dataset    |     *ECML-PKDD Workshop '23*      |     Human Activity Segmentation Challenge@ECML/PKDD’23 🌟     | [Challenge Link](https://ecml-aaltd.github.io/aaltd2023/challenge.html) |
|    Survey     | *IEEE J. Sel.  Areas Commun. '21* | Sequential (Quickest) Change Detection Classical Results and New Directions |                             None                             |
| Visualization |          *IEEE TVCG '21*          | MultiSegVA Using Visual Analytics to Segment Biologging Time Series on Multiple Scales |                             None                             |
|    Survey     | *IEEE J. Sel. Areas Commun. '21*  | Sequential (Quickest) Change Detection Classical Results and New Directions |                             None                             |
|    Survey     |       *Signal Process. '20*       | Selective review of offline change point detection methods 🌟 | [Ruptures](https://centre-borelli.github.io/ruptures-docs/)  |
|  Evaluation   |            *Arxiv '20*            |     An Evaluation of Change Point Detection Algorithms 🌟     | [TCPDBench](https://github.com/alan-turing-institute/TCPDBench) |
|    Survey     |      *Knowl. Inf. Syst. '17*      | A survey of methods for time series change point detection 🌟 |                             None                             |
|  Evaluation   |         *Inf. Syst. '17*          | An evaluation of combinations of lossy compression and change-detection approaches for time-series data |                             None                             |
|    Survey     | *IEEE Trans Hum. Mach. Syst. '16* | Movement Primitive Segmentation for Human Motion Modeling A Framework for Analysis 🌟 |                             None                             |
|    Survey     |            *EAAI '11*             |             A review on time series data mining              |                             None                             |
|    Survey     |            *CSUR '11*             |                   Time-series data mining                    |                             None                             |
|    Dataset    |             *GI '04*              |    Segmenting Motion Capture Data into Distinct Behaviors    | [Website](http://graphics.cs.cmu.edu/projects/segmentation/) 🌟 |
## [David Hallac (Stanford)](https://scholar.google.com/citations?hl=zh-CN&user=o23Mn0cAAAAJ&view_op=list_works&sortby=pubdate)
|                             TYPE                             |             Venue              |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :----------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *KDD Workshop MiLeTS '20*    |    Driver2vec Driver Identification from Automotive Data     |    [Driver2vec](https://github.com/JingboYang/driver2vec)    |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *Adv. Data Anal. Classif. '19* |  Greedy Gaussian segmentation of multivariate time series 🌟  |             [GGS](https://github.com/cvxgrp/GGS)             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          *Arxiv '18*           | MASA: Motif-Aware State Assignment in Noisy Time Series Data | [MASA](https://github.com/snap-stanford/masa?utm_source=catalyzex.com) |
|                        *Ph.D. Thesis*                        |         *ProQuest '18*         | Inferring Structure from Multivariate Time Series Sensor Data |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |           *KDD '17*            | Toeplitz Inverse Covariance-Based Clustering of Multivariate Time Series Data 🌟 |         [TICC](https://github.com/davidhallac/TICC)          |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |           *KDD '17*            |   Network Inference via the Time-Varying Graphical Lasso 🌟   |         [TVGL](https://github.com/davidhallac/TVGL)          |
## [Shaghayegh Gharghabi](https://scholar.google.com/citations?hl=en&user=EITBC1YAAAAJ&view_op=list_works&sortby=pubdate) (from [Eamonn Keogh](https://www.cs.ucr.edu/~eamonn/)'s Lab, UC Riverside)
|                             TYPE                             |   Venue    |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :--------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *DMKD '19* | Domain agnostic online semantic segmentation for multi-dimensional time series 🌟 | [Floss & datasets](https://www.cs.ucr.edu/%7eeamonn/FLOSS/)  |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) | *ICDM '17* | Matrix Profile VIII Domain Agnostic Online Semantic Segmentation at Superhuman Performance Levels 🌟 | [Floss](https://sites.google.com/site/onlinesemanticsegmentation/) |

## [Yasuko Matsubara](https://www.dm.sanken.osaka-u.ac.jp/~yasuko/) & [Yasushi Sakurai](https://www.dm.sanken.osaka-u.ac.jp/~yasushi/) (from [Sakurai & Matsubara Lab](https://www.dm.sanken.osaka-u.ac.jp/))
|                             TYPE                             |     Venue     |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :-----------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |    WWW '23    | Fast and Multi-aspect Mining of Complex Time-stamped Event Streams 🌟 |      [CubeScope](https://github.com/kotaNakm/CubeScope)      |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |   *KDD '22*   | Fast Mining and Forecasting of Co-evolving Epidemiological Data Streams 🌟 |                             None                             |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |  *CIKM '22*   |      Modeling Dynamic Interactions over Tensor Streams       |          [Dismo](https://github.com/kokikwbt/dismo)          |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *CIKM '22*   | Mining Reaction and Diffusion Dynamics in Social Activities 🌟 |                             None                             |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) | *NeurIPS '21* |         SSMF Shifting Seasonal Matrix Factorization          |           [ssmf](https://github.com/kokikwbt/ssmf)           |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |   *KDD '20*   |  Non-Linear Mining of Social Activities in Tensor Streams 🌟  |                             None                             |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |  *ICDM '19*   |      Multi-aspect mining of complex sensor sequences 🌟       |   [CubeMarker](https://github.com/TakatoHonda/CubeMarker)    |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *KDD '19*   | Dynamic Modeling and Forecasting of Time-evolving Data Streams |   [OrbitMap](https://github.com/yasuko-matsubara/orbitmap)   |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *CIKM '19*   |     Automatic Sequential Pattern Mining in Data Streams      |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *KDD '16*   | Regime Shifts in Streams: Real-time Forecasting of Co-evolving Time Sequences | [RegimeCast](https://www.dm.sanken.osaka-u.ac.jp/~yasuko/SRC/regimecast.zip) |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |   *WWW '16*   |       Non-linear mining of competing local activities        | [CompCube](https://www.dm.sanken.osaka-u.ac.jp/~yasuko/SRC/compcube.zip) |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |    WWW '15    | The web as a jungle: Non-linear dynamical systems for co-evolving online activities 🌟 | [Ecoweb & dataset](https://www.dm.sanken.osaka-u.ac.jp/~yasuko/SRC/ecoweb.zip) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *SIGMOD '14*  |  AutoPlait Automatic Mining of Co-evolving Time Sequences 🌟  | [AutoPlait](https://www.dm.sanken.osaka-u.ac.jp/~yasuko/SRC/autoplait.zip) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *ICDM '14*   |    Fast and Exact Monitoring of Co-evolving Data Streams     |                             None                             |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |   *KDD '14*   |  FUNNEL Automatic Mining of Spatially Coevolving Epidemics   | [Funnel](https://www.dm.sanken.osaka-u.ac.jp/~yasuko/SRC/funnel.zip) |

## [Bryan Hooi](https://bhooi.github.io/) (NUS)
|                             TYPE                             |     Venue      |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *TKDE '22*   | Time Series Anomaly Detection with Adversarial Reconstruction Networks 🌟 | [BeatGAN](https://github.com/BGT-M/spartan2-tutorials/blob/master/BeatGAN.ipynb) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *IJCAI '19*   | BeatGAN Anomalous Rhythm Detection using Adversarially Generated Time Series 🌟 |        [BeatGAN](https://github.com/hi-bingo/BeatGAN)        |
|                        *Ph.D. Thesis*                        | *ProQuest '19* | Anomaly Detection in Graphs and Time Series Algorithms and Applications |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *SDM '19*    | Branch and Border Partition Based Change Detection in Multivariate Time Series 🌟 |         [Bnb](https://bhooi.github.io/code/bnb.zip)          |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |   *SDM '19*    | SMF Drift-Aware Matrix Factorization with Seasonal Patterns  |        [smf & dataset](https://bhooi.github.io/smf/)         |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) |   *WWW '17*    | AutoCyclone Automatic Mining of Cyclic Online Activities with Robust Tensor Factorization | [AutoCyclone](http://www.cs.cmu.edu/~tsubasat/code/AutoCyclone.zip) |

## [Liangzhe Chen](https://sites.google.com/view/liangzhechen) & [Anika Tabassum](https://scholar.google.com/citations?hl=zh-CN&user=4Aw4_1sAAAAJ&view_op=list_works&sortby=pubdate) (Virginia Tech, supervised by [B. Aditya Prakash](https://scholar.google.com/citations?hl=zh-CN&user=C-NftTgAAAAJ&view_op=list_works&sortby=pubdate))
|                             TYPE                             |     Venue      |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
|                        *Ph.D. Thesis*                        | *ProQuest '21* | Explainable and Network-Based Approaches for Decision-making in Emergency Management |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *CIKM '21*   |    Actionable Insights in Urban Multivariate Time-series     |         [RaTSS](https://github.com/AdityaLab/RaTSS)          |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *TIST '20*   | Cut-n-Reveal: Time-Series Segmentations with Explanations 🌟  |  [Cut-n-Reveal](https://github.com/AdityaLab/Cut-n-Reveal)   |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *AAAI '18*   |           Automatic Segmentation of Data Sequences           | [DASSA](https://github.com/lzz19980125/awesome-time-series-segmentation-papers/tree/main/DASSA-master) |
|                        *Ph.D. Thesis*                        | *ProQuest '18* |    Segmenting, Summarizing and Predicting Data Sequences     |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *vt.edu '18*  |    Segmentations with Explanations for Outage Analysis 🌟     |                             None                             |

## [Shohreh Deldari](https://scholar.google.com/citations?hl=zh-CN&user=250tnREAAAAJ&view_op=list_works&sortby=pubdate) (from [Cruise research group](https://cruiseresearchgroup.github.io/), RMIT ) &  [Flora D. Salim](https://scholar.google.com/citations?hl=zh-CN&user=Yz35RSYAAAAJ) (UNSW)

|                             TYPE                             |            Venue             |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :--------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |         *Arxiv '23*          | Detecting Change Intervals with Isolation Distributional Kernel 🌟 | [ICD](https://github.com/IsolationKernel/Codes/tree/main/IDK/iCID) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |         *IMWUT '22*          | COCOA Cross Modality Contrastive Learning for Sensor Data 🌟  |    [COCOA](https://github.com/cruiseresearchgroup/COCOA)     |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          *WWW '21*           | Time Series Change Point Detection with Self-Supervised Contrastive Predictive Coding 🌟 |    [TSCP2](https://github.com/cruiseresearchgroup/TSCP2)     |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |         *IMWUT '20*          | ESPRESSO Entropy and ShaPe awaRe timE-Series SegmentatiOn for Processing Heterogeneous Sensor Data | [ESPRESSO](https://github.com/cruiseresearchgroup/ESPRESSO)  |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *Knowl. Inf. Syst. '20*    | Unsupervised online change point detection in high-dimensional time series |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |     *WSDM Workshop '19*      | Inferring Work Routines and Behavior Deviations with Life-logging Sensor Data |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *Pervasive Mob. Comput. '17* | Information gain-based metric for recognizing transitions in human activities 🌟 |  [IGTs](https://github.com/cruiseresearchgroup/IGTS-python)  |

## [Peng Wang](https://scholar.google.com/citations?hl=en&user=fxcAZkoAAAAJ&view_op=list_works&sortby=pubdate) (fudan University)

|                             TYPE                             |    Venue     |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :----------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |  *ICDE '21*  | GRAB: Finding Time Series Natural Structures via A Novel Graph-based Scheme | [GRAB](https://github.com/lzz19980125/awesome-time-series-segmentation-papers/tree/main/GRAB-master) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *SIGMOD '11* |              Finding Semantics in Time Series 🌟              |                             None                             |

## [Arik Ermshaus](https://sites.google.com/view/arikermshaus) (Humboldt-Universität zu Berlin)

|                             TYPE                             |          Venue           |           Paper Title and Paper Interpretation            |                             Code                             |
| :----------------------------------------------------------: | :----------------------: | :-------------------------------------------------------: | :----------------------------------------------------------: |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |       *Arxiv '23*        | Raising the ClaSS of Streaming Time Series Segmentation 🌟 | [Clasp](https://github.com/ermshaua/classification-score-stream) |
|                           Dataset                            | *ECML-PKDD Workshop '23* |   Human Activity Segmentation Challenge@ECML/PKDD’23 🌟    | [Challenge Link](https://ecml-aaltd.github.io/aaltd2023/challenge.html) |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |        *DMKD '23*        |     ClaSP: parameter-free time series segmentation 🌟      | [Clasp](https://sites.google.com/view/ts-parameter-free-clasp/) |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |        *CIKM '21*        |            ClaSP - Time Series Segmentation 🌟             |       [Clasp](https://sites.google.com/view/ts-clasp/)       |

## [Lei Li](https://lileicc.github.io/) (CMU)

|                             TYPE                             |     Venue      |             Paper Title and Paper Interpretation             |                           Code                           |
| :----------------------------------------------------------: | :------------: | :----------------------------------------------------------: | :------------------------------------------------------: |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) | *Neurips '13*  |  MLDS Multilinear Dynamical Systems for Tensor Time Series   |         [mlds](https://github.com/lileicc/mlds)          |
|                        *Ph.D. Thesis*                        | *ProQuest '11* |      Fast Algorithms for Mining Co-evolving Time Series      |                           None                           |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *KDD '09*    | DynaMMo: Mining and Summarization of Coevolving Sequences with Missing Values 🌟 |      [dynammo](https://github.com/lileicc/dynammo)       |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *VLDB '10*   |      Parsimonious Linear Fingerprinting for Time Series      | [pliF](https://lileicc.github.io/software/plif-r345.zip) |

## [Feng Zhou](https://www.f-zhou.com/) (CMU)

|                             TYPE                             |    Venue    |             Paper Title and Paper Interpretation             |                  Code                  |
| :----------------------------------------------------------: | :---------: | :----------------------------------------------------------: | :------------------------------------: |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *TPAMI '12* | Hierarchical Aligned Cluster Analysis for Temporal Clustering of Human Motion 🌟 | [HACA](https://www.f-zhou.com/tc.html) |

## [Chun-Tung Li](https://scholar.google.com/citations?hl=en&user=9qsSJ-0AAAAJ&view_op=list_works&sortby=pubdate) (CUHK)

|                             TYPE                             |                Venue                |             Paper Title and Paper Interpretation             |                       Code                       |
| :----------------------------------------------------------: | :---------------------------------: | :----------------------------------------------------------: | :----------------------------------------------: |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) | *ACM Trans. Comput. Healthcare '20* | mSIMPAD: Efficient and Robust Mining of Successive Similar Patterns of Multiple Lengths in Time Series 🌟 | [mSIMPAD](https://github.com/chuntungli/mSIMPAD) |
|                        *Ph.D. Thesis*                        |           *ProQuest '21*            | Mobile sensing based human stress monitoring for smart health applications |                       None                       |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |           *IEEE MASS '21*           | Repetitive Activity Monitoring from Multivariate Time Series A Generic and Efficient Approach |                       None                       |

## [Tong Hanghang](https://scholar.google.com/citations?hl=en&user=RaINcuUAAAAJ&view_op=list_works&sortby=pubdate) (UIUC)

|                             TYPE                             |   Venue   |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :-------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) | *WWW '21* |                Network of Tensor Time Series                 |          [NET3](https://github.com/baoyujing/NET3)           |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) | *SDM '15* |      Fast Mining of a Network of Coevolving Time Series      | [dcmf](https://github.com/kokikwbt/dcmf/tree/master) (Unofficial) |
| ![spatio-temporal forecasting](https://img.shields.io/badge/-Tensor-blue) | *KDD '15* | Facets: Fast comprehensive mining of coevolving high-order time |  [facets](https://github.com/kokikwbt/facets) (Unofficial)   |

## Others

|                             TYPE                             |             Venue             |             Paper Title and Paper Interpretation             |                             Code                             |
| :----------------------------------------------------------: | :---------------------------: | :----------------------------------------------------------: | :----------------------------------------------------------: |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *Information Fusion '24*    | MultiBEATS Blocks of eigenvalues algorithm for multivariate time series dimensionality reduction 🌟 |  [MultiBEATS](https://github.com/auroragonzalez/multiBEATS)  |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *Information Sciences '24*   | Memetic segmentation based on variable lag aware for multivariate time series 🌟 |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          *TKDE '23*           | Change Point Detection in Multi-channel Time Series via a Time-invariant Representation 🌟 |      [MC-TIRE](https://github.com/caozhenxiang/MC-TIRE)      |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |           *TII '23*           | A Boundary Consistency-Aware Multitask Learning Framework for Joint Activity Segmentation and Recognition With Wearable Sensors | [Coming soom](https://github.com/xspc/Segmentation-Sensor-based-HAR) 🙃 |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          SIGMOD '23           | Time2State: An Unsupervised Framework for Inferring the Latent States in Time Series Data 🌟 |     [Time2State](https://github.com/Lab-ANT/Time2State)      |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *TKDD '23*           |        Modeling Regime Shifts in Multiple Time Series        |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |     *World Wide Web '23*      | Anomaly and change point detection for time series with concept drift |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *EAAI '23*           | PrecTime A deep learning architecture for precise time series segmentation in industrial manufacturing operations |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          *IMWUT '22*          | ColloSSL Collaborative Self-Supervised Learning for Human Activity Recognition 🌟 |     [collossl](https://github.com/akhilmathurs/collossl)     |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          *MSSP '22*           | A multivariate time series segmentation algorithm for analyzing the operating statuses of tunnel boring machines |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |      *Technometrics '22*      | Bayesian Hierarchical Model for Change Point Detection in Multivariate Sequences | [Supplementary Materials](https://doi.org/10.1080/00401706.2021.1927848) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |    *Neurips Workshop '22*     |                 Are uGLAD? Time will tell! 🌟                 |          [tGLAD](https://github.com/Harshs27/tGLAD)          |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *Applied Intelligence '22*   |  Change point detection for compositional multivariate data  |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *ICDM '22*           | Change Detection with Probabilistic Models on Persistence Diagrams |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *EAAI '22*           |   Graft : A graph based time series data mining framework    |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |        *GLOBECOM '22*         | Multi-level Contrast Network for Wearables-based Joint Activity Segmentation and Recognition |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |           ESWA '22            | Real-time Change-Point Detection A deep neural network-based adaptive approach for detecting changes in multivariate time series data |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |  *npj digital medicine '21*   |      U-Sleep: resilient high-frequency sleep staging 🌟       |              [website](https://sleep.ai.ku.dk/)              |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |        *IEEE TSP '21*         | Change Point Detection in Time Series Data Using Autoencoders With a Time-Invariant Representation 🌟 |           [TIRE](https://github.com/deryckt/TIRE)            |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *IJCNN '21*          | A Transferable Technique for Detecting and Localising Segments of Repeating Patterns in Time series |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *IOTJ '21*           | DeepSeg Deep-Learning-Based Activity Segmentation Framework for Activity Recognition Using WiFi |      [DeepSeg](https://github.com/ChunjingXiao/DeepSeg)      |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |  *Information Sciences '21*   | Change-point detection based on adjusted shape context method cost |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |        *IEEE TCYB '20*        | An Online Unsupervised Dynamic Window Method to Track Repeating Patterns From Sensor Data 🌟 |             [FingdingIOR](https://goo.gl/TV1seZ)             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) | *Pattern Recognit. Lett. '20* |     A new approach for optimal time-series segmentation      |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |           *SDM '20*           |      Lag-aware multivariate time-series segmentation 🌟       |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) | *Pattern Recognit. Lett. '20* | Memetic algorithm for multivariate time-series segmentation 🌟 |       [ma_mts](https://github.com/hlim-kist/ma _ mts)        |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |         *ICASSP '20*          |          Modeling Piece-Wise Stationary Time Series          |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |         *Neurips '19*         | U-Time: A Fully Convolutional Network for Time Series Segmentation Applied to Sleep Staging 🌟 |         [U-Time](https://github.com/perslev/U-Time)          |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |     *Neurocomputing '19*      | A hybrid dynamic exploitation barebones particle swarm optimisation algorithm for time series segmentation |   [tssa](https://github.com/ayrna/tssa?tab=readme-ov-file)   |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *TKDE '18*           | BEATS Blocks of Eigenvalues Algorithm for Time series Segmentation 🌟 |       [BEATS](https://github.com/auroragonzalez/BEATS)       |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |          *Arxiv '18*          | Time Series Segmentation through Automatic Feature Learning 🌟 |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) | *Applied Soft Computing '16*  | Change points detection in crime-related time series An on-line fuzzy approach based on a shape space representation |                             None                             |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |          *WACV '16*           | Decomposing Time Series with application to Temporal Segmentation 🌟 |   [Hog1D](https://github.com/AjeyPaiK/HOG1D) (Unofficial)    |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *J. Am. Stat. Assoc. '14*   | A Nonparametric Approach for Multiple Change Point Analysis of Multivariate Data |                             None                             |
| ![univariate time series forecasting](https://img.shields.io/badge/-Univariate-brightgreen) |     *Neural Networks '13*     | Change-point detection in time-series data by relative density-ratio estimation 🌟 | [RuLSIF](https://github.com/anewgithubname/change_detection?tab=readme-ov-file) |
