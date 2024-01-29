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
|                             TYPE                             |     Venue      |             Paper Title and Paper Interpretation             |                           Code                            |
| :----------------------------------------------------------: | :------------: | :----------------------------------------------------------: | :-------------------------------------------------------: |
|                        *Ph.D. Thesis*                        | *ProQuest '21* | Explainable and Network-Based Approaches for Decision-making in Emergency Management |                           None                            |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *CIKM '21*   |    Actionable Insights in Urban Multivariate Time-series     |        [RaTSS](https://github.com/AdityaLab/RaTSS)        |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *TIST '20*   | Cut-n-Reveal: Time-Series Segmentations with Explanations 🌟  | [Cut-n-Reveal](https://github.com/AdityaLab/Cut-n-Reveal) |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |   *AAAI '18*   |           Automatic Segmentation of Data Sequences           |                           None                            |
|                        *Ph.D. Thesis*                        | *ProQuest '18* |    Segmenting, Summarizing and Predicting Data Sequences     |                           None                            |
| ![multivariate time series forecasting](https://img.shields.io/badge/-Multivariate-red) |  *vt.edu '18*  |    Segmentations with Explanations for Outage Analysis 🌟     |                           None                            |
