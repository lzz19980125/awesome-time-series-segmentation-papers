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
