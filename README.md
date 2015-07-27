BAED: Burned Area Extraction and Dating

CODE is stored in baed_code.r.

Section I: algorithm input and outputs

The primary inputs to the BAED algorithm are the normalized multidate Landsat Surface Reflectance products, and time series vegetation index. The BAED algorithm will calculate spectral index difference between successive landsat images, find “core burned” pixels, fit burn probability surface, and covert continuous burn probability surface to binary by maximizing detection rate, predicting expected NDVI range, and dating each burned patch.

The main output of BAED is a burned patch in vector format, with each patch was associated with disturbance time, cloud conditions, and confidence level. A patch without disturbance time found was assigned 0, and regarded as false detections, and therefore should be discarded. Confidence was classified as five levels with 1 as most confident and 5 as least confidence.

Section II: algorithm steps

The processing flow of BAED contains familiar steps of image selection, preprocessing, and determine burned cores, shaping burned patches, and date burned patches.

2.1 Image selection: The selected image can be from three generation of landsat data, inclduing landsat 5(Thematic Mapper), landsat 7 (Enhanced Thematic Mapper plus), and landsat 8 (Operational Land Imager). Six bands from visible, near infrared, shortwave infrared wavelength with 30 meter spatial resolution should be used. In image selection process, the absence of clouds was given highest priority, and then the consistency of seasonality.

2.2 Preprocessing: This includes geometric correction, atmospheric correction, and radiometric normalization. The geometric correction and atmospheric correction can be done through ESPA, but radiometric normalization will be a key step to correct the effects from sun-sensors-geometry, phonology, and random error. Iteratively re-weighted multivariate alteration detection (IR-MAD) technique can be used for this purpose, and it can be done in ENVI (Canty & Nielsen, 2008).

2.3 Determining burned cores: The aim of this step was to identify most possible burned pixels on the landscape to minimize the omission rate of potential burned patches. The determination of “core burned” pixel for burned land mapping was based on spectral indices change induced by fire between two landsat images. The spectral indices used in this step includes: (1) Normalized Difference Vegetation Index (NDVI), Normalized Burn Ratio (NBR), Disturbance Index (DI) due to their sensitivity to fire and wide application by previous studies.

2.4 Shaping the burned areas: This step aimed at growing “core burned” pixels into patches by analyzing the spectral properties of those pixels close to the seed pixels, and applied a region growth algorithm to enlarge the burned patches previously detected. The BAED first to estimate the probability that a pixel was burned based on a boosted regression tree models (BRT) was used. BRT model tend to produce higher binary classification accuracies that other machine-learning approaches (hastie and others 2009). The burn probability map was first converted to a binary burned/unburned image, and converted into polygons, and holes within burned patches were removed.

2.5 Dating the burned patches: The dataset used to date burned patches is based on the MOD13Q1 V5 product, which provides EVI and NDVI at 16 day composite using the maximum value composites (MVC) method at 250-meter spatial resolution. The basis of our method is comparison of expected NDVI with NDVI within burned patches to find change date. Expected NDVI described expected NDVI range in each year for undisturbed vegetation.
