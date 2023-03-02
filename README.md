## DGP and Application of Space State Model(Code with Data) ##
## ※ Notice ※
___
+ **Textbook**
  + State-Space Models with Regime Switching -Classical and Gibbs-Sampling Approaches with Applications
  + By Chang-Jin Kim and Daniel C R. Halbert 
  + Publisher: The MIT Press

+ **Code**
  + This code is based on the work of Professor Young-min Kim, an assistant professor at Jeonbuk National University.

+ **Data**
  + Some data were referred to Professor Chang-jin Kim's data, and others were obtained from KOSIS.

+ **The meaning of each DGP, DGP_Practice, Practice is as follows.**
  + DGP : Data Generating Process (Basic Code)
  + DGP_Practice : Data Generating Process for models to be run in Practice (Modified code for Application)
  + Practice : Modified code and Real Data(Include result_plot)
___

## ※ Index ※
___
+ **Unobserved Component Model.**
  + DGP
    + UC_DGP
  + DGP_Practice
    + UC_DGP_Bivariate_Practice_Version
    + UC_DGP_Uivariate_Practice_Version
  + Practice
    + UC_Uivariate_Add_Period
    + UC_Uivariate_Kims_Data
    
+ **Dynamic Common Factor Model.**
  + DGP
    + DCF_DGP
  + DGP_Practice
    + DCF_DGP_JeonBuk_Practice_Version
    + DCF_DGP_KH-Kang(2012)_Practice_Version
  + Practice
    + DCF_K.H_Kang(2012)
    + DCF_Jeonbuk

+ **Time Varying Parameter Model.**
  + DGP
    + TVP_DGP
  + DGP_Practice
    + TVP_DGP_Korea_Practice_Version
  + Practice
    + TVP_Korea
    + TVP_Korea_Del_Constant
___
## ※ Summary ※
___
## **Unobserved Component Model.**
### Model

$Y_t = N_t + X_t$ <br/>
$N_t = Mu +N_{t-1} + v_t,\quad v_t \sim iidN(0, \sigma^2_v)$<br/>
$X_t = \phi_1 * X_{t-1} + \phi_2 * X_{t-2} + e_t,\quad e_t \sim iidN(0, \sigma^2_e)$

$N_t$ = Trend Component<br/>
$X_t$ = Cyclical Component<br/>
$v_t$, $e_t$ = Independent white noise processes
## Result
![image](https://user-images.githubusercontent.com/109870987/222392656-319789e0-330a-45ef-a1cb-b40f3235ca02.png)

+ Decomposition of trend and cyclical elements in Real GDP.
+ There were cases where the model failed to find a trend and the results showed two cyclical factors.
  + In this case, the problem was solved by removing the constraint of $\sigma^2_v$, the variance of the trend component, or by inducing a low value.

![image](https://user-images.githubusercontent.com/109870987/222394744-53d9b8aa-7c47-41cd-a96d-f44c1edc1b04.png)
+ The result was derived by setting the constraints as shown in the picture in **"paramconst"**, which is a sheet that gives constraints.
