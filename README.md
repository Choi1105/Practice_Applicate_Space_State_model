## DGP and Application of Space State Model(Code with Data) ##
## ※ Notice ※
___
+ **Textbook**
  + [State-Space Models with Regime Switching -Classical and Gibbs-Sampling Approaches with Applications](https://mitpress.mit.edu/9780262535502/state-space-models-with-regime-switching)
  + By Chang-Jin Kim and Daniel C R. Halbert 
  + Publisher: The MIT Press

+ **Code**
  + This code is based on the work of Professor [Young-min Kim](https://kimymecon.weebly.com/matlab.html), an assistant professor at Jeonbuk National University.

+ **Data**
  + Some data were referred to Professor [Chang-jin Kim's data](http://econ.korea.ac.kr/~cjkim/), and others were obtained from [KOSIS](https://kosis.kr/statHtml/statHtml.do?orgId=214&tblId=DT_214N_Z01900&vw_cd=MT_ZTITLE&list_id=214_21404&seqNo=&lang_mode=ko&language=kor&obj_var_id=&itm_id=&conn_path=MT_ZTITLE).

+ **Reference**
  + [Diebold, F. X., and C. L. Li, 2006, Forecasting the term structure of government bond yields,Journal of Econometrics130, pp. 337–364.](https://www.sas.upenn.edu/~fdiebold/papers/paper49/Diebold-Li.pdf)
  + [Kang, K. H., 2012a, Forecasting the term structure of Korean Government bond tields usingthe dynamic Nelson-Siegel class models,Asia-Pacific Journal of Financial Studies41, pp.765–787.](https://onlinelibrary.wiley.com/doi/full/10.1111/ajfs.12000)

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

**Measurement equation**<br/>
$Y_t = H*B_t$ <br/>

**Transition equation**<br/>
$B_t = Mu + F*B_{t-1} + u_t,\quad u_t \sim N(0,Q)$ <br/>

**SS Parameter**<br/>
$C$ = 0  <br/>
$H$ = [1 1 0]<br/>
$R$ = 0  <br/>
$Mu$ = [ $Mu$ ; 0 ; 0 ]<br/>
$F$ = [ 1 0 0 ; 0  $\phi_1\ \phi_2$ ; 0 1 0 ] <br/>
$Q$ = [ $\sigma_v^2$ 0 0 ; 0 $\sigma_e^2$ 0 ; 0 0 0 ]<br/>

## Result
![image](https://user-images.githubusercontent.com/109870987/222392656-319789e0-330a-45ef-a1cb-b40f3235ca02.png)

+ Decomposition of trend and cyclical component in Real GDP.
+ There were cases where the model failed to find a trend component and the results showed two cyclical factors.
  + In this case, the problem was solved by removing the constraint of $\sigma^2_v$, the variance of the trend component, or by inducing a low value.
  + The reason this constraint is acceptable is that trends should have lower variance than cycles.
  + The result was derived by setting the constraints as shown in the picture in **"paramconst"**, which is a sheet that gives constraints.
  
    ![image](https://user-images.githubusercontent.com/109870987/222394744-53d9b8aa-7c47-41cd-a96d-f44c1edc1b04.png)
___
## **Dynamic Common Factor Model.**
### Model 1 (Practice Jeonbuk)

$Y_{1t} = \Gamma_1 * C_{t} + e_{1t},\quad e_{1t} \sim iidN(0, \sigma^2_{e1})$<br/>
$Y_{2t} = \Gamma_2 * C_{t} + e_{2t},\quad e_{2t} \sim iidN(0, \sigma^2_{e2})$<br/>
$Y_{3t} = \Gamma_3 * C_{t} + e_{3t},\quad e_{3t} \sim iidN(0, \sigma^2_{e3})$<br/>
$Y_{4t} = \Gamma_4 * C_{t} + e_{4t},\quad e_{4t} \sim iidN(0, \sigma^2_{e4})$<br/>
$Y_{5t} = \Gamma_5 * C_{t} + e_{5t},\quad e_{5t} \sim iidN(0, \sigma^2_{e5})$<br/>
$Y_{6t} = \Gamma_6 * C_{t} + e_{6t},\quad e_{6t} \sim iidN(0, \sigma^2_{e6})$<br/>
$Y_{7t} = \Gamma_7 * C_{t} + e_{7t},\quad e_{7t} \sim iidN(0, \sigma^2_{e7})$<br/>
<br/>
$where\ \Gamma_1 = 1$<br/>
NO correlation between $e_1t, ... ,e_7t$<br/>

$C_t = Mu + \phi_1 * C_{t-1} + \phi_2 * C_{t-2} + v_t,\quad v_t \sim iidN(0, \sigma^2_v)$<br/>

**Measurement equation**<br/>
$Y_t = C + H*B_t + E_t,\quad E_t \sim iidN(0, R)$ <br/>

**Transition equation**<br/>
$B_t = Mu + F*B_{t-1} + u_t,\quad u_t \sim N(0,Q)$ <br/>

**SS Parameter**<br/>
$C$ = [0 0 0 0 0 0 0]'  <br/>
$H$ = [1 $\Gamma_2\  \Gamma_3\ \Gamma_4\ \Gamma_5\ \Gamma_6\ \Gamma_7\$]'<br/>
$R$ = diag($\sigma^2_{e1}\,\sigma^2_{e2}\,\sigma^2_{e3}\,\sigma^2_{e4}\,\sigma^2_{e5}\,\sigma^2_{e6}\,\sigma^2_{e7}\$)<br/>
$Mu$ = [ $Mu$ ; 0 ; 0 ]<br/>
$F$ = [ $\phi_1\ \phi_2$ 0 ; 1 0 0 ; 0 1 0 ]  <br/>
$Q$ = [ $\sigma_v^2$ 0 0 ; 0 0 0 ; 0 0 0 ]<br/>

## Result
![image](https://user-images.githubusercontent.com/109870987/222414639-c7491f32-85cc-4a86-982d-7414593534d0.png)

+ The result in the first row shows the estimated coincident composite index (Jeonbuk).<br/>
+ The result in the second row is the actual published coincident composite index (Jeonbuk).<br/>
+ The result in row 2, column 2 is the cyclical fluctuation of the coincident composite index, which the trend is removed from the result in row 2, column 1. <br/>
+ If you check the **Log_1Diff_result_Plot**, which is the result of estimation with the log-first-difference data, you can see that the estimation result and the picture of 2nd row, 2nd column are the same. <br/>


