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
  + [Clark, P., 1987,The cyclical component of U.S. economic activity, Quarterly Journal of Economics 102, 797-814](https://www.jstor.org/stable/1884282)
  + [Nelson, C. R., & Plosser, C. R. (1982). Trends and random walks in macroeconomic time series: Some evidence and implications. Journal of Monetary Economics, 10(2), 139–162.](https://www.sciencedirect.com/science/article/abs/pii/0304393282900125)
  + [Burmeister, Edwin, and Kent Wall, "Kalman Filtering Estimation of Unobserved Rational Expectations with an Application to the German Hyperinflation," Journal of Econometrics, November 1982, 20, 255–84.](https://www.sciencedirect.com/science/article/abs/pii/0304407682900215)
  + [Stock, J.H., Watson, M.W., 1991. A probability model of the coincident economic indicators. In: Moore, G., Lahiri, K. (Eds.), The Leading Economic Indicators: New Approaches and Forecasting Records. Cambridge
University Press, Cambridge, pp. 63–89](https://www.nber.org/papers/w2772)
  + [Kim, C.J., Nelson, C.R., 1989. The time-varying parameter model for modeling
changing conditional variance: the case of the Lucas hypothesis. J. Bus. Econ.
Stat. 7 (4), 433–440](https://www.tandfonline.com/doi/abs/10.1080/07350015.1989.10509755)

+ **The meaning of each DGP, DGP_Practice, Practice is as follows.**
  + DGP : Data Generating Process (Basic Code)
  + DGP_Practice : Data Generating Process for models to be run in Practice (Modified code for Application)
  + Practice : Modified code and Real Data(Include result_plot)
  
___

## ※ Index ※
___
## Space State Model
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

## Regime Switching Models
+ **Serially Uncorrelated Data and Switching.**
  + DGP
    + The case of Independent Switching
    + The case of Markov Switching
    
+ **Serially Correlated Data and Switching.**
  + DGP
    + Dealing with the Problem of Unobserved $S_t$

+ **Structure Break.**
  + DGP
    + 2Structure Break
    + 3Structure Break
    + Expected Duration of a Regime in a Markov-Switching Model
      
+ **Kim's Smoothing Algorithm.**
  + DGP
    + Markov-Switching Dependent Model with Kim's Smoothing
    + Markov-Switching Serially Correlated Model with Kim's Smoothing
___
## ※ Summary ※

## **Unobserved Component Model.**
### Model 1.

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

___
### Model 2.(Real_Estate)

$Y_t = N_t + X_t$ <br/>
$N_t = Mu +N_{t-1} + v_t,\quad v_t \sim iidN(0, \sigma^2_v)$<br/>
$X_t = \phi_1 * X_{t-1} + \phi_2 * X_{t-2} + \phi_3 * X_{t-3} + \phi_4 * X_{t-4} + e_t,\quad e_t \sim iidN(0, \sigma^2_e)$

$N_t$ = Trend Component<br/>
$X_t$ = Cyclical Component<br/>
$v_t$, $e_t$ = Independent white noise processes

**Measurement equation**<br/>
$Y_t = H*B_t$ <br/>

**Transition equation**<br/>
$B_t = Mu + F*B_{t-1} + u_t,\quad u_t \sim N(0,Q)$ <br/>

**SS Parameter**<br/>
$C$ = 0  <br/>
$H$ = [1 1 0 0 0]<br/>
$R$ = 0  <br/>
$Mu$ = [ $Mu$ ; 0 ; 0 ; 0 ; 0 ]<br/>
$F$ = [ 1 0 0 0 0 ; 0  $\phi_1\ \phi_2$ $\phi_3\ \phi_4$ ; 0 1 0 0 0 ; 0 0 1 0 0 ; 0 0 0 1 0 ] <br/>
$Q$ = [ $\sigma_v^2$ 0 0 0 0 ; 0 $\sigma_e^2$ 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 0 ; 0 0 0 0 0 ]<br/>
___
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

___
### Model 2 (Practice K.H_Kang(2012))

$Y_{1t} = \Lambda_1 * B_{t} + e_{1t},\quad e_{1t} \sim iidN(0, \sigma^2_{e1})$<br/>
$Y_{2t} = \Lambda_2 * B_{t} + e_{2t},\quad e_{2t} \sim iidN(0, \sigma^2_{e2})$<br/>
$Y_{3t} = \Lambda_3 * B_{t} + e_{3t},\quad e_{3t} \sim iidN(0, \sigma^2_{e3})$<br/>
$Y_{4t} = \Lambda_4 * B_{t} + e_{4t},\quad e_{4t} \sim iidN(0, \sigma^2_{e4})$<br/>
$Y_{5t} = \Lambda_5 * B_{t} + e_{5t},\quad e_{5t} \sim iidN(0, \sigma^2_{e5})$<br/>
$Y_{6t} = \Lambda_6 * B_{t} + e_{6t},\quad e_{6t} \sim iidN(0, \sigma^2_{e6})$<br/>
$Y_{7t} = \Lambda_7 * B_{t} + e_{7t},\quad e_{7t} \sim iidN(0, \sigma^2_{e7})$<br/>
$Y_{8t} = \Lambda_8 * B_{t} + e_{8t},\quad e_{8t} \sim iidN(0, \sigma^2_{e8})$<br/>
$Y_{9t} = \Lambda_9 * B_{t} + e_{9t},\quad e_{9t} \sim iidN(0, \sigma^2_{e9})$<br/>
$Y_{10t} = \Lambda_10 * B_{t} + e_{10t},\quad e_{10t} \sim iidN(0, \sigma^2_{e10})$<br/>

$$\Lambda =
\begin{bmatrix}
1 & \frac{1-e^{-\tau_1\\lambda}}{\tau_1\\lambda} & \frac{1-e^{-\tau_1\\lambda}}{\tau_1\\lambda} - e^{-\tau_1\\lambda}\\
1 & \frac{1-e^{-\tau_2\\lambda}}{\tau_2\\lambda} & \frac{1-e^{-\tau_2\\lambda}}{\tau_2\\lambda} - e^{-\tau_2\\lambda}\\
\vdots  & \vdots  &  \vdots \\
1 & \frac{1-e^{-\tau_N\\lambda}}{\tau_N\\lambda} & \frac{1-e^{-\tau_N\\lambda}}{\tau_N\\lambda}- e^{-\tau_N\\lambda}
\end{bmatrix}$$

$\tau = [3, 6, 9, 12, 18, 24, 30, 36, 60, 120]$  <br/>
$\lambda = 0.15$ <br/>

**Measurement equation**<br/>
$Y_t = C + \Lambda*B_t + E_t,\quad E_t \sim iidN(0, R)$ <br/>

**Transition equation**<br/>
$B_t = [bL_t, bS_t, bC_t]\quad Level, Slope, Curvature$ <br/>
$B_t = Mu + F*B_{t-1} + u_t,\quad u_t \sim N(0,Q)$ <br/>

**SS Parameter**<br/>
$C$ = [0 0 0 0 0 0 0 0 0 0]'  <br/>
$H$ = $\Lambda$<br/>
$R$ = diag($\sigma^2_{e1}\,\sigma^2_{e2}\,\sigma^2_{e3}\,\sigma^2_{e4}\,\sigma^2_{e5}\,\sigma^2_{e6}\,\sigma^2_{e7}\,\sigma^2_{e8}\,\sigma^2_{e9}\,\sigma^2_{e10}\$)<br/>
$Mu$ = [ $Mu1$ ; $Mu2$ ; $Mu3$ ]<br/>
$F$ = [ $\phi_1\$ 0 0 ; 0 $\phi_2$ 0 ; 0 0 $\phi_3$ ]  <br/>
$Q$ = [ $\sigma_v^1$ $Cov_1$ $Cov_2$ ; $Cov_1$ $\sigma_v^2$ $Cov_3$ ; $Cov_2$ $Cov_3$ $\sigma_v^3$ ]<br/>
___
___
## **Time Varying Parameter Model.**
### Model

$I_t = B0_t + B1_t * CPI_{t-1} + B2_t * IAIP_{t-1} + B3_t * I_{t-1} + e_t,\quad e_t \sim iidN(0, \sigma^2_e)$<br/>

$B0_{t} = B0_{t-1} + v_{0t},\quad v_{0t} \sim iidN(0, \sigma^2_{v0})$<br/>
$B1_{t} = B1_{t-1} + v_{1t},\quad v_{1t} \sim iidN(0, \sigma^2_{v1})$<br/>
$B2_{t} = B2_{t-1} + v_{2t},\quad v_{2t} \sim iidN(0, \sigma^2_{v2})$<br/>
$B3_{t} = B3_{t-1} + v_{3t},\quad v_{3t} \sim iidN(0, \sigma^2_{v3})$<br/>

$CPI$, $IAIP$, $I$ are exogenuous variable <br/>

**Measurement equation**<br/>
$Y_t = H*B_t + E_t,\quad E_t \sim iidN(0, R)$ <br/>

**Transition equation**<br/>
$B_{it} = B_{it-1} + u_{it},\quad u_{it} \sim N(0,Q)\quad i = 0,1,2,3$<br/>

**SS Parameter**<br/>
$C$ = 0  <br/>
$H$ = [1  $CPI_{t-1}$  $IAIP_{t-1}$  $I_{t-1}$]<br/>
$R$ = $\sigma^2_{e}$<br/>
$Mu$ = [0 0 0 0]'<br/>
$F$ = eye(4)  <br/>
$Q$ = [ $\sigma_v^0$ 0  0  0 ; 0 $\sigma_v^1$ 0  0 ; 0  0  $\sigma_v^2$  0 ; 0  0  0  $\sigma_v^3$]<br/>

