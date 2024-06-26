
## N2 Fit Function

A function to evaluate the resolving power (RP) given a nitrogen spectrum is available, but still experimental. The function is called `fit_n2` and it tries to automatically fit a nitrogen spectrum and calculate the resolving power. The function accepts two kinds of arguments to retrieve the correct spectra to fit, and the motor and detector name:

* A negative integer number, that refers to which scan you want to evaluate. To evaluate the last scan use
  ```python
  fit_n2(-1)
  ```
* A string representing the unique identifier of a scan
  ```python
  fit_n2('67548d')
  ```

The function tries to fit the spectra, and if it is successful it estimates the RP by calculating the Gaussian contribution to the FWHM of the N2:1s-->𝛑* transition via the fit parameters. Additionally, it calculates the 3rd peak to 1st valley ratio.

It is possible to pass the following arguments to the function:

* `motor`: `string`, the motor name to be used as x-axis
* `detector`: `string`, the detector readings to be used as y-axis
* `print_fit_report`: `boolean`, it will print the complete fit report from the package `lmfit`
* `save_img`: `string`, the absolute path and name to save the plots produced by the fit routine
* `fit`: `bool`, if False, disable the fit routine and plot only the data and the initial guess
* `fix_param`: `bool`, if True, the (&sigma;), gamma (&gamma;) and `skew` parameters of each peak are the same.

Additionally, it is possible to modify the initial guess of the parameters, see the section "Manual modification of the initial parameters for the fit". Here below is the function with all the possible arguments to copy/paste and modify:

```python
fit_n2(scan=-1, motor='pgm', detector='Keithley01', print_fit_report=False, save_img=False, fit=True, vc1='auto', amp_sf=6, sigma=0.02, sigma_min=0.001, sigma_max=0.02, gamma=0.055)
```

### Fit Function: SkewedVoigtModel

To perform the fit the python package [lmfit](https://lmfit.github.io/) is used. Two fit functions are available.

#### All Free Parameters

Pass the following parameter to the fit function:
```python
fix_param=False
```
The fit function is a sum of a straight line and ten Skewed Voigt Functions, as defined in the lmfit package, see the documentation about the [SkewedVoigtModel](https://lmfit.github.io/lmfit-py/builtin_models.html). Each function has five parameters: amplitude (A), center (&mu;), sigma (&sigma;), and gamma (&gamma;), as usual for a Voigt distribution, and adds a new parameter `skew`.

#### Fixed Parameters

Pass the following parameter to the fit function (or do not pass anything, by default it is true):
```python
fix_param=True
```
The fit function is a sum of a straight line and seven Skewed Voigt Functions, as defined in the lmfit package, see the documentation about the [SkewedVoigtModel](https://lmfit.github.io/lmfit-py/builtin_models.html). Each function has five parameters: amplitude (A), center (&mu;), sigma (&sigma;), and gamma (&gamma;), as usual for a Voigt distribution, and adds a new parameter `skew`. The sigma (&sigma;), gamma (&gamma;) and `skew` are the same for all the skewed Voigt functions.

### Automatic Guessing of the Initial Parameters for the Fit

The function tries to find out automatically the best initial parameters for the fit. First of all, normalization to the maximum value of the data is performed (this might create problems if we have an outlier with very high intensity).

* **Centers position (&mu;)**: The function looks for the maximum at the lowest energy in the spectra and assumes it is the first peak. The center of the other peaks is assumed using theoretical values for the peak separation. The fit procedure limits lower and upper bound to +/- 2.355\*2\*&sigma;.
* **Amplitudes (A)**: The amplitude is defined as the intensity of the data at the position of the centers (&mu;) and scaled by a factor of 6.
* **Sigma (&sigma;)**: 0.02 eV.
* **Gamma (&gamma;)**: 0.0563. This value determines the Lorentzian FWHM: FWHM_l = 2\*&gamma;.
* **`skew`**: 0.

### Manual Modification of the Initial Parameters for the Fit

A number of parameters can be modified by passing the following arguments to the function:

* `vc1`: The center of the first peak, can be a `float` or set to 'auto'
* `amp_sf`: Scaling factor for the amplitude, default is 6
* `sigma`: The sigma of the skewed Voigt functions, default 0.02
* `sigma_min`: The lower bound of the sigma value for the fit, default 0.001
* `sigma_max`: The upper bound of the sigma value for the fit, default 0.02
* `gamma`: The gamma parameter, default 0.055

### Estimation of the RP

For the fit, a sum of 11 skewed Voigt functions is assumed. Once the fit is performed the RP is calculated as the ratio of the center and the FWHM of the second peak:

```
RP = v2_&mu; / v2_fwhm_g
```

where the gaussian contribution to the FWHM is calculated as:
```
v2_fwhm_g = 2 * v2_&sigma; * sqrt(2 * ln(2))
```

The problem of this method is that some of the monochromator contributions escape into the Lorentzian shape, the &gamma; parameter. One would have to know the physical lifetime broadening of each line extremely precisely in order to determine the "artificial" monochromator contribution into &gamma;.

### 3rd-Peak to 1st-Valley Ratio

The fit routine returns the v/p ratio.

[Here](http://help.bessy.de/~follath/spektren/nitrogen/simulation.html) is a link to some of R. Ollath's calculations (accessible only within BESSY network).
