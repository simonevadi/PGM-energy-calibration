# Procedure for Absolute Energy Calibration of Plane Grating Monochromators

**Objective:**
To calibrate the energy scale of plane grating monochromators (PGMs) by determining and correcting for angular offsets in the grating angles using measurements at different fixed focus constant (\(c_{ff}\)) values.

**Materials Needed:**
1. Plane grating monochromator (PGM)
2. Beamline capable of changing \(c_{ff}\)
3. Angular encoders with high precision
4. Known gas absorption spectra for calibration
5. Data acquisition system

**Steps:**

1. **Initial Setup:**
   - Ensure the monochromator is operational and connected to the beamline.
   - Verify the functionality of angular encoders and ensure they can read angles with an accuracy of better than 200 per 360°.

2. **Baseline Measurement:**
   - Perform initial energy calibration using the classical method, such as zero order checks and known gas absorption spectra analysis, to establish a reference.

3. **Collect Spectra at Different \(c_{ff}\) Values:**
   - Set the PGM to a specific \(c_{ff}\) value.
   - Record the spectrum of a known gas (e.g., \(N_2\), He, Ar) using the beamline.
   - Repeat the measurements at multiple \(c_{ff}\) values to gather a comprehensive set of spectra.

4. **Identify Energy Shifts:**
   - Compare the spectra obtained at different \(c_{ff}\) values.
   - Note any energy shifts (typically 1-2%) in the spectra, which indicate misalignments in the grating angles.

5. **Apply Fitting Procedure:**
   - Use the grating equation and the following formulas to calculate the true energy diffracted by the monochromator considering angular offsets:
     \[
     Nkl = 2 \cos Y \sin(Y + b) \quad \text{and} \quad Y = \frac{a + b}{2}
     \]
     \[
     \sin a = \frac{lNk}{c_{ff}^2 - 1} \left( \sqrt{c_{ff}^2 + (c_{ff}^2 - 1)^2 (lNk)^2} - 1 \right)
     \]
     \[
     b = -\arccos(c_{ff} \cos a)
     \]
     \[
     l_s = \frac{2 \cos(Y + \Delta Y) \sin(Y + \Delta Y + b + \Delta b)}{Nk}
     \]
   - Perform a least squares fit to the shifted energy positions, using the spectra obtained at different \(c_{ff}\) values to determine angular corrections \(\Delta Y\) and \(\Delta b\).

6. **Insert Angular Offsets:**
   - Adjust the PGM setup by applying the determined angular corrections \(\Delta Y\) and \(\Delta b\).
   - Re-measure the spectra at different \(c_{ff}\) values to verify the corrections.

7. **Validate Calibration:**
   - Perform zero order checks and compare the re-measured spectra with known gas absorption spectra.
   - Ensure that the residual energy scale changes are less than a few \(10^{-4}\), indicating successful calibration.

8. **Documentation and Reporting:**
   - Document the calibration process, including all measurements, calculated corrections, and final verification results.
   - Summarize the calibration results in a table, listing beamline, gas used, grating line density, angular corrections, and residual energy variations.

**Conclusion:**
By following this procedure, you can achieve a precise and accurate energy calibration for PGMs, ensuring that the energy scale is correct across the entire operational range of the monochromator. This method is quicker and can be performed without prior knowledge of the exact energy of the absorption features, making it versatile for various beamline applications.
