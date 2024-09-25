
from __future__ import print_function
import numpy as np
import astropy.constants as const

def compute_snr(telescope, camera, exposure_in_hours, ab_magnitudes):

    diff_limit = 1.22*(500.*0.000000001)*206264.8062/telescope.aperture
    #print 'diff_limit', diff_limit

    # source_magnitudes = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]) * ab_magnitudes

    sky_brightness = np.array([23.807, 25.517, 22.627, 22.307, 21.917, 22.257, 21.757, 21.567, 22.417, 22.537])
    fwhm_psf = 1.22 * camera.pivotwave * 0.000000001 * 206264.8062 / telescope.aperture
    fwhm_psf[fwhm_psf < diff_limit] = fwhm_psf[fwhm_psf < diff_limit] * 0.0 + diff_limit

    sn_box = np.round(3. * fwhm_psf / camera.pixel_size)

    number_of_exposures = np.array([3., 3., 3., 3., 3., 3., 3., 3., 3., 3.])

    desired_exposure_time = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]) * exposure_in_hours * 3600.

    # time_per_exposure = desired_exposure_time / number_of_exposures

    signal_counts = camera.total_qe * desired_exposure_time * camera.ab_zeropoint * camera.aperture_correction * (np.pi) / 4. * \
        (telescope.aperture  * 100.0)**2 * camera.derived_bandpass * 10.**(-0.4*ab_magnitudes)

    # shot_noise_in_signal = signal_counts ** 0.5

    sky_counts = camera.total_qe * desired_exposure_time * camera.ab_zeropoint * np.pi / 4. * (telescope.aperture  * 100.0)**2 * \
        camera.derived_bandpass * 10.**(-0.4*sky_brightness) * (camera.pixel_size * sn_box)**2

    # shot_noise_in_sky = sky_counts ** 0.5

    read_noise = camera.detector_read_noise * sn_box * number_of_exposures**0.5

    dark_noise = sn_box * (camera.dark_current * desired_exposure_time)**0.5

    thermal_counts = desired_exposure_time * C_thermal(telescope, camera.pivotwave, camera.bandpass_r, sn_box**2, camera.pixel_size, \
                                                       camera.total_qe)

    snr = signal_counts / (signal_counts + sky_counts + read_noise**2 + dark_noise**2 + thermal_counts)**0.5

    #print
    #print
    #print 'aperture', telescope.aperture, ' meters'
    #print 'temperature', telescope.temperature, ' K'
    #print 'source_mag', source_magnitudes
    #print 'pivot waves', camera.pivotwave
    #print 'pixels', camera.pixel_size, ' arcsec'
    #print 'ab zeropoints', camera.ab_zeropoint
    #print 'total_qe', camera.total_qe
    #print 'ap corr', camera.aperture_correction
    #print 'bandpass R', camera.bandpass_r
    #print 'derived_bandpass', camera.derived_bandpass
    #print 'read_noise', camera.detector_read_noise
    #print 'dark rate', camera.dark_current
    #print 'sky_brightness', sky_brightness
    #print 'fwhm_psf', fwhm_psf
    #print 'sn_box', sn_box
    #print 'number_of_exposures', number_of_exposures
    #print 'detector_read_noise', camera.detector_read_noise
    #print 'time_per_exp', time_per_exposure
    #print 'signal_counts', signal_counts
    #print 'shot_noise_in_signal', shot_noise_in_signal
    #print 'sky_counts', sky_counts
    #print 'shot_noise_in_sky', shot_noise_in_sky
    #print 'read noise', read_noise
    #print 'dark noise', dark_noise
    #print 'thermal counts ', thermal_counts
    #print 'SNR', snr, np.max(snr)
    #print
    #print

    return snr


def C_thermal(telescope, wavelength_in_nm, R_bandwidth, N_pix, pixel_scale, QE_tot):

    # Needs inputs:
    #   - telescope object
    #   - R value - "resolution" of bandwidth
    #   - QE_tot = pixel QE
    #   - N_pix = number of pixels in photometric extraction box. NOTE: this should be the square of the 1D size
    #                   that is, if you intend a 3x3 box, use N_pix = 9 here.
    #   - pixel_scale = pixel size in arcsec

    wavelength_in_cm = wavelength_in_nm * 1e-7 # convert nm to cm
    bandwidth = wavelength_in_cm / R_bandwidth # still in cm

    h = const.h.to('erg s').value # Planck's constant erg s
    c = const.c.to('cm/s').value # speed of light [cm / s]

    energy_per_photon = h * c / wavelength_in_cm

    D = telescope.aperture * 100. # telescope diameter in cm

    Omega = 2.3504e-11 * pixel_scale**2 * N_pix

    #print 'Planck: ', planck(wavelength_in_cm, telescope.temperature)
    #print 'QE * Planck / E_per_phot', QE_tot * planck(wavelength_in_cm, telescope.temperature) / energy_per_photon
    #print 'E_per_phot', energy_per_photon
    #print 'Omega', Omega

    C_thermal = telescope.ota_emissivity * planck(wavelength_in_cm, telescope.temperature) / energy_per_photon * \
                        (np.pi / 4. * D**2) * QE_tot * Omega * bandwidth

    return C_thermal




def planck(wave,temp):

    h = const.h.to('erg s').value # Planck's constant erg s
    c = const.c.to('cm/s').value # speed of light [cm / s]

    k=1.3806488e-16 # Boltzmann's constant [erg deg K^-1]

    x = 2. * h * c**2 / wave**5

    exponent = (h * c / (wave * k * temp))

    PLANCK = x / (np.exp(exponent)-1.)

    return PLANCK


