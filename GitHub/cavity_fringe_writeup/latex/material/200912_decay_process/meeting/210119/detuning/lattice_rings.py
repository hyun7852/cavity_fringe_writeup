#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-
# Time-stamp: "2020-10-09 17:19:25 sb"

#  file       lattice_rings.py
#  copyright  (c) Sebastian Blatt 2020

import numpy as np
import scipy.special
import srlab.plot.pyplot as plt
import srlab.units as u


AUPOL = 4 * u.pi * u.epsilon0 * u.bohrradius ** 3
SCALAR_POLARIZABILITY_1S0_813 = 278 * AUPOL
SCALAR_POLARIZABILITY_1S0_914 = 253 * AUPOL
SCALAR_POLARIZABILITY_3P0_813 = SCALAR_POLARIZABILITY_1S0_813
SCALAR_POLARIZABILITY_3P0_914 = 215 * AUPOL


# 1D lattice math uses results from the Mathieu equation for an infinitely
# extended cos^2 kx potential. In this approximation, we can get all
# information about the bands (except for their resolution into quasi-momenta)
# directly from special functions. Note that the Mathieu q parameter is simply
# the trap depth in lattice recoils divided by four.


def lattice_mathieu_q_from_trap_depth_recoil(Ut):
    return Ut / 4


def lattice_band_edge_lower(n, q):
    return scipy.special.mathieu_a(n, q) - 2 * q


def lattice_band_edge_upper(n, q):
    return scipy.special.mathieu_b(n + 1, q) - 2 * q


def lattice_band_center(n, q):
    lo = lattice_band_edge_lower(n, q)
    hi = lattice_band_edge_upper(n, q)
    return (lo + hi) / 2


def lattice_trapped_bands_harmonic_approximation(q):
    """If we approximate the lattice as harmonic with a trap frequency
    given by the ground state, we can ask how many harmonic oscillator
    levels would fit into this trap.
    """

    return int(np.round(np.sqrt(q)))


def lattice_trapped_bands(q):
    """If we want to know how many states actually are supported by the
    lattice (note that we omit gravity here), we have to find the
    highest band that still fits completely into the lattice.
    """
    for n in range(int(np.round(q))):
        if lattice_band_edge_upper(n, q) >= 0.0:
            return n
    return -1


def has_units(func):
    """A decorator that does NOTHING, but indicates that the result of
    the wrapped function is dimensionful using the srlab.units module.
    """

    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)

    return wrapper


@has_units
def recoil_frequency(mass, wavelength):
    return u.h / (2 * mass * wavelength ** 2)


@has_units
def lattice_trap_depth_from_trap_frequency(trap_frequency, recoil_frequency):
    return u.h * trap_frequency ** 2 / (4 * recoil_frequency)


@has_units
def lattice_parameters_from_trap_frequency(trap_depth, recoil_frequency):
    Erec = u.h * recoil_frequency
    q = lattice_mathieu_q_from_trap_depth_recoil(float(trap_depth / Erec))
    n_max = lattice_trapped_bands(q) - 1
    return (q, n_max)


def lattice_1d_sideband_energy_range(n, delta_n, qg, qe):
    """Calculate the detuning of the n -> n + delta_n sideband when

    transitioning from a ground state 1D lattice described by Mathieu
    parameter qg to an excited state 1D lattice described by Mathieu
    parameter qe.

    Returns a tuple of (mean energy g, mean energy e, minimal detuning,
    maximal detuning) in units of the lattice recoil frequency.
    """

    Eg_lo = lattice_band_edge_lower(n, qg)
    Eg_hi = lattice_band_edge_upper(n, qg)
    Ee_lo = lattice_band_edge_lower(n + delta_n, qe)
    Ee_hi = lattice_band_edge_upper(n + delta_n, qe)
    return ((Eg_lo + Eg_hi) / 2, (Ee_lo + Ee_hi) / 2, Ee_lo - Eg_hi, Ee_hi - Eg_lo)


def lattice_1d_sideband_list(delta_n, qg, qe):
    rc = []
    for n in range(lattice_trapped_bands(qg)):
        Eg, Ee, delta_lo, delta_hi = lattice_1d_sideband_energy_range(
            n, delta_n, qg, qe
        )
        rc.append(((n, delta_n), Eg, Ee, delta_lo, delta_hi))
    return rc


def lattice_2d_sideband_energy_range(n, m, delta_n, delta_m, qg1, qg2, qe1, qe2):
    """Calculate the detuning of the (n, m) -> (n + delta_n, m +
    delta_m) sideband when transitioning from a ground state
    noninterfering 2D lattice described by Mathieu parameters (qg1, qg2)
    to an excited state noninterfering 2D lattice described by Mathieu
    parameter (qe1, qe2).

    Returns a tuple of (mean energy g, mean energy e, minimal detuning,
    maximal detuning) in units of the lattice recoil frequency.
    """
    # dE_zero = (
    #     lattice_band_center(0, qe1)
    #     + lattice_band_center(0, qe2)
    #     - lattice_band_center(0, qg1)
    #     - lattice_band_center(0, qg2)
    # )

    Eg_lo = lattice_band_edge_lower(n, qg1) + lattice_band_edge_lower(m, qg2)
    Eg_hi = lattice_band_edge_upper(n, qg1) + lattice_band_edge_upper(m, qg2)
    Ee_lo = lattice_band_edge_lower(n + delta_n, qe1) + lattice_band_edge_lower(
        m + delta_m, qe2
    )
    Ee_hi = lattice_band_edge_upper(n + delta_n, qe1) + lattice_band_edge_upper(
        m + delta_m, qe2
    )
    return ((Eg_lo + Eg_hi) / 2, (Ee_lo + Ee_hi) / 2, Ee_lo - Eg_hi, Ee_hi - Eg_lo)


def lattice_2d_sideband_list(delta_n, delta_m, qg1, qg2, qe1, qe2):
    rc = []
    for n in range(lattice_trapped_bands(qg1)):
        for m in range(lattice_trapped_bands(qg2)):
            Eg, Ee, delta_lo, delta_hi = lattice_2d_sideband_energy_range(
                n, m, delta_n, delta_m, qg1, qg2, qe1, qe2
            )
            rc.append(((n, m, delta_n, delta_m), Eg, Ee, delta_lo, delta_hi))
    return rc


def boltzmann_factors_unnormalized(energies, energy_zero, temperature):
    """Get a list of Boltzmann factors the list of `energies`. Subtract
    `energy_zero` from each energy in the list and assume `temperature`.
    Returns a list of *unnormalized* Boltzmann factors.
    """
    return np.exp(-(energies - energy_zero) / temperature)


def boltzmann_factors(energies, energy_zero, temperature):
    """Get a list of Boltzmann factors the list of `energies`. Subtract
    `energy_zero` from each energy in the list and assume `temperature`.
    Returns a list of normalized Boltzmann factors assuming that the
    specified energies are all there are in the system (truncated Boltzmann
    distribution).
    """
    bs = boltzmann_factors_unnormalized(energies, energy_zero, temperature)
    return bs / np.sum(bs)


def boltzmann_factors_from_transition_list(transition_list, temperature):
    energies = []
    for _, Eg, _, _, _ in transition_list:
        energies.append(Eg)
    energies = np.array(energies)
    return boltzmann_factors(energies, np.min(energies), temperature)


def lorentz_lineshape(delta, delta0, fwhm):
    """Peak-normalized Lorentzian lineshape as a function of detuning `delta`,
    centered at `delta0` with full-width at half maximum `fwhm`.
    """
    x = 2 * (delta - delta0) / fwhm
    return 1 / (1 + x ** 2)


def linear_lorentz_response(detuning, fwhm, transition_list, boltzmann_factors):
    """Sum over all transitions in `transition_list` and calculate the
    linear Lorentzian response assuming detuning `detuning` and the same
    full-width at half maximum `fwhm` for each transition. The
    amplitudes of each transition are summed according to the normalized
    Boltzmann factors that match the order of transition_list.
    """

    rc = 0.0
    for (b, data) in zip(boltzmann_factors, transition_list):
        _, _, _, delta_lo, delta_hi = data
        rc += b * lorentz_lineshape(detuning, (delta_lo + delta_hi) / 2, fwhm)
    return rc


def plot_spectrum(
    ax,
    nu_rec,
    transition_list,
    boltzmann_factors,
    clock_laser_fwhm,
    detuning_min,
    detuning_max,
    npoints=1000,
):
    """Assuming linear response, plot an excitation spectrum for a fixed lattice depth."""
    deltas = np.linspace(
        float(detuning_min / nu_rec), float(detuning_max / nu_rec), npoints
    )
    response = np.zeros(deltas.shape)
    for j, delta in enumerate(deltas):
        response[j] = linear_lorentz_response(
            delta, clock_laser_fwhm, transition_list, boltzmann_factors
        )
    ax.plot(deltas * nu_rec.in_units_of("kHz"), response)

    ax.set_xlabel("Detuning from free space (kHz)")
    ax.set_ylabel("Linear response")
    return ax


def gaussian_distribution(xs, ys, sigma_x, sigma_y):
    Xs, Ys = np.meshgrid(xs, ys)
    return np.exp(-(Xs ** 2) / (2 * sigma_x ** 2) - (Ys ** 2) / (2 * sigma_y ** 2))


def lattice_distribution(xs, ys, wx, wy, qg1, qg2):
    Xs, Ys = np.meshgrid(xs, ys)
    Ux = qg1 * np.exp(-2 * (Xs ** 2) / wx ** 2)
    Uy = qg2 * np.exp(-2 * (Ys ** 2) / wy ** 2)
    rc = Ux + Uy
    return rc / np.max(np.max(rc))


def show_image(ax, image, with_contours=True):
    ax.imshow(
        image,
        aspect=1,
        origin="upper",
        interpolation="bilinear",
        cmap="Blues",
        vmin=0,
        vmax=1,
    )
    if with_contours:
        cs = ax.contour(image, origin="upper", colors="k")  # , levels=[0.1, 0.5, 0.95])
        ax.clabel(cs, inline=1)
    return ax


if __name__ == "__main__":

    mass = 88 * u.amu
    wavelength = u.Real("914 nm")

    nu_rec = recoil_frequency(mass, wavelength)
    E_rec = u.h * 

    nu_trap = u.Real("110 kHz")
    trap_depth_g = lattice_trap_depth_from_trap_frequency(nu_trap, nu_rec)
    eta = float(SCALAR_POLARIZABILITY_3P0_914 / SCALAR_POLARIZABILITY_1S0_914)
    # eta = 1
    trap_depth_e = eta * trap_depth_g

    temperature = float(u.kB * u.Real("0.1 uK") / E_rec)
    clock_laser_fwhm = float(u.Real("5 kHz") / nu_rec)

    # These are the lattice parameters (Mathieu-q, trapped bands) in the
    # deepest part of the lattice.
    q_g, n_g_max = lattice_parameters_from_trap_frequency(trap_depth_g, nu_rec)
    q_e, n_e_max = lattice_parameters_from_trap_frequency(trap_depth_e, nu_rec)

    print("q_g",q_g,n_g_max)
    print("q_e",q_g,n_g_max)

    plt.set_defaults(style="default", use_latex=False)
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=[297 / 25.4, 210 / 25.4])
    (ax1, ax2), (ax3, ax4) = axs

    lattice_waist = u.Real("456 um")
    pixel_size = u.Real("5.4 um")
    image_width = 150
    image_height = 150
    atom_sigma_x = 20 * pixel_size
    atom_sigma_y = 30 * pixel_size

    xs = np.arange(-image_width // 2, image_width // 2 - 1) * pixel_size.in_units_of(
        "um"
    )
    ys = np.arange(-image_height // 2, image_height // 2 - 1) * pixel_size.in_units_of(
        "um"
    )

    atom_image = gaussian_distribution(
        xs, ys, atom_sigma_x.in_units_of("um"), atom_sigma_y.in_units_of("um")
    )
    lattice_image = lattice_distribution(
        xs,
        ys,
        lattice_waist.in_units_of("um"),
        lattice_waist.in_units_of("um"),
        q_g,
        q_g,
    )
    ax1.set_title("Initial atomic distribution")
    show_image(ax1, atom_image)
    ax2.set_title("Lattice intensity distribution")
    show_image(ax2, lattice_image)

    sz = atom_image.shape
    linear_response_image = np.zeros(sz)
    for i in range(sz[0]):
        print(f"{i/sz[0]*1e2:.1f}%")
        for j in range(sz[1]):
            qg = lattice_image[i, j] * q_g
            qe = qg * eta
            transition_list = lattice_2d_sideband_list(0, 0, qg, qg, qe, qe)

            bs = boltzmann_factors_from_transition_list(transition_list, temperature)
            linear_response_image[i, j] = linear_lorentz_response(
                100, clock_laser_fwhm, transition_list, bs
            )
    # print("transition list",transition_list)
    ax3.set_title("Linear response image")
    show_image(ax3, linear_response_image, with_contours=False)
    ax4.set_title("Final atom distribution")
    show_image(
        ax4,
        (1 - linear_response_image / np.max(np.max(linear_response_image)))
        * atom_image,
        with_contours=False,
    )

    plt.save_pdf(fig, "lattice_rings.pdf")

# lattice_rings.py ends here
