#!/usr/bin/env python
##############################################################################
#
# (c) 2025 Simon Billinge.
# All rights reserved.
#
# File coded by: Caden Myers, Simon Billinge, and members of the Billinge
#                group.
#
# See GitHub contributions for a more detailed list of contributors.
# https://github.com/diffpy/diffpy.cmipdf/graphs/contributors
#
# See LICENSE.rst for license information.
#
##############################################################################
"""Form factors (characteristic functions) used in PDF nanoshape
fitting.

These are used to calculate the attenuation of the PDF due to a finite
size. For a crystal-like nanoparticle, one can calculate the PDF via
Gnano(r) = f(r) Gcryst(r), where f(r) is the nanoparticle characteristic
function and Gcryst(f) is the crystal PDF.

These functions are meant to be imported and added to a FitContribution
using the 'register_function' method of that class.
"""

__all__ = [
    "spherical_particle",
    "spheroidal_particle",
    "lognormal_spherical_distribution",
    "sheet_particle",
    "shellCF",
    "shellCF2",
    "SASCF",
]

import numpy
from numpy import arctan as atan
from numpy import arctanh as atanh
from numpy import ceil, exp, log, log2, pi, sign, sqrt
from numpy.fft import fftfreq, ifft
from scipy.special import erf

from diffpy.srfit.fitbase.calculator import Calculator


def spherical_particle(radial_dist, p_diameter):
    """Spherical nanoparticle characteristic function.

    Parameters
    ----------
    radial_dist : float or array-like
        The distance of interaction.
    p_diameter : float
        The particle diameter.


    From Kodama et al., Acta Cryst. A, 62, 444-453
    (converted from radius to diameter)
    """
    f = numpy.zeros(numpy.shape(radial_dist), dtype=float)
    if p_diameter > 0:
        x = numpy.array(radial_dist, dtype=float) / p_diameter
        inside = x < 1.0
        xin = x[inside]
        f[inside] = 1.0 - 1.5 * xin + 0.5 * xin * xin * xin
    return f


def spheroidal_particle(radial_dist, r_equatorial, r_polar):
    """Spheroidal characteristic function specified using radii.

    Spheroid with radii (r_equatorial, r_equatorial, r_polar)

    Parameters
    ----------
    radial_dist : float or array-like
        The distance of interaction.
    r_polar : float
        The polar radius of the spheroid.
    r_equatorial
        The equatorial radius of the spheroid.


    Note
    ----
    - `r_equatorial < r_polar` equates to a prolate spheroid
    - `r_equatorial > r_polar` equates to a oblate spheroid
    - `r_equatorial == r_polar` is a sphere
    """
    d_equatorial = 2.0 * r_equatorial
    pelpt = 1.0 * r_polar / r_equatorial
    return _calculate_spheroidal_cf(radial_dist, d_equatorial, pelpt)


def _calculate_spheroidal_cf(r, d_equatorial, axis_ratio):
    """Calculate the spheroidal nanoparticle characteristic function.

    Form factor for ellipsoid with radii
    (d_equatorial/2, d_equatorial/2, axis_ratio*d_equatorial/2)

    Parameters
    ----------
    r : float or array-like
        The distance of interaction
    d_equatorial : float
        The equatorial diameter
    axis_ratio : float
        The ratio of axis lengths


    From Lei et al., Phys. Rev. B, 80, 024118 (2009)
    """
    pelpt = 1.0 * axis_ratio

    if d_equatorial <= 0 or pelpt <= 0:
        return numpy.zeros_like(r)

    # to simplify the equations
    v = pelpt
    d = 1.0 * d_equatorial
    d2 = d * d
    v2 = v * v

    if v == 1:
        return spherical_particle(r, d_equatorial)

    rx = r
    if v < 1:

        r = rx[rx <= v * d_equatorial]
        r2 = r * r
        f1 = (
            1
            - 3 * r / (4 * d * v) * (1 - r2 / (4 * d2) * (1 + 2.0 / (3 * v2)))
            - 3
            * r
            / (4 * d)
            * (1 - r2 / (4 * d2))
            * v
            / sqrt(1 - v2)
            * atanh(sqrt(1 - v2))
        )

        r = rx[numpy.logical_and(rx > v * d_equatorial, rx <= d_equatorial)]
        r2 = r * r
        f2 = (
            (
                3 * d / (8 * r) * (1 + r2 / (2 * d2)) * sqrt(1 - r2 / d2)
                - 3
                * r
                / (4 * d)
                * (1 - r2 / (4 * d2))
                * atanh(sqrt(1 - r2 / d2))
            )
            * v
            / sqrt(1 - v2)
        )

        r = rx[rx > d_equatorial]
        f3 = numpy.zeros_like(r)

        f = numpy.concatenate((f1, f2, f3))

    elif v > 1:

        r = rx[rx <= d_equatorial]
        r2 = r * r
        f1 = (
            1
            - 3 * r / (4 * d * v) * (1 - r2 / (4 * d2) * (1 + 2.0 / (3 * v2)))
            - 3
            * r
            / (4 * d)
            * (1 - r2 / (4 * d2))
            * v
            / sqrt(v2 - 1)
            * atan(sqrt(v2 - 1))
        )

        r = rx[numpy.logical_and(rx > d_equatorial, rx <= v * d_equatorial)]
        r2 = r * r
        f2 = (
            1
            - 3 * r / (4 * d * v) * (1 - r2 / (4 * d2) * (1 + 2.0 / (3 * v2)))
            - 3.0
            / 8
            * (1 + r2 / (2 * d2))
            * sqrt(1 - d2 / r2)
            * v
            / sqrt(v2 - 1)
            - 3
            * r
            / (4 * d)
            * (1 - r2 / (4 * d2))
            * v
            / sqrt(v2 - 1)
            * (atan(sqrt(v2 - 1)) - atan(sqrt(r2 / d2 - 1)))
        )

        r = rx[rx > v * d_equatorial]
        f3 = numpy.zeros_like(r)

        f = numpy.concatenate((f1, f2, f3))

    return f


def lognormal_spherical_distribution(radial_dist, p_diameter, p_sigma):
    """Spherical nanoparticle characteristic function with lognormal
    size distribution.

    Parameters
    ----------
    radial_dist : float or array-like
        The distance of interaction.
    p_diameter : float
        The mean particle diameter.
    p_sigma : float
        The log-normal width of the particle diameter.


    Here, radial_dist is the independent variable, mu is the mean of the
    distribution
    (not of the particle size), and s is the width of the distribution.
    This is
    the characteristic function for the lognormal distribution of particle
    diameter:

    F(r, mu, s) = 0.5*Erfc((-mu-3*s^2+Log(r))/(sqrt(2)*s))
               + 0.25*r^3*Erfc((-mu+Log(r))/(sqrt(2)*s))*exp(-3*mu-4.5*s^2)
               - 0.75*r*Erfc((-mu-2*s^2+Log(r))/(sqrt(2)*s))*exp(-mu-2.5*s^2)

    The expectation value of the distribution gives the average particle
    diameter, p_diameter. The variance of the distribution gives p_sigma^2.
    mu and s can be expressed in terms of these as:

    s^2 = log((p_sigma/p_diameter)^2 + 1)
    mu = log(p_diameter) - s^2/2

    Source unknown
    """
    if p_diameter <= 0:
        return numpy.zeros_like(radial_dist)
    if p_sigma <= 0:
        return spherical_particle(radial_dist, p_diameter)

    sqrt2 = sqrt(2.0)
    s = sqrt(log(p_sigma * p_sigma / (1.0 * p_diameter * p_diameter) + 1))
    mu = log(p_diameter) - s * s / 2
    if mu < 0:
        return numpy.zeros_like(radial_dist)

    return (
        0.5 * erfc((-mu - 3 * s * s + log(radial_dist)) / (sqrt2 * s))
        + 0.25
        * radial_dist
        * radial_dist
        * radial_dist
        * erfc((-mu + log(radial_dist)) / (sqrt2 * s))
        * exp(-3 * mu - 4.5 * s * s)
        - 0.75
        * radial_dist
        * erfc((-mu - 2 * s * s + log(radial_dist)) / (sqrt2 * s))
        * exp(-mu - 2.5 * s * s)
    )


def sheet_particle(r, thickness):
    """Nanosheet characteristic function.

    Parameters
    ----------
    r: float or array-like
        The distance of interaction.
    thickness : float
        The thickness of nanosheet.


    From Kodama et al., Acta Cryst. A, 62, 444-453
    """
    # handle zero or negative thickness.  make it work for scalars and arrays.
    if thickness <= 0:
        return 0 * thickness
    # process scalar r
    if numpy.isscalar(r):
        rv = 1 - 0.5 * r / thickness if r < thickness else 0.5 * thickness / r
        return rv
    # handle array-type r
    ra = numpy.asarray(r)
    lo = ra < thickness
    hi = ~lo
    f = numpy.empty_like(ra, dtype=float)
    f[lo] = 1 - 0.5 * ra[lo] / thickness
    f[hi] = 0.5 * thickness / ra[hi]
    return f


def shellCF(r, radius, thickness):
    """Spherical shell characteristic function.

    Parameters
    ----------
    radius
        Inner radius
    thickness
        Thickness of shell


    outer radius = radius + thickness

    From Lei et al., Phys. Rev. B, 80, 024118 (2009)
    """
    d = 1.0 * thickness
    a = 1.0 * radius + d / 2.0
    return shellCF2(r, a, d)


def shellCF2(r, a, delta):
    """Spherical shell characteristic function.

    Parameters
    ----------
    a
        Central radius
    delta
        Thickness of shell

    outer radius = a + thickness/2


    From Lei et al., Phys. Rev. B, 80, 024118 (2009)
    """
    a = 1.0 * a
    d = 1.0 * delta
    a2 = a**2
    d2 = d**2
    dmr = d - r
    dmr2 = dmr**2

    f = (
        r
        * (
            16 * a * a2
            + 12 * a * d * dmr
            + 36 * a2 * (2 * d - r)
            + 3 * dmr2 * (2 * d + r)
        )
        + 2 * dmr2 * (r * (2 * d + r) - 12 * a2) * sign(dmr)
        - 2 * (2 * a - r) ** 2 * (r * (4 * a + r) - 3 * d2) * sign(2 * a - r)
        + r * (4 * a - 2 * d + r) * (2 * a - d - r) ** 2 * sign(2 * a - d - r)
    )

    f[r > 2 * a + d] = 0

    den = 8.0 * r * d * (12 * a2 + d2)
    zmask = den == 0.0
    vmask = ~zmask
    f[vmask] /= den[vmask]
    f[zmask] = 1
    return f


class SASCF(Calculator):
    """Calculator class for characteristic functions from sas-models.

    This class wraps a sas.models.BaseModel to calculate I(Q) related to
    nanoparticle shape. This I(Q) is inverted to f(r) according to:
    f(r) = 1 / (4 pi r) * SINFT(I(Q)),
    where "SINFT" represents the sine Fourier transform.

    Attributes
    ----------
    _model
        BaseModel object this adapts.
    Managed Parameters
        These depend on the parameters of the BaseModel object held by _model.
        They are created from the 'params' attribute of the BaseModel. If a
        dispersion is set for the BaseModel, the dispersion "width" will be
        accessible under "<parname>_width", where <parname> is the name a
        parameter adjusted by dispersion.
    """

    def __init__(self, name, model):
        """Initialize the generator.

        Parameters
        ----------
        name
            A name for the SASCF
        model
            SASModel object this adapts.
        """
        Calculator.__init__(self, name)

        self._model = model

        from diffpy.srfit.sas.sasparameter import SASParameter

        # Wrap normal parameters
        for parname in model.params:
            par = SASParameter(parname, model)
            self.addParameter(par)

        # Wrap dispersion parameters
        for parname in model.dispersion:
            name = parname + "_width"
            parname += ".width"
            par = SASParameter(name, model, parname)
            self.addParameter(par)

        return

    def __call__(self, r):
        """Calculate the characteristic function from the transform of
        the BaseModel."""

        # Determine q-values.
        # We want very fine r-spacing so we can properly normalize f(r). This
        # equates to having a large qmax so that the Fourier transform is
        # finely spaced. We also want the calculation to be fast, so we pick
        # qmax such that the number of q-points is a power of 2. This allows us
        # to use the fft.
        #
        # The initial dr is somewhat arbitrary, but using dr = 0.01 allows for
        # the f(r) calculated from a particle of diameter 50, over r =
        # arange(1, 60, 0.1) to agree with the spherical_particle with
        # Rw < 1e-4%.
        #
        # We also have to make a q-spacing small enough to compute out to at
        # least the size of the signal.
        raise NotImplementedError(
            "As of release 3.2.0, SAS characteristic functions are not working"
            + " but we hope to have them working again in a future release."
        )

        dr = min(0.01, r[1] - r[0])
        ed = 2 * self._model.calculate_ER()

        # Check for nans. If we find any, then return zeros.
        if numpy.isnan(ed).any():
            y = numpy.zeros_like(r)
            return y

        rmax = max(ed, 2 * r[-1])
        dq = pi / rmax
        qmax = pi / dr
        numpoints = int(2 ** (ceil(log2(qmax / dq))))
        qmax = dq * numpoints

        # Calculate F(q) = q * I(q) from model
        q = fftfreq(int(qmax / dq)) * qmax
        fq = q * self._model.evalDistribution(q)

        # Calculate g(r) and the effective r-points
        rp = fftfreq(numpoints) * 2 * pi / dq
        # Note sine transform = imaginary part of ifft
        gr = ifft(fq).imag

        # Calculate full-fr for normalization
        assert rp[0] == 0.0
        frp = numpy.zeros_like(gr)
        frp[1:] = gr[1:] / rp[1:]

        # Inerpolate onto requested grid, do not use data after jump in rp
        assert numpoints % 2 == 0
        nhalf = numpoints / 2
        fr = numpy.interp(r, rp[:nhalf], gr[:nhalf])
        vmask = r != 0
        fr[vmask] /= r[vmask]

        # Normalize. We approximate fr[0] by using the fact that f(r) is linear
        # at low r. By definition, fr[0] should equal 1.
        fr0 = 2 * frp[2] - frp[1]
        fr /= fr0

        # Fix potential divide-by-zero issue, fr is 1 at r == 0
        fr[~vmask] = 1

        return fr


def erfc(x):
    return 1.0 - erf(x)


# End of file
