"""Collection of functions not tied to a particular object."""
from __future__ import annotations

import typing

import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.stats


@typing.no_type_check
def exponential(domain, scale, decay):
    """Exponential function for fitting/plotting."""
    return scale * np.exp(-decay * domain)


@typing.no_type_check
def odr_exponential(parameters, domain):
    """Exponential function for odr estimates."""
    return parameters[0] * np.exp(-parameters[1] * domain)


def sample_parameters(
    random: np.random.Generator,
    a_fit: float,
    b_fit: float,
    covariance: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """
    function to, where possible, use maximum likelihood estimates of a and b
    parameter plus estimated covariance matrix to directly sample from the
    probability distribution of a and b (assume Gaussian with mean of max
    likelihood estimates and given covariance structure)
    """
    return random.multivariate_normal([a_fit, b_fit], covariance, size=200)


def find_confidence_interval(
    random: np.random.Generator,
    distances: npt.NDArray[np.float64],
    sampled_parameters: npt.NDArray[np.float64],
    bit_threshold: pd.DataFrame | None = None,
    index: pd.Index | None = None,
) -> pd.DataFrame:
    """Gives an empirical estimate of the prediction interval.
    function to take each of the sampled a, b values and use them to sample
    directly from the distribution of scores taking into account the Gaussian
    noise (a function of distance, a, b)."""
    # species x samples
    point_estimates = exponential(
        distances,
        sampled_parameters[:, 0],
        sampled_parameters[:, 1],
    )

    # species x samples
    exp = np.exp(-1 * sampled_parameters[:, 1] * distances)
    variance = sampled_parameters[:, 0] * (1 - exp) * exp

    # if variance is negative, need to ignore those samples
    invalid = variance <= 0
    # set to 0 to handle issues with sqrt and random normal
    variance[invalid] = 0

    # with normal, can't broadcast a shape on matrix inputs, instead
    # perform a standard normal draw and do rescaling manually
    # species x samples x draws
    draws = (
        random.standard_normal(size=(*point_estimates.shape, 200))
        * np.sqrt(variance[..., None])
        + point_estimates[..., None]
    )

    # set invalid values to nan, keep first entry as point estimate
    # this works because the invalid entries have 0 variance.
    draws[invalid, 1:] = np.nan

    # now find mean and std, ignoring nans
    mean = np.nanmean(draws, axis=(1, 2))  # type: ignore[no-untyped-call]
    std = np.nanstd(draws, axis=(1, 2))  # type: ignore[no-untyped-call]

    low, high = scipy.stats.norm.interval(0.99, mean, std)

    # calculate p values analytically from std estimate
    if bit_threshold is not None:
        p_values = scipy.stats.norm.cdf(bit_threshold["bit_threshold"], mean, std)

        return pd.DataFrame(
            {
                "p_values": p_values,
                "low_interval": low,
                "high_interval": high,
            },
            index=index,
        )

    return pd.DataFrame(
        {
            "low_interval": low,
            "high_interval": high,
        },
        index=index,
    )
