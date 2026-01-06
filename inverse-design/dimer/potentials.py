import jax
from jax import config
import jax.numpy as jnp




config.update("jax_enable_x64", True)


def smoothing(r, ron, rcut):
    # Ensure rcut and ron are properly separated
    rcut = jnp.clip(rcut, 1e-6, None)  # Prevent rcut from being too small
    ron = jnp.clip(ron, 1e-6, rcut - 1e-6)  # Ensure ron < rcut

    # Ensure valid distance for the potential
    return jnp.where(
        r < ron,
        1,
        jnp.where(
            r > rcut,
            0,
            jnp.clip((rcut**2 - r**2) ** 2 * (rcut**2 + 2 * r**2 - 3 * ron**2), a_min=0)
            / ((rcut**2 - ron**2) ** 3),
        ),
    )


def morse(r, rmin, rmax, D0, alpha, r0):
    return jnp.where(
        r >= rmax,
        0.0,
        D0 * (jnp.exp(-2 * alpha * (r - r0)) - 2 * jnp.exp(-alpha * (r - r0))),
    )


def morse_x(r, rmin, rmax, D0, alpha, r0, ron):
    """Cut Morse at rmax to match LAMMPS morse pair style."""
    return morse(r, rmin, rmax, D0, alpha, r0)


def morse_x_repulsive(r, rmin, rmax, D0, alpha, r0, ron):
    # FIXME: isn't this just -morse_x(r, rmin, rmax, D0, alpha, r0, ron)
    return -morse(r, rmin, rmax, D0, alpha, r0) * smoothing(r, ron, rmax)


"""
def repulsive(r, rmin, rmax, A, alpha):
    epsilon = 1e-6 
    base = jnp.maximum(rmax - r, epsilon)  

    return jnp.where(r < rmax, (A / (alpha * rmax)) * base**alpha, 1e-10)
"""


def lj_cut(r, epsilon, sigma, rcut):
    """Lennard-Jones 12-6 truncated (no shifting) to mirror LAMMPS lj/cut."""
    r_safe = jnp.maximum(r, 1e-12)
    sr = sigma / r_safe
    sr6 = sr ** 6
    sr12 = sr6 ** 2
    energy = 4.0 * epsilon * (sr12 - sr6)
    return jnp.where(r < rcut, energy, 0.0)


def repulsive(r, epsilon, sigma, rcut):
    """Alias for lj_cut for compatibility with existing call sites."""
    return lj_cut(r, epsilon=epsilon, sigma=sigma, rcut=rcut)
