#!/usr/bin/env python3

import numpy as np
import jax.numpy as jnp
import jax
from jax.test_util import check_grads
#import scipy
from jax import config
config.update("jax_enable_x64", True)

def r_random(seed, xyzmin, xyzmax, n):
    if output_format == "code":
        print(f"    // random r-vectors: seed = {seed}, xyzmin = {xyzmin}, xyzmax = {xyzmax}, n = {n}")
    np.random.seed(seed)
    return np.array(np.random.uniform(xyzmin, xyzmax, 3*n).reshape(n, 3))

# Works, but is not used.
#def spherical_harmonics_numpy(r, lmax):
#    x2 = r[0] * r[0]
#    y2 = r[1] * r[1]
#    z2 = r[2] * r[2]
#    rn = np.sqrt(x2 + y2 + z2)
#    theta = np.arccos(r[2] / rn)
#    phi = np.sign(r[1]) * np.arccos(r[0] / np.sqrt(x2 + y2))
#    if phi < 0:
#        phi = 2 * np.pi + phi
#    result = np.array([])
#    for l in range(0, lmax+1):
#        m_list = np.arange(0, l+1)
#        ylm_complex = jax.scipy.special.sph_harm(m_list, np.array([l]), phi, theta)
#        # Convert to real spherical harmonics.
#        for m in np.arange(-l, l+1):
#            if m < 0:
#                ylm_real = np.sqrt(2) * pow(-1.0, m) * np.imag(ylm_complex[abs(m)])
#            elif m == 0:
#                ylm_real = np.real(ylm_complex[0])
#            else:
#                ylm_real = np.sqrt(2) * pow(-1.0, m) * np.real(ylm_complex[abs(m)])
#            result = np.append(result, ylm_real)
#    return result

def spherical_harmonics_single(r, l, m):
    x2 = r[0] * r[0]
    y2 = r[1] * r[1]
    z2 = r[2] * r[2]
    rn = jnp.sqrt(x2 + y2 + z2)
    theta = jnp.arccos(r[2] / rn)
    phi = jnp.sign(r[1]) * jnp.arccos(r[0] / jnp.sqrt(x2 + y2))
    if phi < 0:
        phi = 2 * jnp.pi + phi
    ylm_complex = jax.scipy.special.sph_harm(jnp.array([abs(m)]), jnp.array([l]), phi, theta)
    if m < 0:
        return (jnp.sqrt(2) * pow(-1.0, m) * jnp.imag(ylm_complex))[0]
    elif m == 0:
        return (jnp.real(ylm_complex))[0]
    else:
        return (jnp.sqrt(2) * pow(-1.0, m) * jnp.real(ylm_complex))[0]

def spherical_harmonics(r, lmax):
    # Seems JAX gets confused if this input array is anything but a raw Python list.
    r = r.tolist()
    result = []
    for l in range(0, lmax+1):
        for m in jnp.arange(-l, l+1):
            result.append(spherical_harmonics_single(r, l, m))
    return jnp.array(result)

def spherical_harmonics_deriv(r, lmax):
    r = r.tolist()
    grad_f = jax.grad(spherical_harmonics_single)
    result = []
    for l in range(0, lmax+1):
        for m in jnp.arange(-l, l+1):
            result.append([float(x) for x in grad_f(r, l, m)])
    return jnp.array(result)

def test_data(r, f, deriv, df=None, delta=1.0E-8, **kwargs):
    fr = jnp.array([f(ir, **kwargs) for ir in r])
    if deriv:
        # User provided derivative function.
        if df is not None:
            if output_format == "code":
                print("    // using analytic gradient (user-defined function)")
            dfr = jnp.array([df(ir, **kwargs) for ir in r])
    else:
        dfr = np.array([[[0.0, 0.0, 0.0] for ilm in range(0, len(fr[0]))] for ir in r])
    if output_format == "code":
        print(f"    // function parameters =  {kwargs}")
        print("    // clang-format off")
    for ri, fri, dfri in zip(r, fr, dfr):
        if output_format == "code":
            #print("    " + " ".join([f"tc->r.push_back({x:23.16E});" for x in ri]))
            rnorm = np.linalg.norm(ri)
            print("    " + f"tc->rnorm.push_back({rnorm:23.16E}); " + " ".join([f"tc->rhat.push_back({x/rnorm:23.16E});" for x in ri]))
            for frij, dfrij in zip(fri, dfri):
                #print(f"    tc->fr.push_back({frij:23.16E});", end='')
                print(f"    tc->ylm.push_back({frij:23.16E});", end='')
                if deriv:
                    #print(" " + " ".join([f"tc->dfr.push_back({x:23.16E});" for x in dfrij]), end=' ')
                    print(" " + " ".join([f"tc->ylmd.push_back({x:23.16E});" for x in dfrij]), end=' ')
                print("")
        else:
            # Not implemented.
            pass
    if output_format == "code":
        print("    // clang-format on")
    return

if __name__ == '__main__':
    output_format="code"
    #test_data(r_random(1, -1.5, 1.5, 5), spherical_harmonics, False, lmax=2)
    #test_data(r_random(1, -1.5, 1.5, 5), spherical_harmonics, True, df=spherical_harmonics_deriv, lmax=3)
    #test_data(r_random(123, -1.5, 1.5, 5), spherical_harmonics, True, df=spherical_harmonics_deriv, lmax=5)
    #test_data(r_random(932, -1.0, 1.0, 5), spherical_harmonics, True, df=spherical_harmonics_deriv, lmax=8)
    test_data(r_random(456, -0.9, 0.9, 5), spherical_harmonics, True, df=spherical_harmonics_deriv, lmax=4)
