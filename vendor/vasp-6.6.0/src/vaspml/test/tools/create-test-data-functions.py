#!/usr/bin/env python3

import numpy as np
import jax.numpy as jnp
import jax
from jax.test_util import check_grads
import scipy
from jax import config
config.update("jax_enable_x64", True)

def x_random(seed, xmin, xmax, n):
    if output_format == "code":
        print(f"    // random x-values: seed = {seed}, xmin = {xmin}, xmax = {xmax}, n = {n}")
    np.random.seed(seed)
    return np.random.uniform(xmin, xmax, n)

def x_linspace(xmin, xmax, n):
    if output_format == "code":
        print(f"    // regular grid x-values: xmin = {xmin}, xmax = {xmax}, n = {n}")
    return np.linspace(xmin, xmax, n)

def x_geomspace(xmin, xmax, n):
    if output_format == "code":
        print(f"    // logarithmic grid x-values: xmin = {xmin}, xmax = {xmax}, n = {n}")
    return np.geomspace(xmin, xmax, n)

def gaussian(x, width):
    return jnp.exp(-width * x**2)

def cutoff_BP(x, rcut):
    if x > rcut: return 0.0
    if x < 0: return 1.0
    return 0.5 * (jnp.cos( x / rcut * jnp.pi) + 1)

def cutoff_MIWA(x, rcut):
    if x > rcut: return 0.0
    if x < 0: return 1.0
    y = 4.0 * x / rcut - 3.0
    if y < -1.0: return 1.0
    if y > 1.0: return 0.0
    return 0.25 * (y**3 - 3.0 * y + 2.0)

def cutoff_POLY2(x, rcut):
    if x > rcut: return 0.0
    if x < 0: return 1.0
    x /= rcut
    return ((15 - 6 * x ) * x - 10) * x**3 + 1

def cutoff_BUMP(x, rcut):
    if x > rcut: return 0.0
    if x < 0: return 1.0
    x /= rcut
    return jnp.exp(1 - 1 / (1 - x**2))

def spherical_bessel(x, n):
    return scipy.special.spherical_jn(n, x)

def spherical_bessel_deriv(x, n):
    return scipy.special.spherical_jn(n, x, derivative=True)

def modified_spherical_bessel(x, n):
    return scipy.special.spherical_in(n, x)

def modified_spherical_bessel_deriv(x, n):
    return scipy.special.spherical_in(n, x, derivative=True)

def exponential_modified_spherical_bessel(x, n, rij, sigma):
    a2 = rij**2 / (2 * sigma**2)
    ab = rij * x / sigma**2
    b2 = x**2 / (2 * sigma**2)
    return jnp.exp(-a2 - b2) * scipy.special.spherical_in(n, ab)

def test_data(x, f, deriv, df=None, delta=1.0E-8, **kwargs):
    fx = np.array([f(ix, **kwargs) for ix in x])
    if deriv:
        try:
            grad_f = jax.grad(f)
            dfx = np.array([grad_f(ix, **kwargs) for ix in x])
            if output_format == "code":
                print(f"    // using analytic gradient (auto-differentiation)")
        except:
            # User provided derivative function.
            if df is not None:
                if output_format == "code":
                    print("    // using analytic gradient (user-defined function)")
                dfx = np.array([df(ix, **kwargs) for ix in x])
            # Numeric derivative with given or default delta.
            elif np.ndim(delta) == 0:
                if output_format == "code":
                    print(f"    // using numeric gradient (central difference) with delta = {delta}")
                fh = [f(ix+delta, **kwargs) for ix in x]
                fl = [f(ix-delta, **kwargs) for ix in x]
                dfx = [(fhi - fli) / (2 * delta) for fhi, fli in zip(fh, fl)]
            # Numeric derivative delta search mode.
            else:
                dfx = []
                for i, idelta in enumerate(delta):
                    fh = [f(ix+idelta, **kwargs) for ix in x]
                    fl = [f(ix-idelta, **kwargs) for ix in x]
                    dfx.append([(fhi - fli) / (2 * idelta) for fhi, fli in zip(fh, fl)])
                with open("numeric-derivatives.dat", "w") as fout:
                    fout.write("# plot for [i=2:*] 'numeric-derivatives.dat' u 1:i w lp lw 2 title 'x-value '.(i-1)\n")
                    print("# plot for [i=2:*] 'numeric-derivatives.dat' u 1:i w lp lw 2 title 'x-value '.(i-1)")
                    dfx = np.transpose(dfx)
                    dfxmin = np.min(dfx, axis=1)
                    dfxspan = np.max(dfx, axis=1) - dfxmin
                    for i in range(len(dfx)):
                        dfx[i] = (dfx[i] - dfxmin[i]) / dfxspan[i]
                    dfx = np.transpose(dfx)
                    for i, idelta in enumerate(delta):
                        fout.write(f"{np.log10(idelta):23.16E} " + " ".join(f"{dfxij:23.16E}" for dfxij in dfx[i]) + "\n")
                return
    else:
        dfx = np.array([0.0 for ix in x])
    if output_format == "code":
        print(f"    // function parameters =  {kwargs}")
        print("    // clang-format off")
    for xi, fxi, dfxi in zip(x, fx, dfx):
        if output_format == "code":
            print(f"    tc->x.push_back({xi:23.16E}); tc->fx.push_back({fxi:23.16E});", end='')
        else:
            print(f"{xi:23.16E} {fxi:23.16E}", end='')
        if deriv:
            if output_format == "code":
                print(f" tc->dfx.push_back({dfxi:23.16E});", end='')
            else:
                print(f" {dfxi:23.16E}", end='')
        print("")
    if output_format == "code":
        print("    // clang-format on")
    return

if __name__ == '__main__':
    output_format="code"
    #test_data(x_random(1, -10, 10, 20), gaussian, True, width=0.1)
    #test_data(x_random(2, -0.1, 0.1, 10), gaussian, True, width=3.0)
    #test_data(x_linspace(-0.5, 6.5, 11), cutoff_BP, True, rcut=6.0)
    #test_data(x_linspace(-0.5, 6.5, 11), cutoff_MIWA, True, rcut=6.0)
    #test_data(x_linspace(-0.5, 6.5, 11), cutoff_POLY2, True, rcut=6.0)
    #test_data(x_linspace(-0.5, 6.5, 11), cutoff_BUMP, True, rcut=6.0)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, delta=np.logspace(-12, -4, num=45), n=0)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, delta=5.0E-8, n=0)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, df=spherical_bessel_deriv, n=0)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, df=spherical_bessel_deriv, n=1)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, df=spherical_bessel_deriv, n=2)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, df=spherical_bessel_deriv, n=3)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, df=spherical_bessel_deriv, n=4)
    #test_data(np.append(x_geomspace(0.01, 0.99, 9), x_linspace(1.0, 21, 11)), spherical_bessel, True, df=spherical_bessel_deriv, n=5)
    #test_data(x_linspace(0.0, 6.0, 13), modified_spherical_bessel, True, df=modified_spherical_bessel_deriv, n=0)
    test_data(x_linspace(0.000001, 6.0, 13), modified_spherical_bessel, True, df=modified_spherical_bessel_deriv, n=1)
    #test_data(x_linspace(0.0, 6.0, 13), modified_spherical_bessel, True, df=modified_spherical_bessel_deriv, n=2)
    #test_data(x_linspace(0.0, 6.0, 13), modified_spherical_bessel, True, df=modified_spherical_bessel_deriv, n=3)
    #test_data(x_linspace(0.0, 6.0, 13), modified_spherical_bessel, True, df=modified_spherical_bessel_deriv, n=4)
    #test_data(x_linspace(0.0, 6.0, 13), modified_spherical_bessel, True, df=modified_spherical_bessel_deriv, n=5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=0, rij=3.0, sigma=0.5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=1, rij=3.0, sigma=0.5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=2, rij=3.0, sigma=0.5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=3, rij=3.0, sigma=0.5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=4, rij=3.0, sigma=0.5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=5, rij=3.0, sigma=0.5)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=0, rij=1.75, sigma=0.8)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=1, rij=1.75, sigma=0.8)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=2, rij=1.75, sigma=0.8)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=3, rij=1.75, sigma=0.8)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=4, rij=1.75, sigma=0.8)
    #test_data(x_linspace(0.0, 6.0, 13), exponential_modified_spherical_bessel, False, n=5, rij=1.75, sigma=0.8)
    #test_data(x_geomspace(1.0E-6, 1.0, 11), exponential_modified_spherical_bessel, False, n=0, rij=0.9, sigma=0.3)
    #test_data(x_geomspace(1.0E-6, 1.0, 11), exponential_modified_spherical_bessel, False, n=1, rij=0.9, sigma=0.3)
    #test_data(x_geomspace(1.0E-6, 1.0, 11), exponential_modified_spherical_bessel, False, n=2, rij=0.9, sigma=0.3)
    #test_data(x_geomspace(1.0E-6, 1.0, 11), exponential_modified_spherical_bessel, False, n=3, rij=0.9, sigma=0.3)
    #test_data(x_geomspace(1.0E-6, 1.0, 11), exponential_modified_spherical_bessel, False, n=4, rij=0.9, sigma=0.3)
    #test_data(x_geomspace(1.0E-6, 1.0, 11), exponential_modified_spherical_bessel, False, n=5, rij=0.9, sigma=0.3)
