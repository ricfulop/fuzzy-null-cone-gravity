"""Regenerate the spectrum on the archived Schwarzschild grid.
Units G = c = hbar = M = 1. Boundary conventions follow
orig3_sec6_semiclassical_einstein_and_continuum.py.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import simpson
from scipy.special import lambertw

DATA = Path(__file__).resolve().parent / 'review_data'


def transmission(j, frequency, samples=2400):
    # Invert the tortoise coordinate r* = r + 2 log(r/2 - 1).
    x = np.linspace(-40.0, 160.0, samples)
    r = 2 * (1 + lambertw(np.exp(x / 2 - 1)).real)
    f = 1 - 2 / r
    L = j * (j + 1)
    potential = f * (L / r**2 - 6 / r**3)
    h = x[1] - x[0]
    w = frequency
    c1 = 1j * L / (2 * w)
    c2 = -L * (L - 2) / (8 * w**2)

    def outgoing(index):
        return np.exp(1j * w * x[index]) * (1 + c1/r[index] + c2/r[index]**2)

    # Numerov recurrence from infinity toward the horizon for all frequencies.
    q = w[None, :]**2 - potential[:, None]
    a = 1 + (h**2 / 12) * q
    b = 2 - (5 * h**2 / 6) * q
    ahead = outgoing(-1)
    current = outgoing(-2)
    for i in range(samples - 2, 0, -1):
        previous = (b[i]*current - a[i+1]*ahead) / a[i-1]
        ahead, current = current, previous
    p0, p1 = np.exp(1j*w*x[0]), np.exp(1j*w*x[1])
    m0, m1 = p0.conj(), p1.conj()
    incident = (current*m1 - ahead*m0) / (p0*m1 - m0*p1)
    return 1 / np.abs(incident)**2


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    w = np.linspace(0.06, 1.60, 155)
    partial = []
    for j in range(2, 7):
        t = transmission(j, w)
        partial.append((2*j+1)/np.pi * w*t / np.expm1(8*np.pi*w))
    curves = np.array(partial).T
    total = curves.sum(axis=1)
    alpha = simpson(total, x=w)
    assert abs(alpha/3.836658052716681e-6 - 1) < 1e-7, alpha
    print('Integrated coefficient', alpha)
    print('Quadrupole fraction', simpson(curves[:, 0], x=w) / alpha)
    np.savetxt(DATA/'spectrum.csv', np.column_stack((w, curves, total)), delimiter=',',
               header='M_omega,j2,j3,j4,j5,j6,total', comments='', fmt='%.14g')


if __name__ == '__main__':
    main()
