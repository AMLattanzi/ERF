#!/usr/bin/env python3
"""
Fit the temporal evolution of the volume-integrated acoustic energy written by
ERF (erf.acoustic_energy_data_log) to an exponential,

    E_a(t) = E_a(t0) * exp(sigma * (t - t0)),

by a least-squares linear fit of ln(E_a) against t. sigma < 0 means the
acoustic energy decays; since E_a is quadratic in the perturbation amplitude,
the amplitude itself evolves as exp(sigma/2 * t).

Log columns: time  Ea_pot  Ea_kin  Ea_tot

Usage:
    python fit_acoustic_energy.py [acoustic_energy.log] [--tmin T] [--tmax T]
                                  [--column tot|pot|kin] [--no-plot] [--save FILE]
"""
import argparse

import numpy as np

COLUMNS = {"pot": 1, "kin": 2, "tot": 3}


def read_log(filename):
    data = np.loadtxt(filename, comments="#", ndmin=2)
    # A restart appends to the log, so keep only the last entry for each time
    _, idx = np.unique(data[::-1, 0], return_index=True)
    return data[len(data) - 1 - idx]


def fit_log_energy(t, E):
    lnE = np.log(E)
    (sigma, intercept), cov = np.polyfit(t, lnE, 1, cov=True)
    resid = lnE - (sigma * t + intercept)
    r2 = 1.0 - np.sum(resid**2) / np.sum((lnE - lnE.mean())**2)
    return sigma, intercept, np.sqrt(cov[0, 0]), r2


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("logfile", nargs="?", default="acoustic_energy.log")
    parser.add_argument("--tmin", type=float, default=None, help="start of fit window [s]")
    parser.add_argument("--tmax", type=float, default=None, help="end of fit window [s]")
    parser.add_argument("--column", choices=COLUMNS.keys(), default="tot",
                        help="energy component to fit (default: tot)")
    parser.add_argument("--no-plot", dest="plot", action="store_false",
                        help="do not show the plot of ln(E) and the linear fit")
    parser.add_argument("--save", default=None, help="save the plot to this file")
    args = parser.parse_args()

    data = read_log(args.logfile)
    t = data[:, 0]
    E = data[:, COLUMNS[args.column]]

    mask = E > 0.0
    if args.tmin is not None:
        mask &= t >= args.tmin
    if args.tmax is not None:
        mask &= t <= args.tmax
    if np.count_nonzero(mask) < 3:
        raise SystemExit("Need at least 3 points with E > 0 in the fit window")

    sigma, intercept, sigma_err, r2 = fit_log_energy(t[mask], E[mask])

    print(f"File            : {args.logfile}")
    print(f"Component       : Ea_{args.column}")
    print(f"Fit window      : [{t[mask][0]:g}, {t[mask][-1]:g}] s  ({np.count_nonzero(mask)} points)")
    print(f"sigma           : {sigma:.6e} +/- {sigma_err:.2e} 1/s   (E ~ exp(sigma t))")
    print(f"amplitude rate  : {sigma / 2:.6e} 1/s   (|p'|,|u'| ~ exp(sigma t / 2))")
    if sigma != 0.0:
        print(f"e-folding time  : {1.0 / abs(sigma):.6e} s")
    print(f"E(t=0) from fit : {np.exp(intercept):.6e}")
    print(f"R^2 of ln(E) fit: {r2:.6f}")

    if args.plot or args.save:
        try:
            import matplotlib
        except ImportError:
            raise SystemExit("--plot/--save need matplotlib (pip install matplotlib)")
        if not args.plot:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7, 4.5))
        pos = E > 0.0
        ax.plot(t[pos], np.log(E[pos]), color="0.6", lw=1.5, label=f"ln(Ea_{args.column})")
        tf = t[mask]
        ax.plot(tf, np.log(E[mask]), "o", color="0.3", ms=3, label="points in fit window")
        ax.plot(tf, sigma * tf + intercept, color="C0", lw=2,
                label=f"fit: sigma = {sigma:.3e} 1/s, R$^2$ = {r2:.4f}")
        ax.axvspan(tf[0], tf[-1], color="C0", alpha=0.08)
        ax.set_xlabel("time [s]")
        ax.set_ylabel("ln(acoustic energy)")
        ax.legend(frameon=False)
        ax.grid(alpha=0.3)
        fig.tight_layout()
        if args.save:
            fig.savefig(args.save, dpi=150)
        if args.plot:
            plt.show()


if __name__ == "__main__":
    main()
