import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_file(fname):
    data   = np.loadtxt(fname, comments="#")
    psi, f, ffp, f_calc, f_err = data.T
    return psi, f, ffp, f_calc, f_err

cubic_file = "intg_ffp_cubic_test.dat"
gauss_file = "intg_ffp_gauss_test.dat"
exp_file   = "intg_ffp_exp_test.dat"
sin_file   = "intg_ffp_sin_test.dat"


def main():

    profiles = ["cubic", "gauss", "exp", "sin"]

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(11,8), constrained_layout=True)

    for ax, name in zip(axs.flat, profiles):
        fname = "intg_ffprime_"+name+"_test.dat"

        psi, f, ffp, f_calc, f_err = read_file(fname)
        xvals = (psi - psi[0]) / (psi[-1] - psi[0])
        xlabel = r"$\psi_N$"

        ax.plot(xvals, f, "k-", lw=2, label="F (exact)")
        ax.plot(xvals, f_calc, "r--", lw=1.5, label="F (reconstruct)")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("F")
        ax.grid(True)

        ax_err = ax.twinx()
        ax_err.plot(xvals, f_err, "b:", lw=1, label="error")
        ax_err.set_ylabel("F_calc - F", color="b")
        ax_err.tick_params(axis="y", colors="b")

        max_err = np.max(np.abs(f_err))
        ax.set_title(f"{name} (max |err| = {max_err:.2e})")

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax_err.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc="best")

    fig.savefig("intg_ffprime_test.png")
    plt.show()

if __name__ == "__main__":
    main()