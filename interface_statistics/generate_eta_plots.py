#!/usr/bin/env python3
import glob
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


def _apply_plot_style():
    """Apply a clean publication-style theme inspired by PDF.ipynb."""
    plt.rcParams.update({
        "figure.figsize": (11, 6),
        "figure.dpi": 130,
        "savefig.dpi": 220,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 12,
        "font.family": "STIXGeneral",
        "mathtext.fontset": "stix",
        "axes.titleweight": "normal",
        "axes.labelweight": "normal",
        "axes.linewidth": 1.0,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "grid.linewidth": 0.8,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
    })


def _style_axes(ax):
    ax.grid(True, which="major")
    ax.minorticks_on()
    ax.grid(True, which="minor", alpha=0.10, linestyle=":", linewidth=0.6)


def _load_table(path):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            content = f.read()
        # Backward compatibility: some files were written with literal "\\n".
        if "\\n" in content:
            content = content.replace("\\n", "\n")
        data = np.loadtxt(io.StringIO(content), comments="#")
    except (OSError, ValueError):
        return np.empty((0, 0))
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return data


def _infer_tp(stats_dir):
    # Try to recover T0_ from the main simulation log in the case folder.
    case_dir = os.path.abspath(os.path.join(stats_dir, os.pardir))
    out_log = os.path.join(case_dir, "out.log")
    if not os.path.exists(out_log):
        return 1.0

    try:
        with open(out_log, "r", encoding="utf-8", errors="replace") as f:
            txt = f.read()
    except OSError:
        return 1.0

    m = re.search(r"T0\s*=\s*([0-9eE+\-.]+)", txt)
    if not m:
        return 1.0

    try:
        tp = float(m.group(1))
    except ValueError:
        return 1.0

    if tp <= 0.0:
        return 1.0
    return tp


def plot_window_stats(stats_file, out_dir, tp):
    if not os.path.exists(stats_file):
        print(f"[eta-plots] Missing file: {stats_file}")
        return False

    data = _load_table(stats_file)
    if data.size == 0 or data.shape[1] < 13:
        print(f"[eta-plots] No numeric rows yet in {stats_file}")
        return False
    # Columns (old, 13):
    # 0:t_end 1:i 2:t_rel 3:win_id 4:is_partial 5:nsnap 6:t_start 7:t_last
    # 8:eta_mean_win 9:eta_sigma_win 10:xi_skew 11:xi_kurt 12:weight_win
    # Extra columns (new, 21):
    # 13:sx_mean_win 14:sx_sigma_win 15:sx_skew 16:sx_kurt
    # 17:sz_mean_win 18:sz_sigma_win 19:sz_skew 20:sz_kurt
    t_rel = data[:, 2]
    t_rel0 = t_rel[0]
    t_dimless = (t_rel - t_rel0)/tp
    eta_mean = data[:, 8]
    eta_sigma = data[:, 9]
    xi_skew = data[:, 10]
    xi_kurt = data[:, 11]
    is_partial = data[:, 4]

    fig, axs = plt.subplots(3, 1, figsize=(11, 10), sharex=True, constrained_layout=True)

    c_mean = "#1f4e79"
    c_sig = "#d95f02"
    c_skew = "#1b9e77"
    c_kurt = "#6a3d9a"

    axs[0].plot(t_dimless, eta_mean, "-o", color=c_mean, label=r"$\langle \eta \rangle_{\rm win}$")
    axs[0].set_ylabel(r"$\langle \eta \rangle_{\rm win}$")
    axs[0].legend(loc="best", frameon=False)
    _style_axes(axs[0])

    axs[1].plot(t_dimless, eta_sigma, "-o", color=c_sig, label=r"$\sigma_{\eta,\rm win}$")
    axs[1].set_ylabel(r"$\sigma_{\eta,\rm win}$")
    axs[1].legend(loc="best", frameon=False)
    _style_axes(axs[1])

    axs[2].plot(t_dimless, xi_skew, "-o", color=c_skew, label=r"$\xi_{\rm skew}$")
    axs[2].plot(t_dimless, xi_kurt, "-o", color=c_kurt, label=r"$\xi_{\rm kurt}$")
    axs[2].axhline(0.0, color="0.3", lw=0.9, ls="--", alpha=0.7)
    axs[2].axhline(3.0, color="0.45", lw=0.9, ls=":", alpha=0.7)

    # Mark partial windows to make interrupted jobs obvious.
    for j, partial in enumerate(is_partial):
        if int(partial) == 1:
            axs[2].axvline(t_dimless[j], color="k", alpha=0.12, lw=0.9)
            axs[2].plot(t_dimless[j], xi_kurt[j], marker="v", color="k", ms=4, alpha=0.45)

    axs[2].set_xlabel(r"$\tau = (t_{\mathrm{rel}} - t_{\mathrm{rel},0})/T_p$")
    axs[2].set_ylabel(r"shape moments")
    axs[2].legend(loc="best", frameon=False, ncol=2)
    axs[2].xaxis.set_major_locator(MaxNLocator(nbins=8))
    _style_axes(axs[2])

    fig.suptitle(r"Windowed interface statistics", fontsize=20)
    out_png = os.path.join(out_dir, "eta_stats_windows.png")
    fig.savefig(out_png, dpi=180)
    plt.close(fig)
    print(f"[eta-plots] Wrote {out_png}")
    return True


def plot_slope_stats(stats_file, out_dir, tp):
    """Plot windowed slope statistics (∂η/∂x and ∂η/∂z) and cross-check γ₁(η) ~ κ(∂η/∂x)."""
    if not os.path.exists(stats_file):
        return False

    data = _load_table(stats_file)
    if data.size == 0 or data.shape[1] < 21:
        print(f"[eta-plots] Slope columns not yet available in {stats_file} (need 21 cols, got {data.shape[1] if data.size else 0})")
        return False

    t_rel = data[:, 2]
    t_dimless = (t_rel - t_rel[0]) / tp
    is_partial = data[:, 4]

    # Elevation moments (for cross-check)
    xi_skew = data[:, 10]   # skewness  of η
    xi_kurt = data[:, 11]   # kurtosis  of η

    # Slope moments: ∂η/∂x
    sx_mean  = data[:, 13]
    sx_sigma = data[:, 14]
    sx_skew  = data[:, 15]
    sx_kurt  = data[:, 16]

    # Slope moments: ∂η/∂z
    sz_mean  = data[:, 17]
    sz_sigma = data[:, 18]
    sz_skew  = data[:, 19]
    sz_kurt  = data[:, 20]

    # ------------------------------------------------------------------ #
    # Figure 1: Mean and std of slopes
    # ------------------------------------------------------------------ #
    fig1, axs1 = plt.subplots(2, 1, figsize=(11, 8), sharex=True, constrained_layout=True)

    axs1[0].plot(t_dimless, sx_mean, "-o", color="#1f4e79", label=r"$\langle s_x \rangle$")
    axs1[0].plot(t_dimless, sz_mean, "-s", color="#d95f02", label=r"$\langle s_z \rangle$")
    axs1[0].axhline(0.0, color="0.4", lw=0.9, ls="--", alpha=0.7)
    axs1[0].set_ylabel(r"mean slope")
    axs1[0].legend(loc="best", frameon=False, ncol=2)
    _style_axes(axs1[0])

    axs1[1].plot(t_dimless, sx_sigma, "-o", color="#1f4e79", label=r"$\sigma_{s_x}$")
    axs1[1].plot(t_dimless, sz_sigma, "-s", color="#d95f02", label=r"$\sigma_{s_z}$")
    axs1[1].set_xlabel(r"$\tau = (t_{\mathrm{rel}} - t_{\mathrm{rel},0})/T_p$")
    axs1[1].set_ylabel(r"slope std")
    axs1[1].legend(loc="best", frameon=False, ncol=2)
    axs1[1].xaxis.set_major_locator(MaxNLocator(nbins=8))
    _style_axes(axs1[1])

    for j, partial in enumerate(is_partial):
        if int(partial) == 1:
            for ax in axs1:
                ax.axvline(t_dimless[j], color="k", alpha=0.12, lw=0.9)

    fig1.suptitle(r"Slope statistics: mean and $\sigma$", fontsize=20)
    out1 = os.path.join(out_dir, "slope_stats_mean_sigma.png")
    fig1.savefig(out1, dpi=180)
    plt.close(fig1)
    print(f"[eta-plots] Wrote {out1}")

    # ------------------------------------------------------------------ #
    # Figure 2: Skewness of slopes
    # ------------------------------------------------------------------ #
    fig2, ax2 = plt.subplots(figsize=(11, 5), constrained_layout=True)
    ax2.plot(t_dimless, sx_skew, "-o", color="#1b9e77", label=r"$\gamma_1(s_x)$")
    ax2.plot(t_dimless, sz_skew, "-s", color="#e7298a", label=r"$\gamma_1(s_z)$")
    ax2.axhline(0.0, color="0.4", lw=0.9, ls="--", alpha=0.7)
    ax2.set_xlabel(r"$\tau = (t_{\mathrm{rel}} - t_{\mathrm{rel},0})/T_p$")
    ax2.set_ylabel(r"skewness")
    ax2.legend(loc="best", frameon=False, ncol=2)
    ax2.xaxis.set_major_locator(MaxNLocator(nbins=8))
    _style_axes(ax2)
    for j, partial in enumerate(is_partial):
        if int(partial) == 1:
            ax2.axvline(t_dimless[j], color="k", alpha=0.12, lw=0.9)
    fig2.suptitle(r"Skewness of surface slopes $\partial\eta/\partial x$ and $\partial\eta/\partial z$", fontsize=18)
    out2 = os.path.join(out_dir, "slope_stats_skewness.png")
    fig2.savefig(out2, dpi=180)
    plt.close(fig2)
    print(f"[eta-plots] Wrote {out2}")

    # ------------------------------------------------------------------ #
    # Figure 3: Kurtosis of slopes
    # ------------------------------------------------------------------ #
    fig3, ax3 = plt.subplots(figsize=(11, 5), constrained_layout=True)
    ax3.plot(t_dimless, sx_kurt, "-o", color="#6a3d9a", label=r"$\kappa(s_x)$")
    ax3.plot(t_dimless, sz_kurt, "-s", color="#ff7f00", label=r"$\kappa(s_z)$")
    ax3.axhline(3.0, color="0.4", lw=0.9, ls=":", alpha=0.8, label="Gaussian = 3")
    ax3.set_xlabel(r"$\tau = (t_{\mathrm{rel}} - t_{\mathrm{rel},0})/T_p$")
    ax3.set_ylabel(r"kurtosis")
    ax3.legend(loc="best", frameon=False, ncol=3)
    ax3.xaxis.set_major_locator(MaxNLocator(nbins=8))
    _style_axes(ax3)
    for j, partial in enumerate(is_partial):
        if int(partial) == 1:
            ax3.axvline(t_dimless[j], color="k", alpha=0.12, lw=0.9)
    fig3.suptitle(r"Kurtosis of surface slopes $\partial\eta/\partial x$ and $\partial\eta/\partial z$", fontsize=18)
    out3 = os.path.join(out_dir, "slope_stats_kurtosis.png")
    fig3.savefig(out3, dpi=180)
    plt.close(fig3)
    print(f"[eta-plots] Wrote {out3}")

    # ------------------------------------------------------------------ #
    # Figure 4: Cross-check γ₁(η) ~ κ(∂η/∂x)
    #   Left axis : γ₁(η) = xi_skew  and  γ₁(sx) = sx_skew
    #   Right axis: κ(η)  = xi_kurt  and  κ(sx)  = sx_kurt
    # ------------------------------------------------------------------ #
    fig4, ax4a = plt.subplots(figsize=(11, 5.5), constrained_layout=True)
    ax4b = ax4a.twinx()

    l1, = ax4a.plot(t_dimless, xi_skew,  "-o",  color="#1b9e77", lw=2.0, label=r"$\gamma_1(\eta)$")
    l2, = ax4a.plot(t_dimless, sx_skew,  "--^", color="#33a02c", lw=1.8, label=r"$\gamma_1(s_x)$")
    ax4a.axhline(0.0, color="0.5", lw=0.8, ls="--", alpha=0.6)
    ax4a.set_xlabel(r"$\tau = (t_{\mathrm{rel}} - t_{\mathrm{rel},0})/T_p$")
    ax4a.set_ylabel(r"skewness  $\gamma_1$", color="#1b9e77")
    ax4a.tick_params(axis="y", labelcolor="#1b9e77")

    l3, = ax4b.plot(t_dimless, xi_kurt,  "-s",  color="#6a3d9a", lw=2.0, label=r"$\kappa(\eta)$")
    l4, = ax4b.plot(t_dimless, sx_kurt,  "--v", color="#9e3d9a", lw=1.8, label=r"$\kappa(s_x)$")
    ax4b.axhline(3.0, color="0.55", lw=0.8, ls=":", alpha=0.7)
    ax4b.set_ylabel(r"kurtosis  $\kappa$", color="#6a3d9a")
    ax4b.tick_params(axis="y", labelcolor="#6a3d9a")

    ax4a.xaxis.set_major_locator(MaxNLocator(nbins=8))
    _style_axes(ax4a)

    handles = [l1, l2, l3, l4]
    ax4a.legend(handles=handles, loc="best", frameon=False, ncol=4)
    for j, partial in enumerate(is_partial):
        if int(partial) == 1:
            ax4a.axvline(t_dimless[j], color="k", alpha=0.12, lw=0.9)

    fig4.suptitle(r"Cross-check: $\gamma_1(\eta)$ vs $\kappa(\partial\eta/\partial x)$", fontsize=18)
    out4 = os.path.join(out_dir, "slope_crosscheck_skew_kurt.png")
    fig4.savefig(out4, dpi=180)
    plt.close(fig4)
    print(f"[eta-plots] Wrote {out4}")

    return True


def plot_pdf_files(stats_dir, out_dir):
    pdf_files = sorted(glob.glob(os.path.join(stats_dir, "eta_pdf_window_*.out")))
    if not pdf_files:
        print("[eta-plots] No eta_pdf_window_*.out files found")
        return False

    latest_png = None
    latest_eta_png = None
    for path in pdf_files:
        data = _load_table(path)
        if data.size == 0 or data.shape[1] < 6:
            print(f"[eta-plots] Skipping empty/incomplete file: {path}")
            continue
        # Columns: 0:eta_center 1:xi_center 2:pdf_xi 3:pdf_gauss 4:count 5:pdf_eta
        xi = data[:, 1]
        pdf_xi = data[:, 2]
        pdf_gauss = data[:, 3]
        eta = data[:, 0]
        pdf_eta = data[:, 5]

        base = os.path.basename(path).replace(".out", "")

        fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
        ax.plot(xi, pdf_xi, color="#0b6e4f", lw=2.2, label=r"$p(\xi)$")
        ax.plot(xi, pdf_gauss, "--", color="black", lw=1.6, label=r"Gaussian")
        ax.fill_between(xi, pdf_xi, 0.0, color="#0b6e4f", alpha=0.10)
        ax.set_xlabel(r"$\xi$")
        ax.set_ylabel(r"PDF")
        ax.set_title(base.replace("eta_pdf_window_", "Window "))
        ax.legend(loc="best", frameon=False)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=8))
        _style_axes(ax)
        out_png = os.path.join(out_dir, f"{base}_xi.png")
        fig.savefig(out_png, dpi=180)
        plt.close(fig)

        fig2, ax2 = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
        ax2.plot(eta, pdf_eta, color="#1f4e79", lw=2.2, label=r"$p(\eta)$")
        ax2.fill_between(eta, pdf_eta, 0.0, color="#1f4e79", alpha=0.10)
        ax2.set_xlabel(r"$\eta$")
        ax2.set_ylabel(r"PDF$(\eta)$")
        ax2.set_title(base.replace("eta_pdf_window_", "Window "))
        ax2.legend(loc="best", frameon=False)
        ax2.xaxis.set_major_locator(MaxNLocator(nbins=8))
        _style_axes(ax2)
        out_png2 = os.path.join(out_dir, f"{base}_eta.png")
        fig2.savefig(out_png2, dpi=180)
        plt.close(fig2)

        latest_png = out_png
        latest_eta_png = out_png2

    if latest_png:
        latest_link = os.path.join(out_dir, "eta_pdf_latest_xi.png")
        try:
            if os.path.islink(latest_link) or os.path.exists(latest_link):
                os.remove(latest_link)
            os.symlink(os.path.basename(latest_png), latest_link)
            print(f"[eta-plots] Updated {latest_link} -> {os.path.basename(latest_png)}")
        except OSError:
            pass

    if latest_eta_png:
        latest_link_eta = os.path.join(out_dir, "eta_pdf_latest_eta.png")
        try:
            if os.path.islink(latest_link_eta) or os.path.exists(latest_link_eta):
                os.remove(latest_link_eta)
            os.symlink(os.path.basename(latest_eta_png), latest_link_eta)
            print(f"[eta-plots] Updated {latest_link_eta} -> {os.path.basename(latest_eta_png)}")
        except OSError:
            pass
    else:
        print("[eta-plots] No valid PDF window files to plot")
        return False

    print(f"[eta-plots] Wrote plots for {len(pdf_files)} PDF windows")
    return True


def main():
    _apply_plot_style()

    stats_dir = sys.argv[1] if len(sys.argv) > 1 else "statistics"
    stats_file = os.path.join(stats_dir, "eta_stats_window.out")
    out_dir = os.path.join(stats_dir, "plots")
    os.makedirs(out_dir, exist_ok=True)

    tp = _infer_tp(stats_dir)
    print(f"[eta-plots] Using Tp = {tp:.10e}")

    ok1 = plot_window_stats(stats_file, out_dir, tp)
    ok2 = plot_pdf_files(stats_dir, out_dir)
    ok3 = plot_slope_stats(stats_file, out_dir, tp)

    if not (ok1 or ok2 or ok3):
        print("[eta-plots] Nothing to plot yet")
        return 0

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
