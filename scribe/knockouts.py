"""
Cross-term knockout experiments for SCRIBE.

Knocks out all off-diagonal interactions in the chi matrix, leaving only the
main diagonal (self-interactions). This tests the importance of cross-mark
interactions for reproducing 3D chromatin structure.
"""

import copy
import csv
import logging
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scribe import utils
from scribe.analysis import (
    get_contactmap,
    get_diagonal,
    get_oe,
    get_pearson,
    get_RMSE,
    get_SCC,
    plot_contactmap,
    plot_diff,
    plot_scatter,
    plot_tri,
)
from scribe.scribe_sim import ScribeSim

logger = logging.getLogger(__name__)


def _load_gthic(gthic):
    """Resolve gthic to a numpy array."""
    if gthic is None:
        return None
    if isinstance(gthic, (str, Path)):
        return np.load(gthic)
    return np.asarray(gthic)


def _save_fig(outpath, dpi=150):
    """Save current figure and close."""
    plt.savefig(outpath, dpi=dpi)
    plt.close()
    logger.info(f"Saved {outpath}")


class CrossTermKnockout:
    """Knockout of all cross-term (off-diagonal) chi interactions.

    Zeroes out all off-diagonal elements of the chi matrix, keeping only
    the main diagonal (self-interaction terms). Runs a single simulation
    and compares to the wildtype.

    Args:
        sim_dir: Path to converged simulation directory (e.g. iteration9).
        output_dir: Where to write knockout results.
        mark_names: Ordered list of mark names matching sequence order.
        gthic: Experimental Hi-C array, or path to a .npy file.
    """

    def __init__(self, sim_dir, output_dir, mark_names, gthic=None):
        self.sim_dir = Path(sim_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.config = utils.load_json(self.sim_dir / "config.json")
        self.seqs = utils.load_sequences_from_dir(self.sim_dir)
        self.marks = list(mark_names)

        assert len(self.marks) == len(
            self.seqs
        ), f"mark_names length ({len(self.marks)}) != sequences length ({len(self.seqs)})"

        self.gthic = _load_gthic(gthic)

    @staticmethod
    def _zero_cross_terms(config):
        """Return a config copy with all off-diagonal chi elements zeroed."""
        config = copy.deepcopy(config)
        chis = np.array(config["chis"])
        config["chis"] = np.diag(np.diag(chis)).tolist()
        return config

    def _wildtype_contacts_path(self):
        return self.sim_dir / "production_out" / "contacts.txt"

    def _ko_contacts_path(self):
        return self.output_dir / "ko_cross_terms" / "production_out" / "contacts.txt"

    def run(self, equilib_sweeps=10_000, production_sweeps=50_000, parallel=7):
        """Run the diagonal-only knockout simulation.

        Args:
            equilib_sweeps: Number of equilibration sweeps.
            production_sweeps: Number of production sweeps.
            parallel: Number of parallel production simulations.
        """
        ko_dir = self.output_dir / "ko_cross_terms"
        contacts_file = ko_dir / "production_out" / "contacts.txt"

        if contacts_file.exists():
            logger.info(f"Skipping — already completed ({contacts_file})")
            return

        if ko_dir.exists():
            logger.info(f"Removing incomplete directory: {ko_dir}")
            shutil.rmtree(ko_dir)

        logger.info("Running cross-term knockout (diagonal-only chi)")

        ko_config = self._zero_cross_terms(self.config)
        ko_config["bead_type_files"] = [f"seq{i+1}.txt" for i in range(len(self.seqs))]

        sim = ScribeSim(
            root=ko_dir,
            config=ko_config,
            seqs=self.seqs,
        )
        sim.run_eq(
            equilibrium_sweeps=equilib_sweeps,
            production_sweeps=production_sweeps,
            parallel_simulations=parallel,
        )

        logger.info(f"Done: cross-term knockout -> {ko_dir}/")

    def analyze(self):
        """Analyze knockout results with comprehensive plots.

        Produces:
            - chi_comparison.png: side-by-side original vs diagonal-only chi
            - diff_cross_terms.png: difference heatmap (ko − wildtype)
            - tri_linear.png: split-triangle contact map (linear scale)
            - tri_log.png: split-triangle contact map (log scale)
            - contactmap_ko_linear.png: knockout contact map (linear)
            - contactmap_ko_log.png: knockout contact map (log)
            - contactmap_wt_linear.png: wildtype contact map (linear)
            - contactmap_wt_log.png: wildtype contact map (log)
            - ps_cross_terms.png: P(s) comparison
            - scatter_ko_vs_wt.png: scatter plot ko vs wildtype
            - scatter_ko_vs_exp.png: scatter plot ko vs experimental
            - oe_comparison.png: observed/expected comparison
            - knockout_results.csv: summary metrics

        Returns:
            dict: {"cross_terms": {"scc_vs_exp": float, "pct_decrease": float, ...}}
        """
        if self.gthic is None:
            raise ValueError("gthic is required for analysis")

        wt_contacts = get_contactmap(str(self._wildtype_contacts_path()))
        ko_contacts = get_contactmap(str(self._ko_contacts_path()))

        # Compute metrics
        wt_scc = get_SCC(wt_contacts, self.gthic)
        ko_scc = get_SCC(ko_contacts, self.gthic)
        ko_vs_wt_scc = get_SCC(ko_contacts, wt_contacts)
        ko_pearson = get_pearson(ko_contacts, self.gthic)
        wt_pearson = get_pearson(wt_contacts, self.gthic)
        ko_rmse = get_RMSE(ko_contacts, self.gthic)
        wt_rmse = get_RMSE(wt_contacts, self.gthic)
        pct_decrease = (wt_scc - ko_scc) / wt_scc * 100

        logger.info(f"Wildtype   SCC={wt_scc:.4f}  Pearson={wt_pearson:.4f}  RMSE={wt_rmse:.4f}")
        logger.info(f"Knockout   SCC={ko_scc:.4f}  Pearson={ko_pearson:.4f}  RMSE={ko_rmse:.4f}")
        logger.info(f"KO vs WT   SCC={ko_vs_wt_scc:.4f}")
        logger.info(f"SCC decrease: {pct_decrease:.2f}%")

        results = {
            "cross_terms": {
                "scc_vs_exp": ko_scc,
                "pct_decrease": pct_decrease,
                "pearson_vs_exp": ko_pearson,
                "rmse_vs_exp": ko_rmse,
                "scc_vs_wt": ko_vs_wt_scc,
            }
        }

        # 1. Chi matrix comparison
        self._plot_chi_comparison()

        # 2. Difference heatmap (ko − wt)
        plot_diff(ko_contacts, wt_contacts, title="ko_cross_terms − wildtype")
        _save_fig(self.output_dir / "diff_cross_terms.png")

        # 3. Split-triangle: ko (upper) vs wt (lower), linear and log
        plot_tri(ko_contacts, wt_contacts, title="upper: knockout | lower: wildtype")
        _save_fig(self.output_dir / "tri_linear.png")

        plot_tri(
            ko_contacts, wt_contacts, title="upper: knockout | lower: wildtype (log)", log=True
        )
        _save_fig(self.output_dir / "tri_log.png")

        # 4. Contact maps: knockout (linear and log)
        plot_contactmap(ko_contacts, title="Cross-term knockout")
        _save_fig(self.output_dir / "contactmap_ko_linear.png")

        plot_contactmap(ko_contacts, title="Cross-term knockout (log)", log=True)
        _save_fig(self.output_dir / "contactmap_ko_log.png")

        # 5. Contact maps: wildtype (linear and log)
        plot_contactmap(wt_contacts, title="Wildtype")
        _save_fig(self.output_dir / "contactmap_wt_linear.png")

        plot_contactmap(wt_contacts, title="Wildtype (log)", log=True)
        _save_fig(self.output_dir / "contactmap_wt_log.png")

        # 6. P(s) comparison
        self._plot_ps(wt_contacts, ko_contacts)

        # 7. Scatter plots
        plot_scatter(ko_contacts, wt_contacts, label1="knockout", label2="wildtype")
        _save_fig(self.output_dir / "scatter_ko_vs_wt.png")

        plot_scatter(ko_contacts, self.gthic, label1="knockout", label2="experimental")
        _save_fig(self.output_dir / "scatter_ko_vs_exp.png")

        # 8. O/E comparison
        self._plot_oe_comparison(wt_contacts, ko_contacts)

        # 9. Results CSV
        self._save_results_csv(results, wt_scc, wt_pearson, wt_rmse)

        return results

    def _plot_ps(self, wt_contacts, ko_contacts):
        """Plot P(s) curves: knockout, wildtype, and experimental."""
        wt_diag = get_diagonal(wt_contacts)
        ko_diag = get_diagonal(ko_contacts)
        exp_diag = get_diagonal(self.gthic)
        x_vec = np.linspace(1 / len(wt_diag), 1, len(wt_diag))

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.semilogy(x_vec, exp_diag, "k--", label="experimental", linewidth=1.5)
        ax.semilogy(x_vec, wt_diag, "b-", label="wildtype", linewidth=1.5)
        ax.semilogy(x_vec, ko_diag, "r-", label="ko_cross_terms", linewidth=1.5)
        ax.set_xlabel("Genomic distance (fraction of region)")
        ax.set_ylabel("Contact probability")
        ax.set_title("P(s): cross-term knockout vs wildtype vs experimental")
        ax.legend()
        lower_bound = max(min(wt_diag.min(), ko_diag.min(), exp_diag.min()) / 5, 1e-7)
        ax.set_ylim(lower_bound, 1)
        fig.tight_layout()
        _save_fig(self.output_dir / "ps_cross_terms.png")

    def _plot_chi_comparison(self):
        """Plot side-by-side: original chi vs diagonal-only chi."""
        chis = np.array(self.config["chis"])
        diag_chis = np.diag(np.diag(chis))
        vmax = np.max(np.abs(chis))

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        for ax, mat, title in zip(
            axes,
            [chis, diag_chis],
            ["Original chi matrix", "Diagonal-only chi matrix"],
        ):
            im = ax.imshow(mat, cmap="bwr", vmin=-vmax, vmax=vmax)
            ax.set_xticks(range(len(self.marks)))
            ax.set_xticklabels(self.marks, rotation=45, ha="right", fontsize=9)
            ax.set_yticks(range(len(self.marks)))
            ax.set_yticklabels(self.marks, fontsize=9)
            ax.set_title(title)
            for i in range(len(self.marks)):
                for j in range(len(self.marks)):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=7)
            fig.colorbar(im, ax=ax)

        fig.tight_layout()
        _save_fig(self.output_dir / "chi_comparison.png")

    def _plot_oe_comparison(self, wt_contacts, ko_contacts):
        """Plot O/E matrices side-by-side: wildtype vs knockout."""
        wt_oe = get_oe(wt_contacts)
        ko_oe = get_oe(ko_contacts)

        fig, axes = plt.subplots(1, 2, figsize=(20, 9))
        for ax, mat, title in zip(
            axes,
            [wt_oe, ko_oe],
            ["Wildtype O/E", "Cross-term knockout O/E"],
        ):
            im = ax.imshow(mat, vmin=0, vmax=2, cmap="bwr")
            ax.set_title(title)
            fig.colorbar(im, ax=ax)

        fig.tight_layout()
        _save_fig(self.output_dir / "oe_comparison.png")

    def _save_results_csv(self, results, wt_scc, wt_pearson, wt_rmse):
        """Save results as CSV and print summary."""
        outpath = self.output_dir / "knockout_results.csv"
        r = results["cross_terms"]

        header = (
            f"\n{'Condition':>14s}  {'SCC':>8s}  {'Pearson':>8s}  {'RMSE':>8s}  {'% SCC dec':>10s}"
        )
        sep = f"{'─' * 14}  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 10}"
        print(header)
        print(sep)
        print(f"{'wildtype':>14s}  {wt_scc:8.4f}  {wt_pearson:8.4f}  {wt_rmse:8.4f}  {'—':>10s}")
        print(
            f"{'cross_terms':>14s}  {r['scc_vs_exp']:8.4f}  {r['pearson_vs_exp']:8.4f}  {r['rmse_vs_exp']:8.4f}  {r['pct_decrease']:10.2f}%"
        )

        with open(outpath, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                ["condition", "scc", "pearson", "rmse", "pct_scc_decrease", "scc_vs_wt"]
            )
            writer.writerow(
                ["wildtype", f"{wt_scc:.4f}", f"{wt_pearson:.4f}", f"{wt_rmse:.4f}", "", ""]
            )
            writer.writerow(
                [
                    "cross_terms",
                    f"{r['scc_vs_exp']:.4f}",
                    f"{r['pearson_vs_exp']:.4f}",
                    f"{r['rmse_vs_exp']:.4f}",
                    f"{r['pct_decrease']:.2f}",
                    f"{r['scc_vs_wt']:.4f}",
                ]
            )

        logger.info(f"Saved {outpath}")
