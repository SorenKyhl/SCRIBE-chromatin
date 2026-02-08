"""
Knockout experiment module for SCRIBE.

Given a converged simulation directory, runs knockout simulations for each
mark (zeroing out that mark's row/column in the chi matrix), then analyzes
the results with publication-quality plots.
"""

import copy
import csv
import logging
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scribe import utils
from scribe.analysis import get_contactmap, get_diagonal, get_SCC, plot_diff
from scribe.scribe_sim import ScribeSim

logger = logging.getLogger(__name__)


class KnockoutExperiment:
    """Manages a full knockout experiment workflow.

    For each mark, creates a modified chi matrix with that mark's interactions
    zeroed out, runs a simulation, and analyzes the impact on Hi-C agreement.

    Args:
        sim_dir: Path to converged simulation directory (must contain config.json
            and sequence files).
        output_dir: Where to write knockout results.
        mark_names: Ordered list of mark names matching sequence order. If None,
            reads from config["mark_names"].
        gthic: Experimental Hi-C array, or path to a .npy file.
    """

    def __init__(self, sim_dir, output_dir, mark_names=None, gthic=None):
        self.sim_dir = Path(sim_dir).resolve()
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.config = utils.load_json(self.sim_dir / "config.json")
        self.seqs = utils.load_sequences_from_dir(self.sim_dir)

        # Resolve mark names
        if mark_names is not None:
            self.marks = list(mark_names)
        elif "mark_names" in self.config:
            self.marks = list(self.config["mark_names"])
        else:
            raise ValueError("mark_names must be provided or present in config['mark_names']")

        assert len(self.marks) == len(
            self.seqs
        ), f"mark_names length ({len(self.marks)}) != sequences length ({len(self.seqs)})"

        # Resolve experimental Hi-C
        if gthic is None:
            self.gthic = None
        elif isinstance(gthic, (str, Path)):
            self.gthic = np.load(gthic)
        else:
            self.gthic = np.asarray(gthic)

    def _zero_out_mark(self, config, mark_idx):
        """Zero out the row and column of the chi matrix for a given mark index."""
        config = copy.deepcopy(config)
        chis = np.array(config["chis"])
        chis[mark_idx, :] = 0.0
        chis[:, mark_idx] = 0.0
        config["chis"] = chis.tolist()
        return config

    def _wildtype_contacts_path(self):
        """Path to wildtype production contacts file."""
        return self.sim_dir / "production_out" / "contacts.txt"

    def run(self, equilib_sweeps=10_000, production_sweeps=50_000, parallel=7):
        """Run knockout simulations for each mark.

        For each mark, zeroes out that mark's row/column in the chi matrix and
        runs an equilibration + production simulation.

        Args:
            equilib_sweeps: Number of equilibration sweeps.
            production_sweeps: Number of production sweeps.
            parallel: Number of parallel production simulations per knockout.
        """
        for idx, mark in enumerate(self.marks):
            ko_dir = self.output_dir / f"ko_{mark}"
            contacts_file = ko_dir / "production_out" / "contacts.txt"

            if contacts_file.exists():
                logger.info(f"Skipping {mark} — already completed ({contacts_file})")
                continue

            # Clean up partial runs so ScribeSim can create the directory fresh
            if ko_dir.exists():
                logger.info(f"Removing incomplete directory: {ko_dir}")
                shutil.rmtree(ko_dir)

            logger.info(f"Running knockout: {mark} (index {idx})")

            ko_config = self._zero_out_mark(self.config, idx)

            # Normalize bead_type_files to match what ScribeSim.setup() writes
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

            logger.info(f"Done: {mark} knockout -> {ko_dir}/")

        logger.info("All knockouts complete.")

    def analyze(self):
        """Analyze knockout results and generate plots.

        Produces:
            - knockout_scc.png: bar chart of % SCC decrease from wildtype
            - diff_{mark}.png: difference heatmap for each knockout vs wildtype
            - ps_{mark}.png: P(s) comparison for each knockout vs wildtype
            - chi_matrix.png: annotated chi matrix heatmap
            - knockout_results.csv: summary table

        Returns:
            dict: {mark: {"scc": float, "pct_decrease": float}}
        """
        if self.gthic is None:
            raise ValueError("gthic is required for analysis")

        # Wildtype contact map
        wt_contacts = get_contactmap(str(self._wildtype_contacts_path()))
        baseline_scc = get_SCC(wt_contacts, self.gthic)
        logger.info(f"Baseline (wildtype) SCC: {baseline_scc:.4f}")

        # Compute knockout SCCs
        results = {}
        for mark in self.marks:
            contacts_path = self.output_dir / f"ko_{mark}" / "production_out" / "contacts.txt"
            ko_contacts = get_contactmap(str(contacts_path))
            ko_scc = get_SCC(ko_contacts, self.gthic)
            pct_decrease = (baseline_scc - ko_scc) / baseline_scc * 100
            results[mark] = {"scc": ko_scc, "pct_decrease": pct_decrease}
            logger.info(f"  {mark:12s}  SCC={ko_scc:.4f}  decrease={pct_decrease:.2f}%")

        # 1. SCC bar chart
        self._plot_scc_bar(results, baseline_scc)

        # 2. Diff heatmaps
        self._plot_diff_heatmaps(wt_contacts)

        # 3. P(s) comparisons
        self._plot_ps_comparisons(wt_contacts)

        # 4. Chi matrix heatmap
        self._plot_chi_matrix()

        # 5. Results CSV
        self._save_results_csv(results, baseline_scc)

        return results

    def _plot_scc_bar(self, results, baseline_scc):
        """Plot SCC % decrease bar chart, sorted by impact."""
        sorted_marks = sorted(results, key=lambda m: results[m]["pct_decrease"], reverse=True)
        pct_values = [results[m]["pct_decrease"] for m in sorted_marks]

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.bar(
            range(len(sorted_marks)), pct_values, color="#4878CF", edgecolor="black", linewidth=0.5
        )
        ax.set_xticks(range(len(sorted_marks)))
        ax.set_xticklabels(sorted_marks, rotation=45, ha="right", fontsize=10)
        ax.set_ylabel("SCC % decrease from wildtype")
        ax.set_title(f"Knockout impact on Hi-C agreement (baseline SCC={baseline_scc:.4f})")
        ax.axhline(0, color="black", linewidth=0.5)
        fig.tight_layout()
        outpath = self.output_dir / "knockout_scc.png"
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        logger.info(f"Saved {outpath}")

    def _plot_diff_heatmaps(self, wt_contacts):
        """Plot knockout − wildtype difference heatmaps."""
        for mark in self.marks:
            contacts_path = self.output_dir / f"ko_{mark}" / "production_out" / "contacts.txt"
            ko_contacts = get_contactmap(str(contacts_path))
            plot_diff(ko_contacts, wt_contacts, title=f"ko_{mark} − wildtype")
            outpath = self.output_dir / f"diff_{mark}.png"
            plt.savefig(outpath, dpi=150)
            plt.close()
            logger.info(f"Saved {outpath}")

    def _plot_ps_comparisons(self, wt_contacts):
        """Plot P(s) curves: each knockout vs wildtype."""
        wt_diag = get_diagonal(wt_contacts)
        x_vec = np.linspace(1 / len(wt_diag), 1, len(wt_diag))

        for mark in self.marks:
            contacts_path = self.output_dir / f"ko_{mark}" / "production_out" / "contacts.txt"
            ko_contacts = get_contactmap(str(contacts_path))
            ko_diag = get_diagonal(ko_contacts)

            fig, ax = plt.subplots(figsize=(6, 4))
            ax.semilogy(x_vec, wt_diag, "k-", label="wildtype", linewidth=1.5)
            ax.semilogy(x_vec, ko_diag, "r-", label=f"ko_{mark}", linewidth=1.5)
            ax.set_xlabel("Genomic distance (fraction of region)")
            ax.set_ylabel("Contact probability")
            ax.set_title(f"P(s): ko_{mark} vs wildtype")
            ax.legend()
            lower_bound = max(min(wt_diag.min(), ko_diag.min()) / 5, 1e-7)
            ax.set_ylim(lower_bound, 1)
            fig.tight_layout()
            outpath = self.output_dir / f"ps_{mark}.png"
            fig.savefig(outpath, dpi=150)
            plt.close(fig)
            logger.info(f"Saved {outpath}")

    def _plot_chi_matrix(self):
        """Plot annotated chi matrix heatmap."""
        chis = np.array(self.config["chis"])
        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(chis, cmap="bwr", vmin=-np.max(np.abs(chis)), vmax=np.max(np.abs(chis)))
        ax.set_xticks(range(len(self.marks)))
        ax.set_xticklabels(self.marks, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(self.marks)))
        ax.set_yticklabels(self.marks, fontsize=9)
        ax.set_title("Chi interaction matrix")

        # Annotate cells with values
        for i in range(len(self.marks)):
            for j in range(len(self.marks)):
                ax.text(j, i, f"{chis[i, j]:.2f}", ha="center", va="center", fontsize=7)

        fig.colorbar(im, ax=ax)
        fig.tight_layout()
        outpath = self.output_dir / "chi_matrix.png"
        fig.savefig(outpath, dpi=150)
        plt.close(fig)
        logger.info(f"Saved {outpath}")

    def _save_results_csv(self, results, baseline_scc):
        """Save results as CSV and print to console."""
        outpath = self.output_dir / "knockout_results.csv"
        sorted_marks = sorted(results, key=lambda m: results[m]["pct_decrease"], reverse=True)

        print(f"\n{'Mark':>14s}  {'SCC':>8s}  {'% decrease':>10s}")
        print(f"{'─'*14}  {'─'*8}  {'─'*10}")
        print(f"{'wildtype':>14s}  {baseline_scc:8.4f}  {'—':>10s}")

        with open(outpath, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["mark", "scc", "pct_decrease"])
            writer.writerow(["wildtype", f"{baseline_scc:.4f}", ""])
            for mark in sorted_marks:
                r = results[mark]
                print(f"{mark:>14s}  {r['scc']:8.4f}  {r['pct_decrease']:10.2f}%")
                writer.writerow([mark, f"{r['scc']:.4f}", f"{r['pct_decrease']:.2f}"])

        logger.info(f"Saved {outpath}")

    def verify_mark_order(self):
        """Verify mark ordering by correlating sequences with ChIP-seq data.

        Correlates each sequence file from sim_dir with named ChIP-seq loaded via
        DataPipeline (all available marks). Asserts the best-matching mark for each
        index matches self.marks.

        Useful for pre-fix simulations where ordering was glob-dependent.
        """
        from scribe.data_pipeline import DataPipeline

        pipeline = DataPipeline(cell="HCT116_auxin", chrom=2, nbeads=len(self.seqs[0]))
        chipseq = pipeline.load_chipseq()

        verified = []
        for i, seq in enumerate(self.seqs):
            best_mark = None
            best_corr = -np.inf
            for mark_name, chip_signal in chipseq.items():
                corr = np.corrcoef(seq, chip_signal)[0, 1]
                if corr > best_corr:
                    best_corr = corr
                    best_mark = mark_name
            verified.append(best_mark)
            logger.info(f"  seq{i+1} -> {best_mark} (r={best_corr:.4f})")

        assert (
            verified == self.marks
        ), f"Mark order mismatch!\n  Expected: {self.marks}\n  Got:      {verified}"
        logger.info("Mark order verified.")
