"""
Module 4 -- Comparative Reporter / Visualisation Suite
Generates publication-quality charts for pathogenicity profiles:
  - Category bar charts (VF / AMR distribution)
  - Comparative heatmaps (sample vs. sample)
  - TPM fold-change plots
  - Risk-level summary dashboards
Uses Plotly for interactive charts and Matplotlib/Seaborn for static exports.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from config import settings
from nanopore_pipeline.classifier.pathogenicity import CategoryProfile
from nanopore_pipeline.utils.logging_config import setup_logging

logger = setup_logging("reporter.visualisation")

RISK_COLOURS = {
    "Critical": "#d62728",
    "High": "#ff7f0e",
    "Medium": "#ffbb33",
    "Low": "#2ca02c",
    "None": "#cccccc",
}


class PathogenicityReporter:
    """Generates visual reports from classification profiles."""

    def __init__(self, output_dir: Path = settings.RESULTS_DIR):
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info("PathogenicityReporter output_dir=%s", output_dir)

    # ── Single-sample bar chart ───────────────────────────────────────────

    def plot_category_bar(
        self,
        profiles: list[CategoryProfile],
        sample_name: str,
        save_html: bool = True,
    ) -> go.Figure:
        categories = [p.category for p in profiles]
        hit_counts = [p.hit_count for p in profiles]
        colours = [RISK_COLOURS.get(p.risk_level, "#999") for p in profiles]
        tpm_values = [p.tpm for p in profiles]

        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=("Hit Count by Category", "TPM by Category"),
            shared_yaxes=True,
        )

        fig.add_trace(
            go.Bar(
                y=categories, x=hit_counts, orientation="h",
                marker_color=colours, name="Hits",
                text=[f"{c} ({r})" for c, r in zip(hit_counts, [p.risk_level for p in profiles])],
                textposition="auto",
            ),
            row=1, col=1,
        )

        fig.add_trace(
            go.Bar(
                y=categories, x=tpm_values, orientation="h",
                marker_color=colours, name="TPM",
                text=[f"{t:.1f}" for t in tpm_values],
                textposition="auto",
            ),
            row=1, col=2,
        )

        fig.update_layout(
            title=f"Pathogenicity Profile: {sample_name}",
            height=max(400, len(categories) * 40 + 200),
            showlegend=False,
        )

        if save_html:
            out = self.output_dir / f"{sample_name}_profile.html"
            fig.write_html(str(out))
            logger.info("Saved bar chart -> %s", out)

        return fig

    # ── Comparative heatmap ───────────────────────────────────────────────

    def plot_comparison_heatmap(
        self,
        comparison: dict,
        sample_a_name: str,
        sample_b_name: str,
        save_html: bool = True,
    ) -> go.Figure:
        categories = list(comparison.keys())
        tpm_a = [comparison[c]["sample_a_tpm"] for c in categories]
        tpm_b = [comparison[c]["sample_b_tpm"] for c in categories]
        fold_changes = [comparison[c]["fold_change_a_vs_b"] for c in categories]

        fig = go.Figure()

        fig.add_trace(go.Heatmap(
            z=[tpm_a, tpm_b],
            x=categories,
            y=[sample_a_name, sample_b_name],
            colorscale="YlOrRd",
            text=[[f"{v:.1f}" for v in tpm_a], [f"{v:.1f}" for v in tpm_b]],
            texttemplate="%{text}",
            colorbar_title="TPM",
        ))

        fig.update_layout(
            title=f"Comparative Pathogenicity: {sample_a_name} vs {sample_b_name}",
            xaxis_title="Pathogenicity Category",
            height=400,
        )

        if save_html:
            out = self.output_dir / f"compare_{sample_a_name}_vs_{sample_b_name}_heatmap.html"
            fig.write_html(str(out))
            logger.info("Saved heatmap -> %s", out)

        return fig

    # ── Fold-change bar chart ─────────────────────────────────────────────

    def plot_fold_change(
        self,
        comparison: dict,
        sample_a_name: str,
        sample_b_name: str,
        save_html: bool = True,
    ) -> go.Figure:
        categories = list(comparison.keys())
        fold_changes = [comparison[c]["fold_change_a_vs_b"] for c in categories]
        colours = ["#d62728" if fc > 2 else "#2ca02c" if fc < 0.5 else "#1f77b4"
                    for fc in fold_changes]

        fig = go.Figure(go.Bar(
            x=categories,
            y=fold_changes,
            marker_color=colours,
            text=[f"{fc:.1f}x" for fc in fold_changes],
            textposition="auto",
        ))

        fig.add_hline(y=1.0, line_dash="dash", line_color="gray",
                       annotation_text="No change")

        fig.update_layout(
            title=f"TPM Fold Change: {sample_a_name} / {sample_b_name}",
            xaxis_title="Category",
            yaxis_title="Fold Change",
            yaxis_type="log",
            height=500,
        )

        if save_html:
            out = self.output_dir / f"fold_change_{sample_a_name}_vs_{sample_b_name}.html"
            fig.write_html(str(out))
            logger.info("Saved fold-change chart -> %s", out)

        return fig

    # ── Risk summary dashboard ────────────────────────────────────────────

    def plot_risk_dashboard(
        self,
        profiles: list[CategoryProfile],
        sample_name: str,
        save_html: bool = True,
    ) -> go.Figure:
        risk_counts = {"Critical": 0, "High": 0, "Medium": 0, "Low": 0}
        for p in profiles:
            risk_counts[p.risk_level] = risk_counts.get(p.risk_level, 0) + 1

        fig = make_subplots(
            rows=1, cols=2,
            specs=[[{"type": "pie"}, {"type": "table"}]],
            subplot_titles=("Risk Distribution", "Category Details"),
        )

        # Pie chart
        labels = list(risk_counts.keys())
        values = list(risk_counts.values())
        pie_colours = [RISK_COLOURS[l] for l in labels]

        fig.add_trace(
            go.Pie(labels=labels, values=values, marker_colors=pie_colours,
                    hole=0.4, textinfo="label+value"),
            row=1, col=1,
        )

        # Table
        fig.add_trace(
            go.Table(
                header=dict(values=["Category", "Genes", "Hits", "TPM", "Risk"]),
                cells=dict(values=[
                    [p.category for p in profiles],
                    [p.unique_genes for p in profiles],
                    [p.hit_count for p in profiles],
                    [f"{p.tpm:.1f}" for p in profiles],
                    [p.risk_level for p in profiles],
                ]),
            ),
            row=1, col=2,
        )

        fig.update_layout(
            title=f"Risk Dashboard: {sample_name}",
            height=500,
        )

        if save_html:
            out = self.output_dir / f"{sample_name}_risk_dashboard.html"
            fig.write_html(str(out))
            logger.info("Saved risk dashboard -> %s", out)

        return fig

    # ── Alignment quality distribution plots ─────────────────────────────

    def plot_alignment_quality(
        self,
        hits_data: list[dict],
        sample_name: str,
        hit_type: str = "VF",
        pident_cutoff: float = 70.0,
        bitscore_cutoff: float = 50.0,
        save_html: bool = True,
    ) -> go.Figure:
        """
        Plot histograms of percent identity and bitscore distributions with
        cutoff lines so the selection thresholds can be justified in the thesis.

        Supervisor note (meeting ~17:00-21:33): "it would be nice to plot it
        somehow or to analyze it to devise a logic that you'd say this alignment
        truly speaks to this virulence gene ... plot the identity, the histogram
        of identity ... you can plot, for example, the distribution and this line
        on the 90 and say, here, see, I'm cutting off this tail and this is real
        stuff."

        Args:
            hits_data: list of dicts with keys 'pident' and 'bitscore'
                       (pass BlastHitRow objects as dicts via __dict__,
                        or pass raw dicts directly).
            sample_name:     Sample identifier used in title and filename.
            hit_type:        "VF" or "AMR" (for plot title labelling).
            pident_cutoff:   Vertical line on the identity histogram (default 70%).
            bitscore_cutoff: Vertical line on the bitscore histogram (default 50).
            save_html:       Write output to results directory.
        """
        if not hits_data:
            logger.warning("plot_alignment_quality: no hits supplied for %s", sample_name)
            return go.Figure()

        pidents = [h["pident"] for h in hits_data]
        bitscores = [h["bitscore"] for h in hits_data]

        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=(
                f"% Identity Distribution (cutoff={pident_cutoff}%)",
                f"Bit Score Distribution (cutoff={bitscore_cutoff})",
            ),
        )

        # Identity histogram
        fig.add_trace(
            go.Histogram(
                x=pidents, nbinsx=40,
                marker_color="#1f77b4", name="% Identity",
                hovertemplate="Identity: %{x:.1f}%<br>Count: %{y}<extra></extra>",
            ),
            row=1, col=1,
        )
        fig.add_vline(
            x=pident_cutoff, line_dash="dash", line_color="red",
            annotation_text=f"cutoff {pident_cutoff}%",
            annotation_position="top right",
            row=1, col=1,
        )

        # Bitscore histogram
        fig.add_trace(
            go.Histogram(
                x=bitscores, nbinsx=40,
                marker_color="#2ca02c", name="Bit Score",
                hovertemplate="Bitscore: %{x:.1f}<br>Count: %{y}<extra></extra>",
            ),
            row=1, col=2,
        )
        fig.add_vline(
            x=bitscore_cutoff, line_dash="dash", line_color="red",
            annotation_text=f"cutoff {bitscore_cutoff}",
            annotation_position="top right",
            row=1, col=2,
        )

        fig.update_layout(
            title=(
                f"Alignment Quality Distribution: {sample_name} ({hit_type})<br>"
                f"<sub>n={len(hits_data)} hits — use these distributions to justify your cutoff thresholds in the thesis</sub>"
            ),
            showlegend=False,
            height=450,
        )
        fig.update_xaxes(title_text="% Identity", row=1, col=1)
        fig.update_xaxes(title_text="Bit Score", row=1, col=2)
        fig.update_yaxes(title_text="Number of hits", row=1, col=1)

        if save_html:
            out = self.output_dir / f"{sample_name}_{hit_type.lower()}_alignment_quality.html"
            fig.write_html(str(out))
            logger.info("Saved alignment quality chart -> %s", out)

        return fig

    # ── JSON export ───────────────────────────────────────────────────────

    def export_json(
        self,
        profiles: list[CategoryProfile],
        sample_name: str,
    ) -> Path:
        data = {
            "sample": sample_name,
            "categories": [
                {
                    "category": p.category,
                    "unique_genes": p.unique_genes,
                    "hit_count": p.hit_count,
                    "mean_identity": round(p.mean_identity, 2),
                    "mean_coverage": round(p.mean_coverage, 2),
                    "tpm": round(p.tpm, 2),
                    "risk_level": p.risk_level,
                    "gene_names": p.gene_names,
                }
                for p in profiles
            ],
        }
        out = self.output_dir / f"{sample_name}_classification.json"
        with open(out, "w") as fh:
            json.dump(data, fh, indent=2)
        logger.info("Exported JSON -> %s", out)
        return out
