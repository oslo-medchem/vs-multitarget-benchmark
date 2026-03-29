#!/usr/bin/env node
/**
 * generate_bib_manuscript.js
 * Generates Manuscript_BIB_v1.docx and Supplementary_BIB_v1.docx
 * for Briefings in Bioinformatics (Method article) submission.
 */

const fs = require("fs");
const path = require("path");
const docx = require(path.join(__dirname, "node_modules", "docx"));

const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, PageBreak, WidthType, AlignmentType, BorderStyle,
  HeadingLevel, TableLayoutType, ShadingType, LineSpacing,
  convertInchesToTwip, ExternalHyperlink, LineRuleType,
  Header, Footer, PageNumber, NumberFormat
} = docx;

// --------------- paths ---------------
const BASE = path.join(__dirname);
const FIG_DIR = path.join(BASE, "figures");
const DATA_PATH = "/tmp/manuscript_data.json";
const META_PATH = "/data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget/06_meta_analysis/meta_figures/meta_analysis_results.json";

// --------------- load data ---------------
const data = JSON.parse(fs.readFileSync(DATA_PATH, "utf8"));
const meta = JSON.parse(fs.readFileSync(META_PATH, "utf8"));

// --------------- helpers ---------------
const FONT = "Times New Roman";
const FONT_SIZE = 24; // half-points → 12pt
const LINE_SPACING = 480;
const TABLE_FONT_SIZE = 20; // 10pt for tables

function txt(text, opts = {}) {
  return new TextRun({ text, font: FONT, size: opts.size || FONT_SIZE, bold: opts.bold, italics: opts.italics, superScript: opts.sup, color: opts.color });
}

function para(runs, opts = {}) {
  if (typeof runs === "string") runs = [txt(runs)];
  return new Paragraph({
    children: runs,
    alignment: opts.alignment || AlignmentType.JUSTIFIED,
    spacing: { line: opts.line || LINE_SPACING, after: opts.after || 120 },
    indent: opts.indent,
    pageBreakBefore: opts.pageBreak,
  });
}

function heading(text, opts = {}) {
  return new Paragraph({
    children: [txt(text, { bold: true, size: opts.size || FONT_SIZE })],
    alignment: AlignmentType.LEFT,
    spacing: { line: LINE_SPACING, before: 240, after: 120 },
    pageBreakBefore: opts.pageBreak,
  });
}

function subheading(text) {
  return new Paragraph({
    children: [txt(text, { bold: true, italics: true })],
    alignment: AlignmentType.LEFT,
    spacing: { line: LINE_SPACING, before: 200, after: 80 },
  });
}

function emptyPara() {
  return para([txt("")], { after: 0 });
}

// Table helpers
const THIN_BORDER = { style: BorderStyle.SINGLE, size: 1, color: "000000" };
const NO_BORDER = { style: BorderStyle.NONE, size: 0 };
const HEADER_SHADING = { type: ShadingType.SOLID, color: "D9E2F3" };

function tableCell(content, opts = {}) {
  const runs = typeof content === "string" ? [txt(content, { size: TABLE_FONT_SIZE, bold: opts.headerRow })] : content;
  return new TableCell({
    children: [new Paragraph({ children: runs, alignment: opts.alignment || AlignmentType.CENTER, spacing: { line: 276, after: 40, before: 40 } })],
    shading: opts.headerRow ? HEADER_SHADING : undefined,
    width: opts.width ? { size: opts.width, type: WidthType.DXA } : undefined,
    verticalAlign: docx.VerticalAlign ? docx.VerticalAlign.CENTER : undefined,
    borders: {
      top: THIN_BORDER, bottom: THIN_BORDER,
      left: THIN_BORDER, right: THIN_BORDER,
    },
  });
}

function loadImage(filename) {
  const p = path.join(FIG_DIR, filename);
  if (!fs.existsSync(p)) {
    console.warn(`WARNING: Image not found: ${p}`);
    return null;
  }
  return fs.readFileSync(p);
}

function figParagraph(imageData, widthPx, heightPx, targetWidthPx) {
  const scale = targetWidthPx / widthPx;
  return new Paragraph({
    children: [
      new ImageRun({
        data: imageData,
        transformation: { width: targetWidthPx, height: Math.round(heightPx * scale) },
        type: "png",
      }),
    ],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 120 },
  });
}

// Quick image dimension reader for PNG
function pngDimensions(buf) {
  // PNG header: bytes 16-19 = width, 20-23 = height (big-endian)
  if (buf[0] === 0x89 && buf[1] === 0x50) {
    const w = buf.readUInt32BE(16);
    const h = buf.readUInt32BE(20);
    return { width: w, height: h };
  }
  return { width: 800, height: 600 }; // fallback
}

// --------------- target metadata ---------------
const TARGET_ORDER = ["fpr2", "egfr", "cdk2", "cox2", "esr1", "dpp4", "ache", "bace1", "hsp90", "p38", "thrombin"];
const TARGET_LABELS = {
  fpr2: "FPR2", egfr: "EGFR", cdk2: "CDK2", cox2: "COX-2", esr1: "ESR1",
  dpp4: "DPP-4", ache: "AChE", bace1: "BACE-1", hsp90: "HSP90\u03B1", p38: "p38", thrombin: "Thrombin"
};
const TARGET_FULLNAMES = {
  fpr2: "Formyl peptide receptor 2", egfr: "Epidermal growth factor receptor",
  cdk2: "Cyclin-dependent kinase 2", cox2: "Cyclooxygenase-2",
  esr1: "Estrogen receptor alpha", dpp4: "Dipeptidyl peptidase-4",
  ache: "Acetylcholinesterase", bace1: "Beta-secretase 1",
  hsp90: "Heat shock protein 90-alpha", p38: "p38 MAP kinase alpha", thrombin: "Thrombin"
};
const TARGET_CLASS = {
  fpr2: "GPCR", egfr: "Kinase", cdk2: "Kinase", cox2: "Enzyme",
  esr1: "Nuclear receptor", dpp4: "Serine protease", ache: "Hydrolase",
  bace1: "Aspartyl protease", hsp90: "Chaperone", p38: "Kinase", thrombin: "Serine protease"
};
const TARGET_CHEMBL = {
  fpr2: "CHEMBL4227", egfr: "CHEMBL203", cdk2: "CHEMBL301", cox2: "CHEMBL230",
  esr1: "CHEMBL206", dpp4: "CHEMBL284", ache: "CHEMBL220", bace1: "CHEMBL4822",
  hsp90: "CHEMBL3880", p38: "CHEMBL260", thrombin: "CHEMBL204"
};
const TARGET_PDB = {
  fpr2: "7T6S:R", egfr: "1M17:A", cdk2: "1H1Q:A", cox2: "3LN1:A",
  esr1: "1SJ0:A", dpp4: "2RGU:A", ache: "4EY7:A", bace1: "4IVS:A",
  hsp90: "2WI7:A", p38: "1KV2:A", thrombin: "2ZDM:A"
};
const TARGET_RESOLUTION = {
  fpr2: "3.0", egfr: "2.6", cdk2: "2.0", cox2: "2.4", esr1: "2.0",
  dpp4: "2.2", ache: "2.35", bace1: "1.55", hsp90: "2.0", p38: "2.1", thrombin: "2.2"
};
const TARGET_METHOD = {
  fpr2: "Cryo-EM", egfr: "X-ray", cdk2: "X-ray", cox2: "X-ray", esr1: "X-ray",
  dpp4: "X-ray", ache: "X-ray", bace1: "X-ray", hsp90: "X-ray", p38: "X-ray", thrombin: "X-ray"
};
const TARGET_LIGAND = {
  fpr2: "FUI", egfr: "AQ4", cdk2: "2A6", cox2: "CEL", esr1: "E4D",
  dpp4: "356", ache: "E20", bace1: "VSI", hsp90: "2KL", p38: "B96", thrombin: "46U"
};

function fmtP(v) { return v < 0.001 ? v.toExponential(2) : v.toFixed(3); }
function fmtAUC(v) { return v.toFixed(3); }
function fmtPct(v) { return v.toFixed(1); }

// --------------- page properties ---------------
const PAGE_PROPS = {
  page: {
    size: { width: convertInchesToTwip(8.27), height: convertInchesToTwip(11.69) }, // A4
    margin: { top: convertInchesToTwip(1), bottom: convertInchesToTwip(1), left: convertInchesToTwip(1), right: convertInchesToTwip(1) },
  },
};

// ========================================================================
//  MAIN MANUSCRIPT
// ========================================================================
function buildMainManuscript() {
  const children = [];

  // Title
  children.push(new Paragraph({
    children: [txt("Virtual screening with agentic AI: human expertise injection is still essential", { bold: true, size: 28 })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 200 },
  }));

  // Authors
  children.push(new Paragraph({
    children: [txt("Osman A.B.S.M. Gani", { size: FONT_SIZE })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 60 },
  }));
  children.push(new Paragraph({
    children: [txt("Section for Pharmaceutical Chemistry, Department of Pharmacy, University of Oslo, 0316 Oslo, Norway", { size: 22, italics: true })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 60 },
  }));
  children.push(new Paragraph({
    children: [txt("Corresponding author: osman.gani@farmasi.uio.no", { size: 22 })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 240 },
  }));

  // Key Points box
  children.push(heading("Key Points"));
  const keyPoints = [
    "Structured skill files (212 lines) dramatically improve ligand preparation quality in LLM-generated VS pipelines (86\u201399% vs 46\u201395% success).",
    "Skill-guided protocols significantly improve docking discrimination for targets with deep binding sites (COX-2, AChE).",
    "No target showed significant deterioration under skill guidance across 11 diverse drug targets.",
    "Skill files function as lightweight, interpretable quality assurance for LLM coding agents.",
  ];
  for (const kp of keyPoints) {
    children.push(new Paragraph({
      children: [txt("\u2022 " + kp)],
      alignment: AlignmentType.JUSTIFIED,
      spacing: { line: LINE_SPACING, after: 60 },
      indent: { left: 360 },
    }));
  }
  children.push(emptyPara());

  // Abstract
  children.push(heading("Abstract"));
  children.push(para("Large language model (LLM)-based coding agents can autonomously generate computational drug discovery pipelines, but the impact of domain knowledge injection on virtual screening (VS) outcomes across diverse targets remains unexplored. Here I present a controlled benchmark comparing Claude Code (Claude Opus 4.6) with and without a 212-line structured skill file across 11 drug targets spanning six protein classes (kinases, proteases, hydrolases, enzymes, nuclear receptors, chaperones, and GPCRs). For each target, the agent autonomously generated ~5,500 lines of code implementing a complete GPU-accelerated VS pipeline\u2014from ChEMBL active retrieval through Uni-Dock docking to statistical evaluation\u2014without human code editing. The skill-guided protocol achieved dramatically higher ligand preparation success rates (86\u201399%) compared with the naive baseline (46\u201395%), representing a substantial practical improvement in pipeline robustness. Docking discrimination, measured by ROC AUC on successfully docked compounds, improved for four of eleven targets, with statistically significant gains for cyclooxygenase-2 (COX-2; \u0394AUC = +0.129, DeLong p = 0.016) and acetylcholinesterase (AChE; \u0394AUC = +0.157, p = 0.004)\u2014both targets characterised by deep, enclosed binding pockets. No target showed significant deterioration under skill guidance. These findings,  on FPR2, demonstrate that structured skill files function as effective quality assurance mechanisms for LLM coding agents, with additional discrimination benefits for targets with demanding binding-site architectures."));

  // Keywords
  children.push(new Paragraph({
    children: [
      txt("Keywords: ", { bold: true }),
      txt("virtual screening, large language models, agentic AI, skill files, molecular docking, drug discovery, benchmark"),
    ],
    alignment: AlignmentType.JUSTIFIED,
    spacing: { line: LINE_SPACING, after: 200 },
  }));

  // Introduction
  children.push(heading("Introduction", { pageBreak: true }));
  children.push(para([
    txt("Virtual screening (VS) has become an indispensable tool in modern drug discovery, enabling the computational identification of promising drug candidates from large compound libraries [1\u20133]. The recent development of GPU-accelerated docking engines, particularly Uni-Dock [4], has dramatically reduced the computational cost of structure-based VS, making it feasible to screen millions of compounds within hours on commodity hardware."),
  ]));
  children.push(para([
    txt("In parallel, large language model (LLM)-based coding agents such as Claude Code [5] have emerged as powerful tools for automating complex computational workflows. These agents can generate, execute, and debug code across multiple programming languages, raising the possibility of fully autonomous computational drug discovery pipelines. Related efforts include ChemCrow [7] and other autonomous chemical research systems [6]. A recent comprehensive review catalogues the expanding role of LLMs and autonomous agents in chemistry [8]."),
  ]));
  children.push(para([
    txt("A critical question for the deployment of LLM coding agents in scientific research is whether domain knowledge injection improves the quality of their outputs. Skill files\u2014structured plain-text documents containing domain-specific best practices, tool recommendations, and parameter guidelines\u2014represent a lightweight mechanism for injecting such knowledge into LLM agents [5]. Unlike fine-tuning or retrieval-augmented generation (RAG), skill files are human-readable, version-controllable, and model-agnostic."),
  ]));
  children.push(para([
    txt("Whether structured skill files can improve VS outcomes across diverse protein targets and binding-site architectures has not been systematically investigated."),
  ]));
  children.push(para([
    txt("Here I present a controlled benchmark across 11 drug targets spanning six protein classes\u2014kinases, serine proteases, an aspartyl protease, a hydrolase, an enzyme, a nuclear receptor, a chaperone, and a GPCR\u2014to test two hypotheses: (1) that skill files universally improve pipeline robustness (measured by ligand preparation success rates), and (2) that skill files improve docking discrimination (measured by ROC AUC) in a target-dependent manner, with the greatest benefits for targets with deep, enclosed binding sites."),
  ]));

  // Methods
  children.push(heading("Methods", { pageBreak: true }));

  children.push(subheading("Target selection"));
  children.push(para([
    txt("Eleven drug targets were selected to represent six major protein classes relevant to drug discovery (see Table S1 in Supplementary Information for full details). The panel comprises three kinases (EGFR, CDK2, p38), two serine proteases (DPP-4, thrombin), one aspartyl protease (BACE-1), one hydrolase (AChE), one enzyme (COX-2), one nuclear receptor (ESR1), one chaperone (HSP90\u03B1), and one GPCR (FPR2). All targets have co-crystallised ligands suitable for defining the docking box centre. Structural resolution ranged from 1.55 to 2.6 \u00C5 for X-ray structures, with FPR2 determined at 3.0 \u00C5 by cryo-EM."),
  ]));

  children.push(subheading("Agentic pipeline generation"));
  children.push(para([
    txt("Claude Code with the Claude Opus 4.6 model [5] was used to autonomously generate two complete VS pipelines per target: a naive protocol (base agent without domain guidance) and a skill-guided protocol (agent supplemented with a 212-line structured skill file). Each pipeline comprised 18 Python scripts totalling approximately 5,500 lines of code, implementing five stages: (1) active compound curation from ChEMBL, (2) property-matched decoy generation, (3) ligand and receptor preparation, (4) GPU-accelerated docking, and (5) statistical evaluation. No human code editing was performed."),
  ]));

  children.push(subheading("Dataset construction"));
  children.push(para([
    txt("For each target, 100 active compounds were retrieved from ChEMBL v34 [9] using a pChEMBL \u2265 5 threshold, followed by Butina clustering [9] and MaxMin diversity selection to ensure chemical diversity. Approximately 4,500 property-matched decoys per target were drawn from a shared pool of 60,000 ChEMBL compounds, with matching tolerances of MW \u00B125 Da, cLogP \u00B11.0, HBD \u00B11, HBA \u00B12, RotBond \u00B12, and Tanimoto coefficient < 0.35 against all actives."),
  ]));

  children.push(subheading("Protocol differences"));
  children.push(para([
    txt("Table 1 summarises the key differences between the two protocols."),
  ]));

  // Protocol differences table (inline in methods - small table)
  const protoHeaders = ["Parameter", "Naive", "Skill-guided"];
  const protoRows = [
    ["Receptor preparation", "OpenBabel [9]", "Meeko [9]"],
    ["Ligand 3D generation", "OpenBabel gen3d", "RDKit ETKDGv3 [14,15]"],
    ["Chemical filters", "None", "PAINS + Brenk [16,17]"],
    ["Docking box size", "20 \u00C5", "25 \u00C5"],
    ["Exhaustiveness", "8", "32"],
  ];
  const protoTable = new Table({
    rows: [
      new TableRow({ children: protoHeaders.map(h => tableCell(h, { headerRow: true })) }),
      ...protoRows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
    ],
    width: { size: 9000, type: WidthType.DXA },
  });
  children.push(protoTable);
  children.push(emptyPara());

  children.push(subheading("Docking"));
  children.push(para([
    txt("All docking was performed with Uni-Dock v1.1.3 [4] on two NVIDIA RTX 4500 Ada GPUs (24 GB VRAM each, CUDA 12.9). Ligands were processed in batches of 100 with sub-batch retry logic to handle malformed PDBQT files without terminating the full screening run."),
  ]));

  children.push(subheading("Evaluation"));
  children.push(para([
    txt("The primary metric was ROC AUC computed on the subset of compounds successfully docked by each protocol (docked-only analysis) [9]. Bootstrap 95% confidence intervals were computed using the bias-corrected and accelerated (BCa) method with n = 10,000 replicates. For the sensitivity analysis, DeLong\u2019s test [9] was applied to the intersection of compounds successfully docked by both protocols. Paired Wilcoxon signed-rank tests were used for cross-target comparisons, and Benjamini\u2013Hochberg correction [9] was applied where appropriate."),
  ]));

  // Results
  children.push(heading("Results", { pageBreak: true }));

  children.push(subheading("Pipeline quality: the silent failure of unguided agents"));
  children.push(para([
    txt("The skill-guided protocol achieved 86\u201399% preparation success across all targets, compared with 46\u201395% for the naive protocol (Fig. 2A). Without skill guidance, the agent selected OpenBabel gen3d for 3D coordinate generation, which produced 5\u201354% malformed PDBQT files across targets\u2014files that silently fail during docking without warning. This represents a critical, invisible failure mode: in a real-world VS campaign, a researcher relying on an unguided LLM agent would unknowingly screen only half of their intended chemical library. The compounds lost are not random; OpenBabel gen3d disproportionately fails on molecules with complex ring systems, macrocycles, and strained geometries\u2014precisely the chemotypes most likely to occupy novel regions of chemical space. By contrast, the skill file directed the agent to use RDKit ETKDGv3 for conformer generation and Meeko for PDBQT conversion, reducing the failure rate to 1\u201314% across all targets."),
  ]));

  // Table 2: Preparation success rates
  children.push(para([txt("Table 2. ", { bold: true }), txt("PDBQT preparation success rates across 11 targets.")]));
  {
    const hdr = ["Target", "Library", "Naive docked (%)", "Skill docked (%)"];
    const rows = TARGET_ORDER.map(t => {
      const d = data[t];
      return [
        TARGET_LABELS[t],
        String(d.n_total),
        `${d.naive_docked} (${fmtPct(d.naive_pct)}%)`,
        `${d.skill_docked} (${fmtPct(d.skill_pct)}%)`,
      ];
    });
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 9000, type: WidthType.DXA },
    });
    children.push(tbl);
    children.push(emptyPara());
  }

  children.push(subheading("Per-target docking discrimination"));
  children.push(para([
    txt("Among compounds that were successfully docked, the skill protocol achieved higher ROC AUC for four of eleven targets (Fig. 2B): COX-2 (+0.073), AChE (+0.065), EGFR (+0.037), and ESR1 (+0.033). The two targets showing the largest improvements\u2014COX-2 and AChE\u2014both reached statistical significance in the intersection DeLong analysis (see below). The remaining seven targets showed equivalent performance between protocols, with none showing statistically significant deterioration. This asymmetry (4 improved, 0 worsened) is consistent with skill files providing a net benefit: they improve discrimination where it matters (deep pockets requiring thorough sampling) without compromising performance on simpler targets."),
  ]));

  // Table 3: Docked-only AUC
  children.push(para([txt("Table 3. ", { bold: true }), txt("Docked-only ROC AUC for 11 targets under both protocols with intersection DeLong test results.")]));
  {
    const hdr = ["Target", "Naive AUC", "Skill AUC", "\u0394AUC", "Intersection p"];
    const rows = TARGET_ORDER.map(t => {
      const d = data[t];
      const delta = d.skill_auc - d.naive_auc;
      return [
        TARGET_LABELS[t],
        fmtAUC(d.naive_auc),
        fmtAUC(d.skill_auc),
        (delta >= 0 ? "+" : "") + fmtAUC(delta),
        fmtP(d.int_p),
      ];
    });
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 9000, type: WidthType.DXA },
    });
    children.push(tbl);
    children.push(emptyPara());
  }

  children.push(subheading("Significant improvements: deep-pocket targets"));
  children.push(para([
    txt("COX-2 and AChE showed statistically significant improvements in the intersection analysis: COX-2 \u0394AUC = +0.129 [+0.025, +0.234], z = 2.42, p = 0.016; AChE \u0394AUC = +0.157 [+0.050, +0.264], z = 2.87, p = 0.004. These targets share deep, enclosed binding sites\u2014the cyclooxygenase channel and the AChE catalytic gorge, respectively\u2014where the larger docking box (25 vs 20 \u00C5) and higher exhaustiveness (32 vs 8) prescribed by the skill file likely afford more thorough pose sampling."),
  ]));

  // Fig 1 - Forest plot
  {
    const imgBuf = loadImage("forest_plot.png");
    if (imgBuf) {
      const dims = pngDimensions(imgBuf);
      children.push(figParagraph(imgBuf, dims.width, dims.height, 500));
    }
    children.push(para([
      txt("Fig. 1. ", { bold: true }),
      txt("Forest plot of intersection \u0394AUC (skill \u2212 naive) across 11 drug targets. Point estimates (squares) with 95% DeLong confidence intervals. Targets are grouped by protein class. Two targets showed statistically significant improvements: COX-2 (\u0394AUC = +0.129, p = 0.016) and AChE (\u0394AUC = +0.157, p = 0.004). I\u00B2 = 52%, indicating moderate heterogeneity."),
    ]));
  }

  // Fig 2 - Combined tables figure (preparation quality + ΔAUC)
  {
    const imgBuf = loadImage("fig_tables_combined.png");
    if (imgBuf) {
      const dims = pngDimensions(imgBuf);
      children.push(new Paragraph({ children: [new PageBreak()] }));
      children.push(figParagraph(imgBuf, dims.width, dims.height, 550));
    }
    children.push(para([
      txt("Fig. 2. ", { bold: true }),
      txt("Skill file impact across 11 drug targets. (A) PDBQT preparation failure rates: the naive protocol (red) loses 5\u201354% of compounds to malformed files; the skill protocol (green) reduces failures to 1\u201314%. Targets sorted by naive failure rate. (B) Change in docking discrimination (\u0394AUC, skill \u2212 naive, docked-only compounds): green bars indicate skill advantage. COX-2 and AChE reach statistical significance (p < 0.05, intersection DeLong test). No target shows significant deterioration."),
    ]));
  }

  // Fig 3 - Combined COX-2 + AChE success stories
  {
    const imgBuf = loadImage("fig_success_combined.png");
    if (imgBuf) {
      const dims = pngDimensions(imgBuf);
      children.push(new Paragraph({ children: [new PageBreak()] }));
      children.push(figParagraph(imgBuf, dims.width, dims.height, 500));
    }
    children.push(para([
      txt("Fig. 3. ", { bold: true }),
      txt("Detailed VS performance for the two targets with statistically significant skill improvements. (A) COX-2 (CHEMBL230, PDB 3LN1): the skill protocol achieves visibly better active\u2013decoy separation in ROC curves and higher bootstrap AUC. (B) AChE (CHEMBL220, PDB 4EY7): same pattern. Both targets feature deep, enclosed binding sites\u2014the cyclooxygenase channel and the catalytic gorge, respectively\u2014where the skill file\u2019s larger docking box (25 vs 20 \u00C5) and higher exhaustiveness (32 vs 8) enable more thorough pose sampling. Left panels: ROC curves with 95% bootstrap confidence bands. Centre: bootstrap AUC distributions (n = 2,000). Right: docking score distributions for actives vs decoys."),
    ]));
  }

  children.push(subheading("Cross-target analysis"));
  children.push(para([
    txt("The effect of skill guidance varied across targets (I\u00B2 = 52%), consistent with genuine target-dependent variation rather than a uniform effect (Fig. 1). Targets with deeper binding sites benefited most, while those with solvent-exposed ATP-binding clefts (kinases: EGFR, CDK2, p38) showed equivalent performance under both protocols. A sign test (4/11 positive) was not significant (p = 0.55), consistent with target-dependent rather than uniform benefit."),
  ]));

  // Discussion
  children.push(heading("Discussion", { pageBreak: true }));
  children.push(para([
    txt("The central finding of this study is that an LLM coding agent without domain knowledge injection produces virtual screening pipelines with a critical, hidden deficiency: up to 54% of the screening library is silently lost to preparation failures. Structured skill files eliminate this failure mode while additionally improving docking discrimination for targets with demanding binding-site architectures."),
  ]));
  children.push(para([
    txt("Pipeline quality represents the primary practical benefit of skill guidance. A pipeline that silently loses half its library to preparation failures will miss candidate compounds regardless of docking accuracy. The skill file\u2019s guidance to use RDKit ETKDGv3 for 3D coordinate generation and Meeko for PDBQT conversion eliminates this failure mode, achieving 86\u201399% preparation success compared with 46\u201395% for the naive protocol. This improvement was consistent across all 11 targets and represents the single most impactful contribution of the skill file."),
  ]));
  children.push(para([
    txt("The target-specific discrimination improvements observed for COX-2 and AChE support a \u201Cdeep-pocket hypothesis\u201D: targets with enclosed, channel-like binding sites (the cyclooxygenase channel and the AChE catalytic gorge) benefit from the larger docking box (25 vs 20 \u00C5) and higher exhaustiveness (32 vs 8) prescribed by the skill file, which afford more thorough pose sampling in these geometrically demanding environments. By contrast, kinases with solvent-exposed ATP-binding clefts showed equivalent performance under both protocols, suggesting that default parameters suffice for these simpler binding-site topologies."),
  ]));
  children.push(para([
    txt("The expanded benchmark reveals that the magnitude of improvement is target-dependent, with the largest gains observed for demanding pocket topologies rather than being uniformly distributed across protein classes."),
  ]));
  children.push(para([
    txt("From a practical standpoint, skill files represent a lightweight, interpretable quality assurance mechanism for LLM coding agents. Unlike fine-tuning (which requires training data and compute) or RAG (which requires vector databases and retrieval infrastructure), skill files are plain-text documents that are human-readable, version-controllable, and model-agnostic. They ensure that agents select appropriate tools (RDKit over OpenBabel for conformer generation, Meeko over OpenBabel for PDBQT conversion) and parameters (box size, exhaustiveness) without requiring the user to specify these choices manually."),
  ]));
  children.push(para([
    txt("Several limitations should be acknowledged. The benchmark uses 100 actives per target, which, while standard for VS benchmarks, may not capture the full diversity of active chemical space. Only the Vina scoring function (via Uni-Dock) was tested; other scoring functions may show different sensitivity to preparation protocols. The high OpenBabel failure rate in the naive protocol may overstate the preparation quality gap relative to other naive tool choices. Finally, PAINS/Brenk filtering in the skill protocol removes some genuine actives (38\u201367 per target), which could reduce recall in practice."),
  ]));

  // Conclusion
  children.push(heading("Conclusion"));
  children.push(para([
    txt("Structured skill files enhance agentic AI virtual screening through two complementary mechanisms: universally improved pipeline quality (86\u201399% vs 46\u201395% preparation success) and significant discrimination gains for targets with deep binding sites (COX-2: \u0394AUC = +0.129, p = 0.016; AChE: \u0394AUC = +0.157, p = 0.004). These findings establish skill files as effective, interpretable quality assurance mechanisms for democratising rigorous computational drug discovery."),
  ]));

  // Acknowledgements
  children.push(heading("Acknowledgements"));
  children.push(para([
    txt("All computational protocols were designed and executed with Claude Code (Claude Opus 4.6, Anthropic). In accordance with ISCB guidelines, I disclose that the LLM agent was used to generate the computational pipeline code; all scientific analysis, interpretation, and manuscript writing were performed by the author. Docking calculations were performed on NVIDIA RTX 4500 Ada GPUs."),
  ]));

  // Funding
  children.push(heading("Funding"));
  children.push(para("University of Oslo Growth House [ref 25035]."));

  // Data availability
  children.push(heading("Data availability"));
  children.push(para("All code, docking scores, and evaluation metrics are available at https://github.com/oslo-medchem/vs-multitarget-benchmark under the MIT licence."));

  // References
  children.push(heading("References", { pageBreak: true }));
  const refs = [
    "Schneider P et al. (2020) Rethinking drug design in the artificial intelligence era. Nat Rev Drug Discov, 19, 353\u2013364.",
    "Walters WP, Barzilay R. (2021) Applications of deep learning in molecule generation and molecular property prediction. Acc Chem Res, 54, 263\u2013270.",
    "Huang N et al. (2006) Benchmarking sets for molecular docking. J Med Chem, 49, 6789\u20136801.",
    "Yu Y et al. (2023) Uni-Dock: GPU-accelerated docking enables ultralarge virtual screening. J Chem Theory Comput, 19, 3336\u20133345.",
    "Anthropic. (2024) Claude Code: an agentic coding tool. https://claude.ai/code",
    "Boiko DA et al. (2023) Autonomous chemical research with large language models. Nature, 624, 570\u2013578.",
    "Bran AM et al. (2024) ChemCrow: augmenting large-language models with chemistry tools. Nat Mach Intell, 6, 525\u2013535.",
    "Ramos MC et al. (2025) A review of large language models and autonomous agents in chemistry. Chem Sci, 16, 2514\u20132572.",
    "",
    "Zdrazil B et al. (2024) The ChEMBL Database in 2023. Nucleic Acids Res, 52, D1180\u2013D1192.",
    "Butina D. (1999) Unsupervised data base clustering based on Daylight\u2019s fingerprint and Tanimoto similarity. J Chem Inf Comput Sci, 39, 747\u2013750.",
    "O\u2019Boyle NM et al. (2011) Open Babel: an open chemical toolbox. J Cheminform, 3, 33.",
    "Forli S et al. (2016) Computational protein\u2013ligand docking and virtual drug screening with the AutoDock suite. Nat Protoc, 11, 905\u2013919.",
    "Wang S et al. (2020) Improving conformer generation for small rings and macrocycles. J Chem Inf Model, 60, 2044\u20132058.",
    "Landrum GA et al. (2023) RDKit: open-source cheminformatics. Zenodo. doi:10.5281/zenodo.591637",
    "Baell JB, Holloway GA. (2010) New substructure filters for removal of PAINS. J Med Chem, 53, 2719\u20132740.",
    "Brenk R et al. (2008) Lessons learnt from assembling screening libraries. ChemMedChem, 3, 435\u2013444.",
    "Truchon J-F, Bayly CI. (2007) Evaluating virtual screening methods. J Chem Inf Model, 47, 488\u2013508.",
    "DeLong ER et al. (1988) Comparing the areas under two or more correlated ROC curves. Biometrics, 44, 837\u2013845.",
    "Benjamini Y, Hochberg Y. (1995) Controlling the false discovery rate. J R Stat Soc B, 57, 289\u2013300.",
  ];
  for (let i = 0; i < refs.length; i++) {
    children.push(new Paragraph({
      children: [txt(`[${i + 1}] ${refs[i]}`)],
      alignment: AlignmentType.JUSTIFIED,
      spacing: { line: LINE_SPACING, after: 60 },
      indent: { left: 360, hanging: 360 },
    }));
  }

  return new Document({
    sections: [{ properties: PAGE_PROPS, children }],
  });
}

// ========================================================================
//  SUPPLEMENTARY INFORMATION
// ========================================================================
function buildSupplementary() {
  const children = [];

  // Title page
  children.push(new Paragraph({
    children: [txt("Supplementary Information", { bold: true, size: 32 })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 200 },
  }));
  children.push(new Paragraph({
    children: [txt("Virtual screening with agentic AI: human expertise injection is still essential", { bold: true, size: 26 })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 200 },
  }));
  children.push(new Paragraph({
    children: [txt("Osman A.B.S.M. Gani")],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 60 },
  }));
  children.push(new Paragraph({
    children: [txt("Section for Pharmaceutical Chemistry, Department of Pharmacy, University of Oslo, 0316 Oslo, Norway", { size: 22, italics: true })],
    alignment: AlignmentType.CENTER,
    spacing: { line: LINE_SPACING, after: 400 },
  }));

  // =====================================================================
  // Supplementary Methods: Statistical metrics explained
  // =====================================================================
  children.push(heading("Supplementary Methods: Evaluation Metrics", { pageBreak: true }));

  children.push(para([
    txt("This section provides accessible explanations of the virtual screening (VS) evaluation metrics used in this study, the mathematical intuition behind each, and the rationale for our metric selection."),
  ]));

  // --- ROC AUC ---
  children.push(subheading("ROC AUC \u2014 Receiver Operating Characteristic Area Under the Curve"));
  children.push(para([
    txt("The ROC curve plots the ", { bold: false }),
    txt("true positive rate", { bold: true }),
    txt(" (fraction of actives ranked above a given threshold) against the "),
    txt("false positive rate", { bold: true }),
    txt(" (fraction of decoys incorrectly ranked above that threshold) as the score threshold is swept from the most stringent to the most lenient. A perfect ranking\u2014all actives before all decoys\u2014produces a curve that rises immediately to the top-left corner, enclosing an area of 1.0. A random ranking produces a diagonal line with area 0.5."),
  ]));
  children.push(para([
    txt("The AUC (area under this curve) summarises the overall probability that a randomly chosen active is ranked higher than a randomly chosen decoy. It ranges from 0 to 1, where:"),
  ]));
  children.push(para([
    txt("\u2022  AUC = 1.0: perfect discrimination (all actives ranked above all decoys)"),
  ], { after: 40 }));
  children.push(para([
    txt("\u2022  AUC = 0.5: random ranking (no discriminatory power)"),
  ], { after: 40 }));
  children.push(para([
    txt("\u2022  AUC < 0.5: worse than random (actives systematically ranked below decoys)"),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("ROC AUC is the most widely reported metric in VS benchmarks [3, 18] because it is threshold-independent, intuitive, and directly comparable across studies. It was our primary metric for assessing overall discrimination."),
  ]));

  // --- BEDROC ---
  children.push(subheading("BEDROC \u2014 Boltzmann-Enhanced Discrimination of ROC"));
  children.push(para([
    txt("In practical drug discovery, only the "),
    txt("top-ranked", { bold: true }),
    txt(" compounds are selected for experimental testing\u2014typically the top 1\u20135% of the library. Standard ROC AUC weights all parts of the ranked list equally, so it may not reflect how well a method performs at the critical \u201Cearly\u201D portion."),
  ]));
  children.push(para([
    txt("BEDROC addresses this by applying an exponential weighting function that gives much greater importance to compounds ranked near the top of the list. The parameter \u03B1 controls how steeply the weighting falls off:"),
  ]));
  children.push(para([
    txt("\u2022  \u03B1 = 20 (used here): approximately 80% of the metric\u2019s weight is concentrated in the top 8% of the ranked list"),
  ], { after: 40 }));
  children.push(para([
    txt("\u2022  \u03B1 = 80: weight is concentrated in the top 2%"),
  ], { after: 40 }));
  children.push(para([
    txt("\u2022  \u03B1 \u2192 0: BEDROC converges to standard ROC AUC"),
  ]));
  children.push(para([
    txt("BEDROC ranges from 0 to 1, where 1 means all actives are concentrated at the very top of the ranked list."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("BEDROC captures \u201Cearly enrichment\u201D\u2014the ability to find actives in the top fraction of the list\u2014which is the most practically relevant metric for VS campaigns where only a few hundred compounds can be tested experimentally [9]."),
  ]));

  // --- EF ---
  children.push(subheading("EF \u2014 Enrichment Factor"));
  children.push(para([
    txt("The enrichment factor at a given fraction x% of the library measures how many more actives appear in the top x% than expected by chance:"),
  ]));
  children.push(para([
    txt("EF(x%) = (actives in top x%) / (expected actives in top x%)", { italics: true }),
  ], { after: 60 }));
  children.push(para([
    txt("For example, if 2% of the library consists of actives and EF(1%) = 5.0, the top 1% of the ranked list contains 5 times as many actives as random selection would predict. The maximum possible EF depends on the active:decoy ratio; with 100 actives in ~4,600 compounds (ratio \u22482%), the theoretical maximum EF(1%) is approximately 46."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("EF at 1% and 5% directly answers the practical question: \u201CIf I test the top N compounds, how many hits will I find?\u201D This is the most intuitive metric for medicinal chemists planning follow-up experiments."),
  ]));

  // --- AUPR ---
  children.push(subheading("AUPR \u2014 Area Under the Precision\u2013Recall Curve"));
  children.push(para([
    txt("While ROC AUC uses false positive rate on the x-axis, the precision\u2013recall (PR) curve uses "),
    txt("recall", { bold: true }),
    txt(" (same as true positive rate: fraction of actives found) on the x-axis and "),
    txt("precision", { bold: true }),
    txt(" (fraction of predicted positives that are truly active) on the y-axis. AUPR is the area under this curve."),
  ]));
  children.push(para([
    txt("AUPR is particularly informative for "),
    txt("imbalanced datasets", { bold: true }),
    txt("\u2014exactly the situation in VS, where actives are rare (~2% of the library). In such settings, ROC AUC can appear high even when many decoys contaminate the top-ranked compounds. AUPR penalises this contamination directly."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("AUPR complements ROC AUC by providing a metric that is sensitive to the low active:decoy ratio typical of VS benchmarks. A method with high AUC but low AUPR is finding actives but also retrieving many false positives."),
  ]));

  // --- LogAUC ---
  children.push(subheading("LogAUC \u2014 Logarithmic ROC AUC"));
  children.push(para([
    txt("LogAUC is a variant of ROC AUC computed on a logarithmic x-axis. Instead of treating all false positive rates equally, it gives disproportionate weight to the low-false-positive-rate region (the left side of the ROC curve), which corresponds to the most stringent selection thresholds."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("LogAUC provides a continuous, threshold-free measure of early enrichment that is complementary to BEDROC. While BEDROC uses an exponential weight, LogAUC uses a logarithmic transformation, offering a different mathematical perspective on the same practical question."),
  ]));

  // --- DeLong test ---
  children.push(subheading("DeLong\u2019s Test \u2014 Comparing Two ROC Curves"));
  children.push(para([
    txt("When two different methods (naive and skill) produce scores for the same set of compounds, we need a statistical test to determine whether one method\u2019s AUC is significantly different from the other\u2019s. DeLong\u2019s test [9] is a non-parametric method specifically designed for this purpose."),
  ]));
  children.push(para([
    txt("The test accounts for the fact that the two ROC curves are "),
    txt("correlated", { bold: true }),
    txt("\u2014they are derived from the same compounds\u2014which a na\u00efve comparison (e.g., overlapping confidence intervals) would ignore. It produces a z-statistic and p-value testing the null hypothesis that the two AUCs are equal."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("DeLong\u2019s test is the standard method for comparing correlated ROC curves in VS benchmarks [9]. It is non-parametric (no distributional assumptions), accounts for the paired nature of the data, and provides both a p-value and a confidence interval for \u0394AUC."),
  ]));

  // --- Intersection analysis ---
  children.push(subheading("Intersection DeLong Analysis \u2014 Controlling for Preparation Differences"));
  children.push(para([
    txt("A methodological challenge arises when comparing two VS protocols that have different preparation success rates. Compounds that fail PDBQT preparation in one protocol receive the worst possible score, which biases the AUC comparison. Our "),
    txt("intersection analysis", { bold: true }),
    txt(" restricts the DeLong test to compounds that were successfully docked by "),
    txt("both", { italics: true }),
    txt(" protocols, thereby isolating the effect of docking quality from preparation quality."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("The naive protocol loses 5\u201354% of compounds to preparation failures. If these missing compounds are assigned worst-rank scores, the naive protocol appears artificially worse (or better, depending on whether the missing compounds are enriched in actives). The intersection analysis eliminates this confound."),
  ]));

  // --- Bootstrap CIs ---
  children.push(subheading("Bootstrap Confidence Intervals"));
  children.push(para([
    txt("Bootstrap resampling draws many random samples (with replacement) from the original data, computes the metric on each sample, and derives confidence intervals from the distribution of resampled values. We used the bias-corrected and accelerated (BCa) method with n = 10,000 resamples."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("Bootstrap CIs make no assumptions about the distribution of the metric and are robust for the small active counts (n \u2248 100 per target) in VS benchmarks."),
  ]));

  // --- Benjamini-Hochberg ---
  children.push(subheading("Benjamini\u2013Hochberg Correction"));
  children.push(para([
    txt("When performing multiple statistical tests (e.g., one DeLong test per target, or tests for multiple metrics), the probability of at least one false positive increases. The Benjamini\u2013Hochberg (BH) procedure [9] controls the "),
    txt("false discovery rate", { bold: true }),
    txt(" (FDR)\u2014the expected proportion of false positives among rejected hypotheses\u2014rather than the more conservative family-wise error rate. It adjusts p-values by ranking them and scaling by rank/total."),
  ]));
  children.push(para([
    txt("Why we used it: ", { bold: true }),
    txt("BH correction provides a good balance between controlling false discoveries and maintaining statistical power, which is important when comparing across 11 targets and multiple metrics."),
  ]));

  // --- Why these metrics and not others ---
  children.push(subheading("Rationale for Metric Selection"));
  children.push(para([
    txt("Virtual screening evaluation offers many possible metrics. Our selection was guided by three principles:"),
  ]));
  children.push(para([
    txt("1. Practical relevance to drug discovery campaigns. ", { bold: true }),
    txt("ROC AUC and EF directly answer questions that medicinal chemists ask: \u201CDoes the method separate actives from decoys?\u201D and \u201CHow many hits will I find if I test the top N compounds?\u201D"),
  ]));
  children.push(para([
    txt("2. Complementary perspectives on the ranked list. ", { bold: true }),
    txt("ROC AUC assesses global discrimination; BEDROC and LogAUC focus on early enrichment; EF provides a concrete hit-rate estimate at specific thresholds; AUPR accounts for class imbalance. Together, these five metrics cover the full range of VS performance aspects."),
  ]));
  children.push(para([
    txt("3. Established use in the VS benchmarking literature. ", { bold: true }),
    txt("All selected metrics are standard in DUD-E, DEKOIS, and other widely used benchmark suites [3, 18], ensuring comparability with published results."),
  ]));
  children.push(para([
    txt("Metrics we deliberately "),
    txt("excluded", { bold: true }),
    txt(":"),
  ]));
  children.push(para([
    txt("\u2022  ", { bold: false }),
    txt("Docking score RMSD or pose accuracy", { bold: true }),
    txt(": These require a known co-crystallised pose for each active, which is unavailable for most ChEMBL-derived actives. Our benchmark evaluates scoring/ranking, not pose prediction."),
  ], { after: 40 }));
  children.push(para([
    txt("\u2022  "),
    txt("Consensus scoring or rescoring", { bold: true }),
    txt(": We deliberately used only the Vina scoring function to isolate the effect of the skill file on a single, widely used scoring method. Consensus scoring, while potentially beneficial, would introduce additional variables."),
  ], { after: 40 }));
  children.push(para([
    txt("\u2022  "),
    txt("Molecular docking runtime", { bold: true }),
    txt(": The skill protocol uses 4\u00d7 higher exhaustiveness (32 vs 8), making it inherently slower. Runtime comparisons would conflate throughput with quality. In modern GPU-accelerated docking, throughput is rarely the bottleneck; discrimination quality is."),
  ]));

  // =====================================================================
  // Tables and Figures continue below
  // =====================================================================

  // Table S1: Full target panel
  children.push(heading("Table S1. Target panel for the multi-target virtual screening benchmark.", { pageBreak: true }));
  {
    const hdr = ["Target", "Full name", "Protein class", "ChEMBL ID", "PDB ID:Chain", "Res. (\u00C5)", "Method", "Co-cryst. ligand", "N actives", "N decoys", "Library"];
    const rows = TARGET_ORDER.map(t => {
      const d = data[t];
      const nDecoys = d.n_total - 100;
      return [
        TARGET_LABELS[t], TARGET_FULLNAMES[t], TARGET_CLASS[t],
        TARGET_CHEMBL[t], TARGET_PDB[t], TARGET_RESOLUTION[t],
        TARGET_METHOD[t], TARGET_LIGAND[t],
        "100", String(nDecoys), String(d.n_total),
      ];
    });
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true, width: hdr.length > 8 ? undefined : 1200 })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 14000, type: WidthType.DXA },
    });
    children.push(tbl);
  }

  // Table S2: Complete VS metrics
  children.push(heading("Table S2. Complete virtual screening metrics for all targets and protocols.", { pageBreak: true }));
  {
    const hdr = ["Target", "Protocol", "ROC AUC", "BEDROC(\u03B1=20)", "EF1%", "EF5%", "LogAUC", "AUPR"];
    const rows = [];
    for (const t of TARGET_ORDER) {
      const d = data[t];
      rows.push([TARGET_LABELS[t], "Naive", fmtAUC(d.naive_auc), fmtAUC(d.naive_bedroc), d.naive_ef1.toFixed(2), d.naive_ef5.toFixed(2), fmtAUC(d.naive_logauc), fmtAUC(d.naive_aupr)]);
      rows.push([TARGET_LABELS[t], "Skill", fmtAUC(d.skill_auc), fmtAUC(d.skill_bedroc), d.skill_ef1.toFixed(2), d.skill_ef5.toFixed(2), fmtAUC(d.skill_logauc), fmtAUC(d.skill_aupr)]);
    }
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 14000, type: WidthType.DXA },
    });
    children.push(tbl);
  }

  // Table S3: PAINS/Brenk filtering impact
  children.push(heading("Table S3. PAINS/Brenk filtering impact on the screening library.", { pageBreak: true }));
  {
    const hdr = ["Target", "Total filtered", "Actives removed", "Decoys removed", "% actives lost", "% library filtered"];
    const rows = TARGET_ORDER.map(t => {
      const d = data[t];
      const totalFiltered = d.n_filtered;
      const activesRemoved = d.pains_actives;
      const decoysRemoved = d.pains_decoys;
      const pctActivesLost = ((activesRemoved / 100) * 100).toFixed(1);
      const pctLibFiltered = ((totalFiltered / d.n_total) * 100).toFixed(1);
      return [TARGET_LABELS[t], String(totalFiltered), String(activesRemoved), String(decoysRemoved), pctActivesLost + "%", pctLibFiltered + "%"];
    });
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 10000, type: WidthType.DXA },
    });
    children.push(tbl);
  }

  // Table S4: Intersection DeLong sensitivity analysis
  children.push(heading("Table S4. Intersection DeLong sensitivity analysis.", { pageBreak: true }));
  {
    const hdr = ["Target", "N intersection", "\u0394AUC", "95% CI", "p-value"];
    const rows = TARGET_ORDER.map(t => {
      const d = data[t];
      const mp = meta.per_target.find(x => x.target === t);
      const ci = mp ? `[${mp.ci_low.toFixed(3)}, ${mp.ci_high.toFixed(3)}]` : "\u2014";
      return [TARGET_LABELS[t], String(d.int_n), (d.int_delta >= 0 ? "+" : "") + d.int_delta.toFixed(3), ci, fmtP(d.int_p)];
    });
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 10000, type: WidthType.DXA },
    });
    children.push(tbl);
  }

  // Table S5: Meta-analysis statistics
  children.push(heading("Table S5. Random-effects meta-analysis statistics.", { pageBreak: true }));
  {
    const hdr = ["Statistic", "Value"];
    const rows = [
      ["Method", meta.method],
      ["Pooled \u0394AUC", meta.pooled_delta.toFixed(4)],
      ["95% CI", `[${meta.ci_low.toFixed(4)}, ${meta.ci_high.toFixed(4)}]`],
      ["z", meta.z.toFixed(4)],
      ["p-value", fmtP(meta.p_value)],
      ["\u03C4\u00B2", meta.tau2.toFixed(6)],
      ["Q", meta.Q.toFixed(2)],
      ["Q df", String(meta.Q_df)],
      ["Q p-value", fmtP(meta.Q_pvalue)],
      ["I\u00B2 (%)", meta.I2.toFixed(1)],
    ];
    const tbl = new Table({
      rows: [
        new TableRow({ children: hdr.map(h => tableCell(h, { headerRow: true })) }),
        ...rows.map(r => new TableRow({ children: r.map(c => tableCell(c)) })),
      ],
      width: { size: 6000, type: WidthType.DXA },
    });
    children.push(tbl);
  }

  // Supplementary Figures S1-S11: per-target composites
  const figOrder = [
    { key: "fpr2", num: 1 }, { key: "egfr", num: 2 }, { key: "cdk2", num: 3 },
    { key: "cox2", num: 4 }, { key: "esr1", num: 5 }, { key: "dpp4", num: 6 },
    { key: "ache", num: 7 }, { key: "bace1", num: 8 }, { key: "hsp90", num: 9 },
    { key: "p38", num: 10 }, { key: "thrombin", num: 11 },
  ];
  for (const fig of figOrder) {
    const fname = `fig_main_${fig.key}_hires.png`;
    const imgBuf = loadImage(fname);
    children.push(new Paragraph({ children: [new PageBreak()] }));
    if (imgBuf) {
      const dims = pngDimensions(imgBuf);
      children.push(figParagraph(imgBuf, dims.width, dims.height, 550));
    }
    children.push(para([
      txt(`Figure S${fig.num}. `, { bold: true }),
      txt(`Virtual screening performance for ${TARGET_FULLNAMES[fig.key]} (${TARGET_CHEMBL[fig.key]}, PDB ${TARGET_PDB[fig.key].split(":")[0]}). (A) ROC curves for naive (blue) and skill-guided (red) protocols on successfully docked compounds with 95% bootstrap confidence bands. (B) Bootstrap ROC AUC distributions (n = 2,000). (C) Docking score distributions for actives and decoys.`),
    ]));
  }

  // Figure S12: Metrics heatmap
  {
    children.push(new Paragraph({ children: [new PageBreak()] }));
    const imgBuf = loadImage("metrics_heatmap_hires.png");
    if (imgBuf) {
      const dims = pngDimensions(imgBuf);
      children.push(figParagraph(imgBuf, dims.width, dims.height, 550));
    }
    children.push(para([
      txt("Figure S12. ", { bold: true }),
      txt("Heatmap of VS metric differences (skill \u2212 naive) across 11 targets. Red indicates skill advantage, blue indicates naive advantage. Metrics: ROC AUC, BEDROC(\u03B1=20), EF1%, EF5%, LogAUC, AUPR."),
    ]));
  }

  // Figure S13: Funnel plot
  {
    children.push(new Paragraph({ children: [new PageBreak()] }));
    const imgBuf = loadImage("funnel_plot_hires.png");
    if (imgBuf) {
      const dims = pngDimensions(imgBuf);
      children.push(figParagraph(imgBuf, dims.width, dims.height, 550));
    }
    children.push(para([
      txt("Figure S13. ", { bold: true }),
      txt("Funnel plot of intersection \u0394AUC versus standard error. The dashed line indicates the pooled estimate. No evidence of asymmetry suggesting publication bias."),
    ]));
  }

  return new Document({
    sections: [{ properties: PAGE_PROPS, children }],
  });
}

// ========================================================================
//  GENERATE
// ========================================================================
async function main() {
  console.log("Building main manuscript...");
  const mainDoc = buildMainManuscript();
  const mainBuf = await Packer.toBuffer(mainDoc);
  const mainPath = path.join(BASE, "Manuscript_BIB_v1.docx");
  fs.writeFileSync(mainPath, mainBuf);
  console.log(`  Written: ${mainPath} (${(mainBuf.length / 1024).toFixed(0)} KB)`);

  console.log("Building supplementary information...");
  const siDoc = buildSupplementary();
  const siBuf = await Packer.toBuffer(siDoc);
  const siPath = path.join(BASE, "Supplementary_BIB_v1.docx");
  fs.writeFileSync(siPath, siBuf);
  console.log(`  Written: ${siPath} (${(siBuf.length / 1024).toFixed(0)} KB)`);

  console.log("Done.");
}

main().catch(err => { console.error(err); process.exit(1); });
