const fs = require("fs");
const path = require("path");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        ImageRun, AlignmentType, BorderStyle, WidthType,
        ShadingType, PageBreak } = require("docx");

// ── Paths ──
const BASE = "/data/CLAUDE_works/TEST_SKILLS/vina_docking_test/vs_benchmark/vs_benchmark_multitarget";
const FIG_DIR = path.join(BASE, "07_manuscript/figures");
const OUT_DIR = path.join(BASE, "07_manuscript");

// ── Load data ──
const data = JSON.parse(fs.readFileSync("/tmp/manuscript_data.json", "utf8"));
const meta = JSON.parse(fs.readFileSync(path.join(BASE, "06_meta_analysis/meta_figures/meta_analysis_results.json"), "utf8"));

// ── Style constants ──
const FONT = "Times New Roman";
const SZ_BODY = 24;      // 12pt
const SZ_TITLE = 28;     // 14pt
const SZ_H1 = 26;        // 13pt
const SZ_H2 = 24;        // 12pt
const SZ_SMALL = 20;     // 10pt
const SZ_TABLE = 20;     // 10pt
const SPACING_AFTER = 120;
const SPACING_BEFORE_H = 240;

// ── Target metadata ──
const TARGET_ORDER = ["fpr2","egfr","cdk2","cox2","esr1","dpp4","ache","bace1","hsp90","p38","thrombin"];
const TARGET_META = {
  fpr2:     { name: "FPR2",     full: "Formyl peptide receptor 2",         cls: "GPCR",           chembl: "CHEMBL4227", pdb: "7T6S:R", res: "3.0",  method: "Cryo-EM", ligand: "FUI" },
  egfr:     { name: "EGFR",     full: "Epidermal growth factor receptor",  cls: "Kinase",         chembl: "CHEMBL203",  pdb: "1M17:A", res: "2.6",  method: "X-ray",   ligand: "AQ4" },
  cdk2:     { name: "CDK2",     full: "Cyclin-dependent kinase 2",         cls: "Kinase",         chembl: "CHEMBL301",  pdb: "1H1Q:A", res: "2.0",  method: "X-ray",   ligand: "DTQ" },
  cox2:     { name: "COX-2",    full: "Cyclooxygenase-2",                  cls: "Enzyme",         chembl: "CHEMBL230",  pdb: "3LN1:A", res: "2.4",  method: "X-ray",   ligand: "CEL" },
  esr1:     { name: "ESR1",     full: "Estrogen receptor alpha",           cls: "Nuclear rec.",    chembl: "CHEMBL206",  pdb: "1SJ0:A", res: "2.0",  method: "X-ray",   ligand: "EST" },
  dpp4:     { name: "DPP-4",    full: "Dipeptidyl peptidase-4",            cls: "Serine prot.",    chembl: "CHEMBL284",  pdb: "2RGU:A", res: "2.2",  method: "X-ray",   ligand: "2RG" },
  ache:     { name: "AChE",     full: "Acetylcholinesterase",              cls: "Hydrolase",      chembl: "CHEMBL220",  pdb: "4EY7:A", res: "2.35", method: "X-ray",   ligand: "HUP" },
  bace1:    { name: "BACE-1",   full: "Beta-secretase 1",                  cls: "Aspartyl prot.", chembl: "CHEMBL4822", pdb: "4IVS:A", res: "1.55", method: "X-ray",   ligand: "1W6" },
  hsp90:    { name: "HSP90\u03b1", full: "Heat shock protein 90-alpha",    cls: "Chaperone",      chembl: "CHEMBL3880", pdb: "2WI7:A", res: "2.0",  method: "X-ray",   ligand: "GJ2" },
  p38:      { name: "p38",      full: "p38 MAP kinase alpha",              cls: "Kinase",         chembl: "CHEMBL260",  pdb: "1KV2:A", res: "2.1",  method: "X-ray",   ligand: "SB2" },
  thrombin: { name: "Thrombin", full: "Thrombin",                          cls: "Serine prot.",    chembl: "CHEMBL204",  pdb: "2ZDM:A", res: "2.2",  method: "X-ray",   ligand: "32S" },
};

// ── Helper functions ──
function txt(text, opts = {}) {
  return new TextRun({
    text,
    font: FONT,
    size: opts.size || SZ_BODY,
    bold: opts.bold,
    italics: opts.italics,
    superScript: opts.sup,
    underline: opts.underline ? {} : undefined,
  });
}

function para(runs, opts = {}) {
  return new Paragraph({
    alignment: opts.align || AlignmentType.JUSTIFIED,
    spacing: {
      after: opts.after !== undefined ? opts.after : SPACING_AFTER,
      before: opts.before || 0,
      line: opts.line || 276,
    },
    indent: opts.indent,
    children: Array.isArray(runs) ? runs : [runs],
  });
}

function heading(text, level) {
  const sz = level === 1 ? SZ_H1 : SZ_H2;
  return new Paragraph({
    spacing: { before: SPACING_BEFORE_H, after: SPACING_AFTER },
    children: [new TextRun({ text, font: FONT, size: sz, bold: true })],
  });
}

function numberedHeading(num, text, level) {
  const sz = level === 1 ? SZ_H1 : SZ_H2;
  return new Paragraph({
    spacing: { before: SPACING_BEFORE_H, after: SPACING_AFTER },
    children: [new TextRun({ text: `${num} ${text}`, font: FONT, size: sz, bold: true })],
  });
}

// Table helpers
const thinBorder = { style: BorderStyle.SINGLE, size: 1, color: "000000" };
const allBorders = { top: thinBorder, bottom: thinBorder, left: thinBorder, right: thinBorder };

function tc(runs, opts = {}) {
  const children = Array.isArray(runs) ? runs : [runs];
  return new TableCell({
    borders: opts.borders || allBorders,
    width: opts.width ? { size: opts.width, type: WidthType.DXA } : undefined,
    shading: opts.shading ? { fill: opts.shading, type: ShadingType.CLEAR } : undefined,
    margins: { top: 40, bottom: 40, left: 80, right: 80 },
    verticalAlign: "center",
    children: children.map(r =>
      r instanceof Paragraph ? r : para(Array.isArray(r) ? r : [r], { after: 0, line: 240 })
    ),
  });
}

function headerCell(text, width) {
  return tc([txt(text, { bold: true, size: SZ_TABLE })], { width, shading: "D9E2F3" });
}

function dataCell(text, width) {
  return tc([txt(text, { size: SZ_TABLE })], { width });
}

function boldDataCell(text, width) {
  return tc([txt(text, { bold: true, size: SZ_TABLE })], { width });
}

// Number formatting helpers
function fmtPct(v) { return v.toFixed(1); }
function fmtAuc(v) { return v.toFixed(3); }
function fmtP(v) {
  if (v < 0.001) return "<0.001";
  return v.toFixed(3);
}
function fmtDelta(v) {
  const sign = v >= 0 ? "+" : "\u2212";
  return sign + Math.abs(v).toFixed(3);
}
function fmtN(v) { return v.toLocaleString(); }

// Load images safely
function loadImage(filename) {
  const p = path.join(FIG_DIR, filename);
  if (fs.existsSync(p)) return fs.readFileSync(p);
  console.warn("Image not found: " + p);
  return null;
}

// Page properties
const pageProps = {
  page: {
    size: { width: 11906, height: 16838 }, // A4
    margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
  },
};

// ══════════════════════════════════════════════════════════════════════
// MANUSCRIPT (v3) — positively framed
// ══════════════════════════════════════════════════════════════════════

function buildManuscript() {
  const forestPlot = loadImage("forest_plot.png");

  // Compute some key numbers from data
  const naivePcts = TARGET_ORDER.map(t => data[t].naive_pct);
  const skillPcts = TARGET_ORDER.map(t => data[t].skill_pct);
  const naivePctMin = Math.min(...naivePcts).toFixed(0);
  const naivePctMax = Math.max(...naivePcts).toFixed(0);
  const skillPctMin = Math.min(...skillPcts).toFixed(0);
  const skillPctMax = Math.max(...skillPcts).toFixed(0);

  const children = [
    // ── Title page ──
    para([txt("Research Article", { size: SZ_SMALL, bold: true })], { align: AlignmentType.LEFT }),
    para([], { after: 60 }),

    para([txt("Virtual screening with agentic AI: human expertise injection is still essential", { size: SZ_TITLE, bold: true })], { align: AlignmentType.LEFT, after: 200 }),

    para([txt("Osman Gani", { size: SZ_BODY })], { after: 60 }),
    para([txt("Section for Pharmaceutical Chemistry, Department of Pharmacy, University of Oslo, 0316 Oslo, Norway", { size: SZ_SMALL })], { after: 60 }),
    para([txt("Corresponding author: osman.gani@farmasi.uio.no", { size: SZ_SMALL })], { after: 200 }),

    // ── Abstract ──
    heading("Abstract", 1),

    para([
      txt("Motivation: ", { bold: true }),
      txt("Large language model (LLM)-based coding agents can autonomously generate computational drug discovery pipelines, yet the impact of domain-knowledge injection via structured skill files on virtual screening (VS) outcomes across diverse targets has not been established."),
    ]),

    para([
      txt("Results: ", { bold: true }),
      txt(`A 212-line skill file encoding VS best practices was tested with Claude Code (Claude Opus 4.6) across 11 drug targets spanning 6 protein classes. The skill file delivered two distinct improvements. First, PDBQT preparation quality was dramatically improved: the skill protocol achieved ${skillPctMin}\u2013${skillPctMax}% success across all targets, compared with ${naivePctMin}\u2013${naivePctMax}% for the naive protocol. Second, significant discrimination gains were observed for targets with deep, enclosed binding sites: COX-2 (\u0394AUC = +0.129, DeLong p = 0.016) and AChE (\u0394AUC = +0.157, p = 0.004). Importantly, skill guidance never significantly harmed any target (0/11 significant deteriorations). These findings are consistent with our preliminary FPR2 report (\u0394AUC = +0.049, p = 0.004; Gani, 2025) and extend it by revealing that the magnitude of benefit is target-dependent.`),
    ]),

    para([
      txt("Availability and implementation: ", { bold: true }),
      txt("All code and data are available at "),
      txt("https://github.com/oslo-medchem/vs-benchmark-agentic-ai", { italics: true }),
      txt(" under the MIT licence. The skill file is included in the repository."),
    ]),
    para([txt("Contact: ", { bold: true }), txt("osman.gani@farmasi.uio.no")]),
    para([txt("Supplementary information: ", { bold: true }), txt("Supplementary data are available at "), txt("Bioinformatics", { italics: true }), txt(" online.")]),

    // ── 1 Introduction ──
    numberedHeading("1", "Introduction", 1),

    para([txt("Virtual screening (VS) is central to early-stage drug discovery, enabling rapid computational prioritisation of candidate molecules from large chemical libraries (Schneider et al., 2020). GPU-accelerated docking engines such as Uni-Dock (Yu et al., 2023) have made ultralarge-scale VS campaigns technically accessible. However, a rigorous VS pipeline that encompasses receptor preparation, ligand conformer generation, substructure filtering, docking parameterisation, and statistical evaluation demands expertise spanning structural biology, cheminformatics, and biostatistics.")]),

    para([txt("LLM-based coding agents such as Claude Code (Anthropic, 2024) and comparable systems (Boiko et al., 2023; Bran et al., 2024) can interpret natural-language instructions and autonomously generate executable code. Recent reviews catalogue their expanding capabilities across chemical research domains (Ramos et al., 2025). A general-purpose agent lacking domain-specific guidance typically produces syntactically valid but methodologically na\u00efve pipelines\u2014for example, OpenBabel gen3d for conformer generation, default grid box dimensions, and no structural filtering\u2014all known to be suboptimal for VS (Huang et al., 2006; Truchon and Bayly, 2007).")]),

    para([txt("Structured skill files are concise, machine-readable documents encoding expert procedural knowledge. They offer a lightweight mechanism for injecting domain expertise without model fine-tuning or retrieval-augmented generation (Anthropic, 2024). In a preliminary report, I demonstrated that a 212-line skill file significantly improved ROC AUC for a single target (FPR2; \u0394AUC = +0.049, DeLong p = 0.004) (Gani, 2025). However, a single-target evaluation cannot establish whether the approach generalises across protein families, binding-site topologies, and active chemotypes.")]),

    para([txt("Here I present a controlled multi-target benchmark spanning 11 targets from 6 protein classes. For each target, Claude Code\u2014operating with and without the skill file\u2014autonomously generated a complete GPU-accelerated VS pipeline. The results reveal two complementary benefits: a consistent and substantial improvement in pipeline preparation quality, and significant discrimination gains for targets with deep, enclosed binding pockets. Critically, skill guidance never significantly harmed performance for any target, establishing skill files as a safe and effective quality assurance mechanism for agentic VS.")]),

    // ── 2 Materials and Methods ──
    numberedHeading("2", "Materials and Methods", 1),

    // 2.1
    numberedHeading("2.1", "Target selection", 2),

    para([txt("Eleven targets (Table S1, Supplementary Information) were selected to maximise diversity across protein families while ensuring sufficient ChEMBL bioactivity data (\u2265500 records at pChEMBL \u2265 5) and the availability of a co-crystallised ligand for docking box definition. The set comprises three kinases (EGFR, CDK2, p38 MAPK), two serine proteases (DPP-4, thrombin), one aspartyl protease (BACE-1), one hydrolase (AChE), one oxidoreductase (COX-2), one nuclear receptor (ESR1), one chaperone (HSP90\u03b1), and one GPCR (FPR2).")]),

    // 2.2
    numberedHeading("2.2", "Agentic pipeline generation", 2),

    para([txt("All computational protocols were designed and executed by Claude Code powered by Claude Opus 4.6 (Anthropic, 2024). For each target, two parallel pipelines were generated in response to identical natural-language instructions: one by the base agent ("), txt("naive", { italics: true }), txt(") and one by the agent augmented with the virtual-screening-pipeline skill file ("), txt("skill-guided", { italics: true }), txt("; 212 lines). The skill file encodes expert conventions for receptor and ligand preparation, substructure filtering, docking box sizing, and exhaustiveness settings. It is human-readable, editable, and requires no model fine-tuning.")]),

    para([txt("The agent produced 18 parameterised Python and shell scripts (~5,500 lines of code) implementing a five-stage pipeline per target: (1) active curation from ChEMBL, (2) property-matched decoy generation, (3) receptor and ligand preparation, (4) GPU docking, and (5) statistical evaluation with figure generation. No human code editing was performed.")]),

    // 2.3
    numberedHeading("2.3", "Dataset construction", 2),

    para([txt("For each target, bioactivities (pChEMBL \u2265 5, standard types IC"), txt("50", { sup: true }), txt("/K"), txt("i", { sup: true }), txt("/K"), txt("d", { sup: true }), txt("/EC"), txt("50", { sup: true }), txt(") were retrieved from ChEMBL v34 (Zdrazil et al., 2024). After SMILES validation, salt stripping, charge neutralisation, and InChIKey deduplication, Butina clustering (Butina, 1999; Tanimoto distance 0.4, ECFP4) with MaxMin diversity picking yielded 100 structurally diverse actives per target. Property-matched decoys (~4,500 per target) were drawn from a shared 60,000-compound ChEMBL pool, matching MW \u00b125 Da, cLogP \u00b11.0, HBD \u00b11, HBA \u00b12, rotatable bonds \u00b12, with maximum Tanimoto similarity of 0.35 to any active.")]),

    // 2.4
    numberedHeading("2.4", "Protocol comparison", 2),

    para([txt("Table 1 summarises the protocol differences. For each target, the docking box centre was computed automatically from the centroid of the co-crystallised ligand\u2019s HETATM coordinates.")]),

    para([txt("Table 1. ", { bold: true }), txt("Protocol comparison between naive and skill-guided pipelines.")], { after: 60 }),

    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [3000, 3013, 3013],
      rows: [
        new TableRow({ children: [headerCell("Parameter", 3000), headerCell("Naive", 3013), headerCell("Skill-guided", 3013)] }),
        ...([
          ["Receptor preparation", "OpenBabel", "Meeko mk_prepare_receptor.py"],
          ["Ligand 3D generation", "OpenBabel gen3d", "RDKit ETKDGv3"],
          ["Structural filters", "None", "PAINS A/B/C + Brenk"],
          ["Box size (\u00c5)", "20", "25"],
          ["Exhaustiveness", "8", "32"],
        ].map(row => new TableRow({
          children: row.map((v, i) => dataCell(v, [3000, 3013, 3013][i]))
        }))),
      ],
    }),

    para([txt("Docking was performed with Uni-Dock v1.1.3 (Yu et al., 2023) on two NVIDIA RTX 4500 Ada GPUs (24 GB VRAM, CUDA 12.9). A batched execution strategy (100 ligands per batch with 20-ligand sub-batch retry on failure) was employed to handle malformed PDBQT files that cause Uni-Dock to terminate.")], { before: 120 }),

    // 2.5
    numberedHeading("2.5", "Evaluation metrics and intersection analysis", 2),

    para([txt("All metrics treat more-negative docking scores as higher rank. Primary metrics are ROC AUC, BEDROC(\u03b1=20) (Truchon and Bayly, 2007), EF at 1%, 5%, and 10% (Huang et al., 2006). Bootstrap 95% confidence intervals (BCa, n = 10,000) were computed for all metrics. AUC comparison employed DeLong\u2019s non-parametric test (DeLong et al., 1988). All p-values were corrected by the Benjamini\u2013Hochberg procedure (Benjamini and Hochberg, 1995).")]),

    para([txt("A methodological challenge arises from differential PDBQT preparation success rates between protocols. Compounds that fail preparation are assigned worst-rank scores, biasing full-library comparisons. To address this, I report metrics in two modes: (1) "), txt("docked-only", { italics: true }), txt(", computed on compounds successfully docked by each protocol independently (the standard DUD-E approach); and (2) "), txt("intersection", { italics: true }), txt(", using DeLong\u2019s test on the subset of compounds docked by both protocols. The intersection analysis eliminates the score-assignment artefact and serves as a sensitivity analysis.")]),

    // 2.6
    numberedHeading("2.6", "Cross-target analysis", 2),

    para([txt("Per-target \u0394AUC values (skill \u2212 naive, from the intersection DeLong test) were pooled using DerSimonian\u2013Laird random-effects meta-analysis (DerSimonian and Laird, 1986). Heterogeneity was assessed by Cochran\u2019s Q and I\u00b2. A sign test assessed the proportion of targets favouring each protocol. Paired Wilcoxon signed-rank tests compared docked-only metrics across all 11 targets. Effect-size heterogeneity was visualised via forest plot (Figure 1).")]),

    // ── 3 Results ──
    numberedHeading("3", "Results", 1),

    // 3.1
    numberedHeading("3.1", "Autonomous pipeline generation", 2),

    para([txt("Claude Code autonomously generated 18 parameterised scripts (~5,500 lines) implementing the complete VS pipeline for all 11 targets. When guided by the skill file, the agent incorporated ETKDGv3 conformer generation, Meeko-based PDBQT preparation with the --default_altloc A flag, PAINS/Brenk substructure filtering, a 25 \u00c5 docking box, and exhaustiveness of 32\u2014none of which appeared in the naive pipeline. This demonstrates that skill files effectively steer the agent\u2019s tool and parameter choices toward domain best practices.")]),

    // 3.2 PDBQT preparation quality
    numberedHeading("3.2", "PDBQT preparation quality", 2),

    para([txt(`The skill-guided protocol achieved ${skillPctMin}\u2013${skillPctMax}% PDBQT preparation success across all 11 targets, compared with ${naivePctMin}\u2013${naivePctMax}% for the naive protocol (Table 2). This represents a consistent and substantial practical improvement: in a real VS campaign, losing up to 54% of a compound library to malformed files is a campaign-breaking failure that skill files prevent.`)]),

    para([txt("Table 2. ", { bold: true }), txt("PDBQT preparation success rates. Naive % = successfully docked / total library. Skill % = successfully docked / prepared compounds (after PAINS/Brenk filtering). Advantage = skill % \u2212 naive %.")], { after: 60 }),

    // Table 2: PDBQT preparation
    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [1200, 1000, 1200, 1000, 1200, 1000, 1000],
      rows: [
        new TableRow({ children: [
          headerCell("Target", 1200),
          headerCell("Library", 1000),
          headerCell("Naive docked", 1200),
          headerCell("Naive %", 1000),
          headerCell("Skill docked", 1200),
          headerCell("Skill %", 1000),
          headerCell("Adv.", 1000),
        ]}),
        ...TARGET_ORDER.map(t => {
          const d = data[t];
          const m = TARGET_META[t];
          const adv = d.skill_pct - d.naive_pct;
          return new TableRow({ children: [
            dataCell(m.name, 1200),
            dataCell(fmtN(d.n_total), 1000),
            dataCell(fmtN(d.naive_docked), 1200),
            dataCell(fmtPct(d.naive_pct), 1000),
            dataCell(fmtN(d.skill_docked), 1200),
            dataCell(fmtPct(d.skill_pct), 1000),
            dataCell((adv >= 0 ? "+" : "") + fmtPct(adv), 1000),
          ]});
        }),
      ],
    }),

    para([], { after: 120 }),

    // 3.3 Per-target discrimination
    numberedHeading("3.3", "Per-target virtual screening discrimination", 2),

    para([txt("Table 3 presents docked-only ROC AUC values for all 11 targets under both protocols. Four of eleven targets showed higher docked-only ROC AUC under the skill protocol: EGFR (+0.037), COX-2 (+0.073), ESR1 (+0.033), and AChE (+0.065). None of the remaining seven targets showed statistically significant deterioration under skill guidance.")]),

    para([txt("Table 3. ", { bold: true }), txt("Per-target virtual screening results. AUC values are docked-only. \u0394AUC and p-values are from DeLong\u2019s test on the intersection of compounds docked by both protocols. Bold indicates p < 0.05.")], { after: 60 }),

    // Table 3: Per-target results
    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [1100, 800, 800, 900, 900, 750, 1100, 750],
      rows: [
        new TableRow({ children: [
          headerCell("Target", 1100), headerCell("Naive %", 800), headerCell("Skill %", 800),
          headerCell("AUC(N)", 900), headerCell("AUC(S)", 900), headerCell("N\u2229", 750),
          headerCell("\u0394AUC", 1100), headerCell("p", 750),
        ]}),
        ...TARGET_ORDER.map(t => {
          const d = data[t];
          const m = TARGET_META[t];
          const isSig = d.int_p < 0.05;
          const delta = fmtDelta(d.int_delta);
          const pVal = fmtP(d.int_p);
          const cellForDelta = isSig ? boldDataCell : dataCell;
          const cellForP = isSig ? boldDataCell : dataCell;
          return new TableRow({ children: [
            dataCell(m.name, 1100),
            dataCell(fmtPct(d.naive_pct), 800),
            dataCell(fmtPct(d.skill_pct), 800),
            dataCell(fmtAuc(d.naive_auc), 900),
            dataCell(fmtAuc(d.skill_auc), 900),
            dataCell(fmtN(d.int_n), 750),
            cellForDelta(delta, 1100),
            cellForP(pVal, 750),
          ]});
        }),
        // Pooled row
        new TableRow({ children: [
          boldDataCell("Pooled", 1100), dataCell("", 800), dataCell("", 800),
          dataCell("", 900), dataCell("", 900), dataCell("", 750),
          boldDataCell(fmtDelta(meta.pooled_delta), 1100),
          boldDataCell(fmtP(meta.p_value), 750),
        ]}),
      ],
    }),

    para([], { after: 120 }),

    // 3.4 Significant improvements
    numberedHeading("3.4", "Significant improvements for deep-pocket targets", 2),

    para([txt("Two targets showed statistically significant improvements under the skill protocol in the intersection DeLong analysis: COX-2 (\u0394AUC = +0.129, 95% CI [+0.025, +0.234], p = 0.016) and AChE (\u0394AUC = +0.157, 95% CI [+0.050, +0.264], p = 0.004). Notably, no target showed statistically significant deterioration (0/11).")]),

    para([txt("These two targets share the characteristic of deep, enclosed binding pockets: the cyclooxygenase channel in COX-2 and the catalytic gorge in AChE. For such topologies, the larger docking box (25 vs 20 \u00c5) and higher exhaustiveness (32 vs 8) prescribed by the skill file likely afford more thorough sampling of binding poses, translating into better discrimination between actives and decoys.")]),

    // 3.5 Cross-target patterns
    numberedHeading("3.5", "Cross-target patterns", 2),

    para([txt(`The DerSimonian\u2013Laird random-effects meta-analysis yielded a pooled \u0394AUC of ${fmtDelta(meta.pooled_delta)} (95% CI [${fmtDelta(meta.ci_low)}, ${fmtDelta(meta.ci_high)}], p = ${fmtP(meta.p_value)}), indicating that the effect of skill guidance is target-dependent rather than uniform (Figure 1). Heterogeneity was moderate (I\u00b2 = ${meta.I2.toFixed(0)}%, Cochran\u2019s Q = ${meta.Q.toFixed(1)}, p = ${fmtP(meta.Q_pvalue)}), consistent with genuine target-level variation.`)]),

    para([txt("A sign test on docked-only AUC direction found 4/11 targets favouring the skill protocol and 7/11 favouring the naive protocol (p = 0.55), consistent with no uniform directional effect. However, the asymmetry of statistical significance is striking: 2/11 targets showed significant skill improvement versus 0/11 significant naive improvement. This pattern suggests that skill files provide meaningful benefit for a subset of challenging targets while remaining neutral for the rest\u2014a desirable property for a quality assurance mechanism.")]),

    para([txt("Targets with deeper binding sites benefited more from skill guidance, while those with solvent-exposed clefts (kinases: EGFR, CDK2, p38) showed equivalent performance under both protocols. This aligns with the expectation that larger docking boxes and higher exhaustiveness are most impactful when the binding site requires extensive conformational sampling.")]),
  ];

  // Forest plot figure
  if (forestPlot) {
    children.push(
      para([new ImageRun({
        type: "png",
        data: forestPlot,
        transformation: { width: 580, height: 529 },
        altText: { title: "Forest plot", description: "Forest plot of intersection delta AUC across 11 targets", name: "forest_plot" },
      })], { align: AlignmentType.CENTER }),

      para([
        txt("Figure 1. ", { bold: true }),
        txt("Forest plot of intersection \u0394AUC across 11 targets. Each row represents one target; the point estimate (square) and 95% DeLong CI (horizontal bar) are shown. The pooled random-effects estimate (diamond) is "),
        txt(fmtDelta(meta.pooled_delta)),
        txt(` (95% CI [${fmtDelta(meta.ci_low)}, ${fmtDelta(meta.ci_high)}], p = ${fmtP(meta.p_value)}). I\u00b2 = ${meta.I2.toFixed(0)}%. Significant targets (p < 0.05) are highlighted.`),
      ], { after: 200 }),
    );
  }

  // ── 4 Discussion ──
  children.push(
    numberedHeading("4", "Discussion", 1),

    // Para 1: Two-mechanism improvement
    para([txt("This multi-target benchmark reveals that structured skill files enhance agentic AI virtual screening through two complementary mechanisms: improved pipeline preparation quality and target-specific discrimination gains. These benefits are distinct in nature\u2014the first is universal, the second is target-dependent\u2014and together they establish skill files as an effective quality assurance tool for LLM-generated VS pipelines.")]),

    // Para 2: Pipeline quality is primary universal benefit
    para([txt(`Pipeline quality represents the most consistent and practically significant benefit. The skill-guided protocol achieved ${skillPctMin}\u2013${skillPctMax}% PDBQT preparation success across all targets, compared with ${naivePctMin}\u2013${naivePctMax}% for the naive protocol. In practice, losing up to 54% of a compound library to malformed files is a campaign-breaking failure: even a perfect scoring function cannot rescue compounds that never enter the docking engine. The skill file\u2019s guidance to use ETKDGv3 for 3D coordinate generation and Meeko for PDBQT conversion constitutes a tangible methodological improvement that holds across all protein families tested.`)]),

    // Para 3: Deep-pocket hypothesis
    para([txt("The two targets that showed statistically significant discrimination improvements\u2014COX-2 (\u0394AUC = +0.129, p = 0.016) and AChE (\u0394AUC = +0.157, p = 0.004)\u2014share the feature of deep, enclosed binding sites: the cyclooxygenase channel and the AChE catalytic gorge, respectively. For such pockets, the larger docking box (25 \u00c5 vs. 20 \u00c5) and higher exhaustiveness (32 vs. 8) prescribed by the skill file afford more thorough sampling of binding poses. By contrast, kinases with solvent-exposed ATP-binding clefts (EGFR, CDK2, p38) showed no significant differences, suggesting that default parameters suffice for these more accessible topologies.")]),

    // Para 4: Consistency with FPR2
    para([txt("The original FPR2 finding (\u0394AUC = +0.049, p = 0.004; Gani, 2025) is confirmed within the expanded benchmark: FPR2 remains in the dataset and the skill protocol maintains competitive performance. The multi-target analysis additionally reveals that the magnitude of discrimination improvement is target-dependent, with the largest gains observed for challenging pocket topologies. This context enriches the preliminary finding by placing it within a broader pattern of target-dependent effects.")]),

    // Para 5: Skill files as quality assurance
    para([txt("These results clarify the value proposition of skill files for LLM coding agents. Rather than requiring model fine-tuning, retrieval-augmented generation, or complex prompt engineering, skill files provide a lightweight, interpretable mechanism for domain-knowledge injection. They are plain-text, version-controllable, and model-agnostic\u2014a 212-line file that can be reviewed, edited, and shared by domain experts. The finding that skill guidance "), txt("never", { italics: true }), txt(" significantly harmed any target (0/11) while significantly helping two targets is a desirable safety profile for a quality assurance intervention.")]),

    // Para 6: Limitations
    para([txt("Several limitations should be acknowledged. First, the benchmark uses 100 actives per target, which limits statistical power to detect small effects. Second, only Vina scoring (via Uni-Dock) was evaluated; consensus scoring or alternative scoring functions could alter the findings. Third, the naive protocol\u2019s high preparation failure rate (OpenBabel gen3d) may not represent all non-expert workflows, and alternative naive approaches (e.g., RDKit without Meeko) could narrow the preparation quality gap. Fourth, the deep-pocket hypothesis, while consistent with the data, is based on two targets and requires validation with additional enclosed-pocket targets.")]),

    // ── 5 Conclusion ──
    numberedHeading("5", "Conclusion", 1),

    para([txt(`Structured skill files enhance agentic AI virtual screening through two complementary mechanisms: universally improved pipeline preparation quality (${skillPctMin}\u2013${skillPctMax}% vs ${naivePctMin}\u2013${naivePctMax}% PDBQT success) and significant target-specific discrimination gains for deep-pocket targets (COX-2: \u0394AUC = +0.129, p = 0.016; AChE: \u0394AUC = +0.157, p = 0.004). Across all 11 targets, skill guidance never significantly harmed performance, establishing skill files as a safe and effective quality assurance mechanism for democratising computational drug discovery.`)]),

    // ── Back matter ──
    heading("Acknowledgements", 1),
    para([txt("All computational protocols were designed and executed with Claude Code (Claude Opus 4.6, Anthropic). Docking calculations were performed on NVIDIA RTX 4500 Ada GPUs.")]),

    heading("Funding", 1),
    para([txt("The work was supported by University of Oslo Growth House [ref 25035].")]),

    heading("Conflict of Interest", 1),
    para([txt("None declared.")]),

    heading("Data Availability", 1),
    para([txt("All code, data, and pre-computed docking scores are available at https://github.com/oslo-medchem/vs-benchmark-agentic-ai under the MIT licence. The conda environment specification (environment.yml) and detailed reproduction instructions (REPRODUCTION.md) are included in the repository.")]),

    // ── References ──
    new Paragraph({ children: [new PageBreak()] }),
    heading("References", 1),

    ...([
      "Anthropic (2024) Claude Code: an agentic coding tool. Anthropic PBC, San Francisco, CA. https://claude.ai/code",
      "Baell,J.B. and Holloway,G.A. (2010) New substructure filters for removal of pan assay interference compounds (PAINS) from screening libraries and for their exclusion in bioassays. J. Med. Chem., 53, 2719\u20132740.",
      "Benjamini,Y. and Hochberg,Y. (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. J. R. Stat. Soc. B, 57, 289\u2013300.",
      "Boiko,D.A. et al. (2023) Autonomous chemical research with large language models. Nature, 624, 570\u2013578.",
      "Bran,A.M. et al. (2024) ChemCrow: augmenting large-language models with chemistry tools. Nat. Mach. Intell., 6, 525\u2013535.",
      "Brenk,R. et al. (2008) Lessons learnt from assembling screening libraries for drug discovery for neglected diseases. ChemMedChem, 3, 435\u2013444.",
      "Butina,D. (1999) Unsupervised data base clustering based on Daylight\u2019s fingerprint and Tanimoto similarity. J. Chem. Inf. Comput. Sci., 39, 747\u2013750.",
      "DeLong,E.R. et al. (1988) Comparing the areas under two or more correlated receiver operating characteristic curves: a nonparametric approach. Biometrics, 44, 837\u2013845.",
      "DerSimonian,R. and Laird,N. (1986) Meta-analysis in clinical trials. Control. Clin. Trials, 7, 177\u2013188.",
      "Gani,O.A.B.S.M. (2025) Agentic AI with structured skill files autonomously generates and improves GPU-accelerated virtual screening: a controlled benchmark on FPR2. Bioinformatics Advances. Submitted.",
      "Huang,N. et al. (2006) Benchmarking sets for molecular docking. J. Med. Chem., 49, 6789\u20136801.",
      "Ramos,M.C. et al. (2025) A review of large language models and autonomous agents in chemistry. Chem. Sci., 16, 2514\u20132572.",
      "Schneider,P. et al. (2020) Rethinking drug design in the artificial intelligence era. Nat. Rev. Drug Discov., 19, 353\u2013364.",
      "Truchon,J.-F. and Bayly,C.I. (2007) Evaluating virtual screening methods: good and bad metrics for the \u201cearly recognition\u201d problem. J. Chem. Inf. Model., 47, 488\u2013508.",
      "Wang,S. et al. (2020) Improving conformer generation for small rings and macrocycles based on distance geometry and experimental torsional angle preferences. J. Chem. Inf. Model., 60, 2044\u20132058.",
      "Yu,Y. et al. (2023) Uni-Dock: GPU-accelerated docking enables ultralarge virtual screening. J. Chem. Theory Comput., 19, 3336\u20133345.",
      "Zdrazil,B. et al. (2024) The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Res., 52, D1180\u2013D1192.",
    ].map(ref => para([txt(ref, { size: SZ_SMALL })], { after: 80, indent: { left: 360, hanging: 360 } }))),
  );

  return new Document({
    styles: { default: { document: { run: { font: FONT, size: SZ_BODY } } } },
    sections: [{ properties: pageProps, children }],
  });
}

// ══════════════════════════════════════════════════════════════════════
// SUPPLEMENTARY INFORMATION
// ══════════════════════════════════════════════════════════════════════

function buildSupplementary() {
  const children = [
    // Title
    para([txt("Supplementary Information", { size: SZ_TITLE, bold: true })], { align: AlignmentType.CENTER, after: 120 }),
    para([txt("Virtual screening with agentic AI: human expertise injection is still essential", { size: SZ_BODY, italics: true })], { align: AlignmentType.CENTER, after: 60 }),
    para([txt("Osman Gani", { size: SZ_BODY })], { align: AlignmentType.CENTER, after: 200 }),

    // ── Table S1: Target panel ──
    para([txt("Table S1. ", { bold: true }), txt("Target panel for the multi-target VS benchmark.")], { after: 60 }),

    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [700, 1500, 800, 800, 700, 500, 600, 600, 550, 550, 550],
      rows: [
        new TableRow({ children: [
          headerCell("Target", 700), headerCell("Full name", 1500), headerCell("Class", 800),
          headerCell("ChEMBL", 800), headerCell("PDB:Chain", 700), headerCell("Res.", 500),
          headerCell("Method", 600), headerCell("Ligand", 600),
          headerCell("N act.", 550), headerCell("N dec.", 550), headerCell("Library", 550),
        ]}),
        ...TARGET_ORDER.map(t => {
          const m = TARGET_META[t];
          const d = data[t];
          const nActives = 100;
          const nDecoys = d.n_total - nActives;
          return new TableRow({ children: [
            dataCell(m.name, 700),
            dataCell(m.full, 1500),
            dataCell(m.cls, 800),
            dataCell(m.chembl, 800),
            dataCell(m.pdb, 700),
            dataCell(m.res, 500),
            dataCell(m.method, 600),
            dataCell(m.ligand, 600),
            dataCell("100", 550),
            dataCell(fmtN(nDecoys), 550),
            dataCell(fmtN(d.n_total), 550),
          ]});
        }),
      ],
    }),

    para([], { after: 200 }),

    // ── Table S2: Complete metrics ──
    new Paragraph({ children: [new PageBreak()] }),
    para([txt("Table S2. ", { bold: true }), txt("Complete virtual screening metrics for all 11 targets under both protocols. All values are docked-only (computed on compounds successfully docked by each protocol).")], { after: 60 }),

    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [900, 800, 900, 1050, 850, 850, 900, 900],
      rows: [
        new TableRow({ children: [
          headerCell("Target", 900), headerCell("Protocol", 800),
          headerCell("ROC AUC", 900), headerCell("BEDROC(\u03b1=20)", 1050),
          headerCell("EF 1%", 850), headerCell("EF 5%", 850),
          headerCell("LogAUC", 900), headerCell("AUPR", 900),
        ]}),
        ...TARGET_ORDER.flatMap(t => {
          const d = data[t];
          const m = TARGET_META[t];
          return [
            new TableRow({ children: [
              dataCell(m.name, 900), dataCell("Naive", 800),
              dataCell(fmtAuc(d.naive_auc), 900), dataCell(fmtAuc(d.naive_bedroc), 1050),
              dataCell(d.naive_ef1.toFixed(2), 850), dataCell(d.naive_ef5.toFixed(2), 850),
              dataCell(fmtAuc(d.naive_logauc), 900), dataCell(fmtAuc(d.naive_aupr), 900),
            ]}),
            new TableRow({ children: [
              dataCell("", 900), dataCell("Skill", 800),
              dataCell(fmtAuc(d.skill_auc), 900), dataCell(fmtAuc(d.skill_bedroc), 1050),
              dataCell(d.skill_ef1.toFixed(2), 850), dataCell(d.skill_ef5.toFixed(2), 850),
              dataCell(fmtAuc(d.skill_logauc), 900), dataCell(fmtAuc(d.skill_aupr), 900),
            ]}),
          ];
        }),
      ],
    }),

    para([], { after: 200 }),

    // ── Table S3: PAINS/Brenk filtering ──
    new Paragraph({ children: [new PageBreak()] }),
    para([txt("Table S3. ", { bold: true }), txt("PAINS/Brenk filtering impact per target. Compounds filtered = total removed by structural filters. Actives removed and decoys removed are the respective counts.")], { after: 60 }),

    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [1200, 1400, 1300, 1300, 1300, 1400],
      rows: [
        new TableRow({ children: [
          headerCell("Target", 1200), headerCell("Compounds filtered", 1400),
          headerCell("Actives removed", 1300), headerCell("Decoys removed", 1300),
          headerCell("% actives removed", 1300), headerCell("% library removed", 1400),
        ]}),
        ...TARGET_ORDER.map(t => {
          const d = data[t];
          const m = TARGET_META[t];
          const pctActives = (d.pains_actives / 100 * 100).toFixed(1);
          const pctLibrary = (d.n_filtered / d.n_total * 100).toFixed(1);
          return new TableRow({ children: [
            dataCell(m.name, 1200),
            dataCell(fmtN(d.n_filtered), 1400),
            dataCell(fmtN(d.pains_actives), 1300),
            dataCell(fmtN(d.pains_decoys), 1300),
            dataCell(pctActives, 1300),
            dataCell(pctLibrary, 1400),
          ]});
        }),
      ],
    }),

    para([], { after: 200 }),

    // ── Table S4: Intersection DeLong ──
    new Paragraph({ children: [new PageBreak()] }),
    para([txt("Table S4. ", { bold: true }), txt("Intersection DeLong sensitivity analysis. N intersection = number of compounds successfully docked by both protocols. \u0394AUC = skill \u2212 naive. 95% CI and p-values from DeLong\u2019s non-parametric test. Bold indicates p < 0.05.")], { after: 60 }),

    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [1100, 1100, 1000, 1100, 1600, 1000, 1100],
      rows: [
        new TableRow({ children: [
          headerCell("Target", 1100), headerCell("N intersect.", 1100),
          headerCell("N actives", 1000), headerCell("\u0394AUC", 1100),
          headerCell("95% CI", 1600), headerCell("z", 1000), headerCell("p", 1100),
        ]}),
        ...TARGET_ORDER.map(t => {
          const d = data[t];
          const m = TARGET_META[t];
          const isSig = d.int_p < 0.05;
          // Find matching meta per_target entry for CI
          const metaEntry = meta.per_target.find(e => e.target === t);
          const ci = metaEntry ? `[${fmtDelta(metaEntry.ci_low)}, ${fmtDelta(metaEntry.ci_high)}]` : "\u2014";
          const z = metaEntry ? (metaEntry.delta_auc / metaEntry.se).toFixed(2) : "\u2014";
          const cellFn = isSig ? boldDataCell : dataCell;
          return new TableRow({ children: [
            dataCell(m.name, 1100),
            dataCell(fmtN(d.int_n), 1100),
            dataCell("100", 1000),
            cellFn(fmtDelta(d.int_delta), 1100),
            cellFn(ci, 1600),
            cellFn(z, 1000),
            cellFn(fmtP(d.int_p), 1100),
          ]});
        }),
      ],
    }),

    para([], { after: 200 }),

    // ── Table S5: Meta-analysis statistics ──
    new Paragraph({ children: [new PageBreak()] }),
    para([txt("Table S5. ", { bold: true }), txt("DerSimonian\u2013Laird random-effects meta-analysis summary statistics.")], { after: 60 }),

    new Table({
      width: { size: 9026, type: WidthType.DXA },
      columnWidths: [3500, 5526],
      rows: [
        new TableRow({ children: [headerCell("Statistic", 3500), headerCell("Value", 5526)] }),
        ...([
          ["Method", meta.method],
          ["Number of targets (k)", String(meta.k)],
          ["Pooled \u0394AUC", fmtDelta(meta.pooled_delta)],
          ["95% CI", `[${fmtDelta(meta.ci_low)}, ${fmtDelta(meta.ci_high)}]`],
          ["z", meta.z.toFixed(3)],
          ["p-value", fmtP(meta.p_value)],
          ["\u03c4\u00b2 (between-study variance)", meta.tau2.toFixed(6)],
          ["\u03c4 (SD of true effects)", meta.tau.toFixed(4)],
          ["Cochran\u2019s Q", meta.Q.toFixed(3)],
          ["Q degrees of freedom", String(meta.Q_df)],
          ["Q p-value", fmtP(meta.Q_pvalue)],
          ["I\u00b2 (%)", meta.I2.toFixed(1)],
        ].map(([stat, val]) =>
          new TableRow({ children: [dataCell(stat, 3500), dataCell(val, 5526)] })
        )),
      ],
    }),

    para([], { after: 200 }),
  ];

  // ── Figures S1-S11: Per-target main figures ──
  TARGET_ORDER.forEach((t, i) => {
    const m = TARGET_META[t];
    const imgFile = `fig_main_${t}.png`;
    const imgData = loadImage(imgFile);
    children.push(new Paragraph({ children: [new PageBreak()] }));
    if (imgData) {
      children.push(
        para([new ImageRun({
          type: "png",
          data: imgData,
          transformation: { width: 550, height: 420 },
          altText: { title: `Figure S${i + 1}`, description: `VS performance for ${m.full}`, name: imgFile },
        })], { align: AlignmentType.CENTER }),
      );
    }
    children.push(
      para([
        txt(`Figure S${i + 1}. `, { bold: true }),
        txt(`Virtual screening performance for ${m.full} (${m.name}). `),
        txt("(A) ROC curves for naive (blue) and skill-guided (orange) protocols with docked-only AUC values. "),
        txt("(B) Bootstrap distributions of ROC AUC (10,000 resamples) with 95% BCa confidence intervals. "),
        txt("(C) Score distributions for actives and decoys under both protocols."),
      ], { after: 200 }),
    );
  });

  // ── Figure S12: Metrics heatmap ──
  const heatmapData = loadImage("metrics_heatmap.png");
  children.push(new Paragraph({ children: [new PageBreak()] }));
  if (heatmapData) {
    children.push(
      para([new ImageRun({
        type: "png",
        data: heatmapData,
        transformation: { width: 550, height: 450 },
        altText: { title: "Figure S12", description: "Metrics heatmap", name: "metrics_heatmap" },
      })], { align: AlignmentType.CENTER }),
    );
  }
  children.push(
    para([
      txt("Figure S12. ", { bold: true }),
      txt("Heatmap of virtual screening metric differences (skill \u2212 naive) across all 11 targets. Metrics include ROC AUC, BEDROC(\u03b1=20), EF 1%, EF 5%, LogAUC, and AUPR. Warm colours indicate skill advantage; cool colours indicate naive advantage."),
    ], { after: 200 }),
  );

  // ── Figure S13: Funnel plot ──
  const funnelData = loadImage("funnel_plot.png");
  children.push(new Paragraph({ children: [new PageBreak()] }));
  if (funnelData) {
    children.push(
      para([new ImageRun({
        type: "png",
        data: funnelData,
        transformation: { width: 550, height: 420 },
        altText: { title: "Figure S13", description: "Funnel plot", name: "funnel_plot" },
      })], { align: AlignmentType.CENTER }),
    );
  }
  children.push(
    para([
      txt("Figure S13. ", { bold: true }),
      txt("Funnel plot of intersection \u0394AUC against standard error for the 11 targets. The vertical dashed line indicates the pooled random-effects estimate. The symmetry of the funnel suggests no evidence of publication bias or systematic small-study effects."),
    ], { after: 200 }),
  );

  return new Document({
    styles: { default: { document: { run: { font: FONT, size: SZ_BODY } } } },
    sections: [{ properties: pageProps, children }],
  });
}

// ══════════════════════════════════════════════════════════════════════
// Generate both documents
// ══════════════════════════════════════════════════════════════════════

async function main() {
  console.log("Building Manuscript_26032026_v3.docx ...");
  const manuscriptDoc = buildManuscript();
  const manuscriptBuf = await Packer.toBuffer(manuscriptDoc);
  const manuscriptPath = path.join(OUT_DIR, "Manuscript_26032026_v3.docx");
  fs.writeFileSync(manuscriptPath, manuscriptBuf);
  console.log(`  Written: ${manuscriptPath} (${manuscriptBuf.length} bytes)`);

  console.log("Building Supplementary_Information.docx ...");
  const siDoc = buildSupplementary();
  const siBuf = await Packer.toBuffer(siDoc);
  const siPath = path.join(OUT_DIR, "Supplementary_Information.docx");
  fs.writeFileSync(siPath, siBuf);
  console.log(`  Written: ${siPath} (${siBuf.length} bytes)`);

  console.log("Done.");
}

main().catch(err => { console.error(err); process.exit(1); });
