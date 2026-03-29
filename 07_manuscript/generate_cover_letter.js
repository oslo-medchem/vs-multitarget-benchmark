const fs = require("fs");
const { Document, Packer, Paragraph, TextRun, AlignmentType } = require("docx");

const FONT = "Times New Roman";
const SZ = 24; // 12pt
const SZ_SM = 22; // 11pt

function txt(text, opts = {}) {
  return new TextRun({ text, font: FONT, size: opts.size || SZ, bold: opts.bold, italics: opts.italics });
}

function para(runs, opts = {}) {
  return new Paragraph({
    alignment: opts.align || AlignmentType.LEFT,
    spacing: { after: opts.after !== undefined ? opts.after : 200, line: 320 },
    children: Array.isArray(runs) ? runs : [runs],
  });
}

const doc = new Document({
  sections: [{
    properties: {
      page: {
        size: { width: 11906, height: 16838 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children: [
      // Date and address
      para([txt("27 March 2026")], { after: 100 }),
      para([], { after: 60 }),
      para([txt("The Editor")]),
      para([txt("Briefings in Bioinformatics", { italics: true })]),
      para([txt("Oxford University Press")], { after: 200 }),

      // Salutation
      para([txt("Dear Editor,")], { after: 200 }),

      // Subject line
      para([
        txt("Re: Submission of Method Article \u2014 ", { bold: true }),
        txt("\u201CVirtual screening with agentic AI: human expertise injection is still essential\u201D", { bold: true, italics: true }),
      ], { after: 200 }),

      // Para 1: Context and timeliness
      para([
        txt("I am pleased to submit the above manuscript for consideration as a Method Article in "),
        txt("Briefings in Bioinformatics", { italics: true }),
        txt(". The deployment of agentic AI\u2014large language model (LLM)-based coding agents that autonomously write and execute code\u2014in computational drug discovery has accelerated dramatically in the past year. Tools such as Claude Code, ChemCrow, and Coscientist can now generate complete virtual screening (VS) pipelines from natural-language instructions alone, promising to democratise structure-based drug discovery for researchers without deep computational expertise. However, the critical question of whether these autonomous agents produce "),
        txt("scientifically rigorous", { italics: true }),
        txt(" pipelines\u2014or merely syntactically valid but methodologically na\u00efve ones\u2014has not been systematically addressed."),
      ]),

      // Para 2: What we show
      para([
        txt("This study provides the first controlled, multi-target benchmark demonstrating that "),
        txt("human expert knowledge remains a prerequisite for high-quality agentic VS", { bold: true }),
        txt(". I show that augmenting an LLM coding agent with a concise structured skill file (212 lines of expert-curated best practices) yields two measurable improvements across 11 diverse drug targets spanning six protein classes: (i) dramatically better ligand preparation quality (86\u201399% success vs 46\u201395% for the unguided agent), and (ii) statistically significant gains in docking discrimination for targets with deep, enclosed binding sites (COX-2: \u0394AUC = +0.129, "),
        txt("p", { italics: true }),
        txt(" = 0.016; AChE: \u0394AUC = +0.157, "),
        txt("p", { italics: true }),
        txt(" = 0.004). Crucially, without the skill file, the agent loses up to 54% of its screening library to malformed ligand files\u2014a silent failure mode that would be invisible to a non-expert user. No target showed significant deterioration under skill guidance."),
      ]),

      // Para 3: Significance and fit
      para([
        txt("These findings carry a timely message for the bioinformatics community: while agentic AI lowers the barrier to entry for computational drug discovery, it does not eliminate the need for domain expertise. Rather, that expertise must be "),
        txt("encoded and injected", { italics: true }),
        txt(" into the agent\u2019s workflow\u2014in our case, as a lightweight, interpretable, version-controllable plain-text skill file. This mechanism is model-agnostic, requires no fine-tuning or infrastructure, and can be authored and peer-reviewed by domain experts. The manuscript introduces a novel intersection-based DeLong analysis that eliminates evaluation confounds arising from differential preparation failure rates\u2014a methodological contribution relevant to any benchmark comparing pipelines with heterogeneous preprocessing."),
      ]),

      // Para 4: Fit to journal
      para([
        txt("I believe this work is well suited for "),
        txt("Briefings in Bioinformatics", { italics: true }),
        txt(" because it addresses a rapidly emerging challenge at the intersection of AI and bioinformatics: how to ensure that autonomous AI agents produce reliable scientific outputs. The study combines rigorous statistical methodology (DerSimonian\u2013Laird random-effects meta-analysis, bootstrap confidence intervals, Benjamini\u2013Hochberg correction) with a practical drug discovery application across a diverse target panel."),
      ]),

      // Para 5: Declarations
      para([
        txt("The manuscript reports original research that has not been published elsewhere and is not under consideration at another journal. In accordance with ISCB guidelines, I disclose that the LLM coding agent (Claude Code, Claude Opus 4.6) was used to generate the computational pipeline code; all scientific design, analysis, interpretation, and manuscript writing were performed by the author. All code and data are publicly available under the MIT licence."),
      ]),

      // Closing
      para([txt("Thank you for considering this submission. I look forward to your response.")], { after: 300 }),

      para([txt("Yours sincerely,")], { after: 400 }),

      para([txt("Osman A.B.S.M. Gani, PhD", { bold: true })]),
      para([txt("Section for Pharmaceutical Chemistry")], { after: 40 }),
      para([txt("Department of Pharmacy")], { after: 40 }),
      para([txt("University of Oslo")], { after: 40 }),
      para([txt("0316 Oslo, Norway")], { after: 40 }),
      para([txt("osman.gani@farmasi.uio.no")]),
    ],
  }],
});

Packer.toBuffer(doc).then(buffer => {
  fs.writeFileSync("Cover_Letter_BIB.docx", buffer);
  console.log("Written: Cover_Letter_BIB.docx (" + buffer.length + " bytes)");
});
