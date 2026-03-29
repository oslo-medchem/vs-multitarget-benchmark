#!/usr/bin/env node
/**
 * Rebuild manuscript from extracted paragraphs with figures embedded
 * and Tables 2/3 removed.
 */
const fs = require("fs");
const path = require("path");
const docx = require(path.join(__dirname, "node_modules", "docx"));
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
        ImageRun, PageBreak, WidthType, AlignmentType, BorderStyle,
        ShadingType } = docx;

const FONT = "Times New Roman";
const SZ = 24;
const SZ_SM = 20;
const LINE_SPACING = 480;

// Load extracted paragraphs
const paragraphs = JSON.parse(fs.readFileSync("/tmp/manuscript_paragraphs.json", "utf8"));
const table1Data = JSON.parse(fs.readFileSync("/tmp/manuscript_table1.json", "utf8"));

// Load figures
const figFiles = {
    "Fig. 1.": { file: "figures/forest_plot.png", width: 480, height: 437 },
    "Fig. 2.": { file: "figures/fig_tables_combined.png", width: 580, height: 225 },
    "Fig. 3.": { file: "figures/fig_success_combined.png", width: 500, height: 445 },
};

// Skip these paragraphs (Table 2/3 captions and content absorbed into Fig 2)
const SKIP_PATTERNS = [
    "Table 2.",
    "Table 3.",
];

function shouldSkip(text) {
    for (const pat of SKIP_PATTERNS) {
        if (text.trim().startsWith(pat)) return true;
    }
    return false;
}

function makePara(text, opts = {}) {
    const runs = [];
    if (opts.boldAll) {
        runs.push(new TextRun({ text, font: FONT, size: SZ, bold: true }));
    } else if (opts.italicAll) {
        runs.push(new TextRun({ text, font: FONT, size: SZ, italics: true }));
    } else {
        runs.push(new TextRun({ text, font: FONT, size: SZ }));
    }
    return new Paragraph({
        alignment: opts.align || AlignmentType.JUSTIFIED,
        spacing: { after: opts.after !== undefined ? opts.after : 200, line: LINE_SPACING },
        children: runs,
    });
}

function makeRichPara(rawRuns, opts = {}) {
    const runs = rawRuns.map(([text, bold, italic]) => {
        return new TextRun({
            text,
            font: FONT,
            size: SZ,
            bold: bold || false,
            italics: italic || false,
        });
    });
    return new Paragraph({
        alignment: opts.align || AlignmentType.JUSTIFIED,
        spacing: { after: opts.after !== undefined ? opts.after : 200, line: LINE_SPACING },
        children: runs,
    });
}

function makeTable1() {
    const border = { style: BorderStyle.SINGLE, size: 1, color: "000000" };
    const borders = { top: border, bottom: border, left: border, right: border };
    const colWidths = [3000, 3013, 3013];
    
    const rows = table1Data.map((row, idx) => {
        return new TableRow({
            children: row.map((cell, ci) => {
                const isHeader = idx === 0;
                return new TableCell({
                    borders,
                    width: { size: colWidths[ci], type: WidthType.DXA },
                    shading: isHeader ? { fill: "D9E2F3", type: ShadingType.CLEAR } : undefined,
                    margins: { top: 40, bottom: 40, left: 80, right: 80 },
                    children: [new Paragraph({
                        children: [new TextRun({
                            text: cell,
                            font: FONT,
                            size: SZ_SM,
                            bold: isHeader,
                        })],
                        spacing: { line: 240 },
                    })],
                });
            }),
        });
    });
    
    return new Table({
        width: { size: 9026, type: WidthType.DXA },
        columnWidths: colWidths,
        rows,
    });
}

function makeFigure(figKey) {
    const fig = figFiles[figKey];
    if (!fig || !fs.existsSync(fig.file)) return [];
    
    const imgBuf = fs.readFileSync(fig.file);
    const elements = [];
    
    elements.push(new Paragraph({
        alignment: AlignmentType.CENTER,
        spacing: { before: 200, after: 100, line: LINE_SPACING },
        children: [new ImageRun({
            type: "png",
            data: imgBuf,
            transformation: { width: fig.width, height: fig.height },
            altText: { title: figKey, description: figKey, name: figKey },
        })],
    }));
    
    return elements;
}

// Build document
const children = [];

for (const p of paragraphs) {
    const text = p.text.trim();
    
    // Skip Table 2/3 captions
    if (shouldSkip(text)) continue;
    
    // Insert figure BEFORE its caption
    for (const figKey of Object.keys(figFiles)) {
        if (text.startsWith(figKey) || (p.bold && text.includes(figKey))) {
            const figElements = makeFigure(figKey);
            children.push(...figElements);
            break;
        }
    }
    
    // Insert Table 1 after its caption
    if (text.startsWith("Table 1")) {
        children.push(makeRichPara(p.raw_runs));
        children.push(makeTable1());
        continue;
    }
    
    // Regular paragraph - preserve formatting from runs
    children.push(makeRichPara(p.raw_runs));
}

const doc = new Document({
    sections: [{
        properties: {
            page: {
                size: { width: 11906, height: 16838 },
                margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
            },
        },
        children,
    }],
});

Packer.toBuffer(doc).then(buffer => {
    fs.writeFileSync("Manuscript_BIB_v1.docx", buffer);
    console.log(`Written: Manuscript_BIB_v1.docx (${(buffer.length/1024).toFixed(0)} KB)`);
});
