import { mkdtempSync, readFileSync, rmSync, writeFileSync, existsSync } from "node:fs";
import { tmpdir } from "node:os";
import { basename, join } from "node:path";
import { pathToFileURL } from "node:url";
import { spawnSync } from "node:child_process";

const root = String.raw`C:\Users\ky\OneDrive\Source Shared\DyGist\Data\Excitations\AnlzRoutine`;
const outputFolder = "[06.21].98.FullTime.PrflModlEvol.Tailess.sansMask";
const edgePath = String.raw`C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe`;

const sources = [
    { id: 96, label: "96 sansMask kept tail", folder: "[06.21].96.FullTime.PrflModlEvol.KeptTail.sansMask" },
    { id: 95, label: "95 avecMask kept tail", folder: "[06.21].95.FullTime.PrflModlEvol.KeptTail.avecMask" },
    { id: 98, label: "98 sansMask tailess", folder: outputFolder },
    { id: 97, label: "97 avecMask tailess", folder: "[06.21].97.FullTime.PrflModlEvol.Tailess.avecMask" },
];

const groups = {
    CFNM: [
        [5.311, 95],
        [5.313, 82],
        [5.316, 52],
        [5.318, 80],
        [5.322, 67],
        [5.325, 96],
        [5.326, 68],
        [5.328, 50],
        [5.332, 81],
        [5.333, 51],
        [5.336, 79],
        [5.338, 53],
    ],
    NTRC: [
        [5.314, 29],
        [5.316, 28],
        [5.318, 27],
        [5.322, 26],
        [5.326, 25],
        [5.332, 61],
        [5.336, 62],
        [5.340, 63],
        [5.343, 64],
    ],
};

const pageWidthPx = 1600;
const stackGapPx = 0;
const rowGapPx = 10;
const marginPx = 10;

function formatIb(ib) {
    return ib.toFixed(3);
}

function svgName(tag, ib, runid) {
    return `prfl_modl_evol_[${tag}_${formatIb(ib)}_r${runid}].svg`;
}

function svgDimensions(svgText, path) {
    const svgTag = svgText.match(/<svg\b[^>]*>/i)?.[0];
    if (!svgTag) {
        throw new Error(`No <svg> tag found in ${path}`);
    }
    const width = svgTag.match(/\bwidth="([\d.]+)"/i)?.[1];
    const height = svgTag.match(/\bheight="([\d.]+)"/i)?.[1];
    if (!width || !height) {
        throw new Error(`No numeric width/height found in ${path}`);
    }
    return { width: Number(width), height: Number(height) };
}

function readSvgAsDataUrl(path) {
    const svg = readFileSync(path, "utf8");
    const dimensions = svgDimensions(svg, path);
    const dataUrl = `data:image/svg+xml;base64,${Buffer.from(svg, "utf8").toString("base64")}`;
    return { dataUrl, dimensions };
}

function buildHtml(tag, entries) {
    const imageWidthPx = pageWidthPx - 2 * marginPx;
    const rows = [];
    let pageHeightPx = 2 * marginPx;

    for (const [ib, runid] of entries) {
        const figureHtml = [];
        let stackHeightPx = 0;
        let expectedRatio = null;

        for (const source of sources) {
            const path = join(root, source.folder, svgName(tag, ib, runid));
            if (!existsSync(path)) {
                throw new Error(`Missing source SVG: ${path}`);
            }
            const { dataUrl, dimensions } = readSvgAsDataUrl(path);
            const ratio = dimensions.height / dimensions.width;
            expectedRatio ??= ratio;
            if (Math.abs(ratio - expectedRatio) > 0.001) {
                throw new Error(`Unexpected aspect ratio in ${path}`);
            }
            const heightPx = imageWidthPx * ratio;
            stackHeightPx += heightPx + stackGapPx;
            figureHtml.push(`
                <div class="figure">
                    <img src="${dataUrl}" width="${imageWidthPx}" height="${heightPx}" alt="${basename(path)}">
                </div>`);
        }

        stackHeightPx -= stackGapPx;
        const rowHeightPx = Math.ceil(stackHeightPx);
        pageHeightPx += rowHeightPx + rowGapPx;
        rows.push(`
            <section class="row" style="height: ${rowHeightPx}px">
                <div class="stack">
                    ${figureHtml.join("\n")}
                </div>
            </section>`);
    }

    pageHeightPx -= rowGapPx;

    return `<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>${tag} profile modulation comparison</title>
<style>
@page {
    size: ${pageWidthPx}px ${Math.ceil(pageHeightPx)}px;
    margin: 0;
}
* {
    box-sizing: border-box;
}
html,
body {
    width: ${pageWidthPx}px;
    min-height: ${Math.ceil(pageHeightPx)}px;
    margin: 0;
    background: white;
    color: #111;
    font-family: Arial, Helvetica, sans-serif;
}
body {
    padding: ${marginPx}px;
}
.row {
    margin-bottom: ${rowGapPx}px;
    break-inside: avoid;
}
.row:last-child {
    margin-bottom: 0;
}
.stack {
    width: ${imageWidthPx}px;
}
.figure {
    width: ${imageWidthPx}px;
    margin: 0 0 ${stackGapPx}px 0;
}
.figure img {
    display: block;
    width: ${imageWidthPx}px;
    height: auto;
}
</style>
</head>
<body>
${rows.join("\n")}
</body>
</html>`;
}

function printPdf(tag) {
    if (!existsSync(edgePath)) {
        throw new Error(`Microsoft Edge was not found at ${edgePath}`);
    }

    const tmpDir = mkdtempSync(join(tmpdir(), `excitation-${tag}-`));
    const htmlPath = join(tmpDir, `comparison.${tag}.html`);
    const pdfPath = join(root, outputFolder, `comparison.[${tag}].pdf`);

    try {
        writeFileSync(htmlPath, buildHtml(tag, groups[tag]), "utf8");
        const result = spawnSync(edgePath, [
            "--headless",
            "--disable-gpu",
            "--no-pdf-header-footer",
            `--print-to-pdf=${pdfPath}`,
            pathToFileURL(htmlPath).href,
        ], { encoding: "utf8" });

        if (result.status !== 0) {
            throw new Error(`Edge failed for ${tag}:\n${result.stderr || result.stdout}`);
        }
        if (!existsSync(pdfPath)) {
            throw new Error(`Edge reported success but did not create ${pdfPath}`);
        }
        console.log(`Wrote ${pdfPath}`);
    } finally {
        rmSync(tmpDir, { recursive: true, force: true });
    }
}

for (const tag of Object.keys(groups)) {
    printPdf(tag);
}
