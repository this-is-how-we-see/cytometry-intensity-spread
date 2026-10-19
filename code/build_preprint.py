"""Build manuscript/preprint.pdf from manuscript/preprint.md: markdown to HTML with the figures
placed above their legends, then headless Chrome prints the HTML to PDF.
Usage: python3 build_preprint.py <project_dir>"""
import re, subprocess, sys
from pathlib import Path
import markdown

CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
FIGS = {1: "fig1_scaling_cost", 2: "fig2_where_the_loss_falls", 3: "fig3_representations", 4: "fig4_between_site", 5: "fig5_hipc_transform", 6: "fig6_cofactor"}
CSS = """
@page { size: Letter; margin: 2.2cm 2.2cm 2.4cm; @bottom-center { content: counter(page); } }
body { font-family: "Charter", "Georgia", serif; font-size: 10.5pt; line-height: 1.45; color: #111; }
h1 { font-size: 17pt; line-height: 1.25; margin: 0 0 .6em; }
h2 { font-size: 12.5pt; margin: 1.6em 0 .5em; }
h3 { font-size: 11pt; margin: 1.2em 0 .4em; font-style: italic; font-weight: normal; }
p { margin: 0 0 .7em; text-align: left; hyphens: auto; }
a { color: #1a4f8b; text-decoration: none; }
table { border-collapse: collapse; width: 100%; font-family: "Helvetica Neue", Arial, sans-serif; font-size: 7.6pt; margin: .4em 0 1.2em; page-break-inside: auto; }
th, td { border-bottom: .5pt solid #bbb; padding: 2.5pt 4pt; text-align: left; vertical-align: top; }
th { border-bottom: .8pt solid #555; font-weight: 600; }
tr { page-break-inside: avoid; }
figure { margin: 1em 0 1.4em; page-break-inside: avoid; }
figure img { width: 100%; }
figure.narrow img { width: 72%; }
figcaption { font-size: 9pt; line-height: 1.35; margin-top: .4em; }
ol li { margin-bottom: .25em; font-size: 9pt; }
"""


def main(root):
    root = Path(root); m = root / "manuscript"
    md = (m / "preprint.md").read_text()
    def fig(match):
        n = int(match.group(1)); body = match.group(0)
        cls = ' class="narrow"' if n == 2 else ""
        cap = markdown.markdown(body).removeprefix("<p>").removesuffix("</p>")
        return f'<figure{cls}><img src="figures/{FIGS[n]}.png" alt="Figure {n}"><figcaption>{cap}</figcaption></figure>'
    md = re.sub(r"^\*\*Figure (\d)\..*$", fig, md, flags=re.M)
    html = markdown.markdown(md, extensions=["tables", "sane_lists"])
    out = m / "preprint.html"
    out.write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>Intensity spread is information</title><style>{CSS}</style></head><body>{html}</body></html>")
    pdf = m / "preprint.pdf"
    subprocess.run([CHROME, "--headless=new", "--disable-gpu", "--no-pdf-header-footer", f"--print-to-pdf={pdf}", out.resolve().as_uri()],
                   check=True, capture_output=True, timeout=180)
    print(pdf, pdf.stat().st_size, "bytes")


if __name__ == "__main__":
    main(sys.argv[1])
