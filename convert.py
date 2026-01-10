#!/usr/bin/env python3

import os
import sys
from pathlib import Path
from typing import Optional, Dict, Any

import jupytext
import nbformat
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor

INPUT_DIR = "./content"
OUTPUT_DIR = "./html"
FILE_PATTERN = "*.md" 
EXECUTE_NOTEBOOKS = True
EXECUTION_TIMEOUT = 300

NBCONVERT_KWARGS = {
    "template_name": "lab",
    #"theme": "dark"
}

JUPYTEXT_KWARGS = {}

def convert_md_to_html() -> None:
    input_path = Path(INPUT_DIR).resolve()
    output_path = Path(OUTPUT_DIR).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_path}")

    output_path.mkdir(parents=True, exist_ok=True)
    md_files = sorted(input_path.glob(FILE_PATTERN))
    exporter = HTMLExporter(**NBCONVERT_KWARGS)
    
    if EXECUTE_NOTEBOOKS:
        executor = ExecutePreprocessor(timeout=EXECUTION_TIMEOUT, kernel_name="python3")

    for md_file in md_files:
        try:
            print(f"Converting: {md_file.name}")
            notebook = jupytext.read(str(md_file), **JUPYTEXT_KWARGS)

            if EXECUTE_NOTEBOOKS:
                print(f"  Executing notebook...")
                try:
                    notebook, _ = executor.preprocess(notebook, {"metadata": {"path": str(input_path)}})
                except Exception as e:
                    print(f"  Warning: Notebook execution failed: {e}", file=sys.stderr)

            (body, resources) = exporter.from_notebook_node(notebook)

            # Add reference to external style.css right before </body> so it overrides nbconvert styles
            html_content = body.replace("</body>", '<link rel="stylesheet" href="style.css">\n</body>')

            # Add reference to external helpers script
            script_ref = '<script src="helpers.js"></script>'
            html_content = html_content.replace("</html>", script_ref + "\n</html>")

            html_file = output_path / md_file.with_suffix(".html").name
            html_file.write_text(html_content)

            print(f" {html_file}")

        except Exception as e:
            print(f"Error converting {md_file.name}: {e}", file=sys.stderr)


def main():
    convert_md_to_html()


if __name__ == "__main__":
    main()
