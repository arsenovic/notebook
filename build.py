#!/usr/bin/env python3
"""
Convert markdown notebooks to executable HTML with optional code execution.

What it does:
    - Executes code cells (configurable timeout)
    - Injects static/style.css for custom styling
    - Injects static/helpers.js for interactive buttons (skipped for index.md)
    - Copies img/ folder to output directory
    - Integrated CI using gh-pages. see  .github/workflows/convert-and-deploy.yml
    
Usage:
    python build.py              # Convert all .md files in ./content
    python build.py filename.md  # Convert specific file
    
Outputs HTML to ./html directory.
"""
import os
import sys
import shutil
from pathlib import Path
from typing import Optional, Dict, Any

import jupytext
import nbformat
from nbconvert import HTMLExporter
from nbconvert.preprocessors import ExecutePreprocessor

INPUT_DIR = "./content"  # Source directory for markdown files
OUTPUT_DIR = "./html"  # Destination directory for HTML output
FILE_PATTERN = "*.md"  # File glob pattern to match
EXECUTE_NOTEBOOKS = True  # Run code cells during conversion
EXECUTION_TIMEOUT = 300  # Seconds allowed per cell execution
EXCLUDE_HELPERS_JS = ["index.md"]# Files to exclude from adding JavaScript helper buttons
NBCONVERT_KWARGS = {
    "template_name": "lab",
    #"theme": "dark"
    #"exclude_input": False,
    #"exclude_output": False,
}
JUPYTEXT_KWARGS = {}

def convert_md_to_html(filename: Optional[str] = None) -> None:
    input_path = Path(INPUT_DIR).resolve()
    output_path = Path(OUTPUT_DIR).resolve()

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_path}")

    output_path.mkdir(parents=True, exist_ok=True)
    
    # If filename is provided, only convert that file; otherwise convert all
    if filename:
        md_files = [input_path / filename]
        if not md_files[0].exists():
            raise FileNotFoundError(f"File not found: {md_files[0]}")
    else:
        md_files = sorted(input_path.glob(FILE_PATTERN))
    
    exporter = HTMLExporter(**NBCONVERT_KWARGS)
    
    if EXECUTE_NOTEBOOKS:
        executor = ExecutePreprocessor(timeout=EXECUTION_TIMEOUT, kernel_name="python3")

    # Extract nbconvert styles from the first conversion
    nbconvert_styles = None

    for md_file in md_files:
        try:
            print(f"{md_file.name}: converting, ", end="")
            notebook = jupytext.read(str(md_file), **JUPYTEXT_KWARGS)

            if EXECUTE_NOTEBOOKS:
                print(f"executing, ", end="")
                try:
                    notebook, _ = executor.preprocess(notebook, {"metadata": {"path": str(input_path)}})
                except Exception as e:
                    print(f"  Warning: Notebook execution failed: {e}", file=sys.stderr)

            (body, resources) = exporter.from_notebook_node(notebook)

            # Add reference to external style.css in the head for highest precedence
            html_content = body.replace("</head>", '<link rel="stylesheet" href="static/style.css">\n</head>')

            # Add reference to external helpers script (skip for files in EXCLUDE_HELPERS_JS)
            if md_file.name not in EXCLUDE_HELPERS_JS:
                script_ref = '<script src="static/helpers.js"></script>'
                html_content = html_content.replace("</html>", script_ref + "\n</html>")

            html_file = output_path / md_file.with_suffix(".html").name
            html_file.write_text(html_content)

            print(f"done. {html_file}")

        except Exception as e:
            print(f"Error converting {md_file.name}: {e}", file=sys.stderr)


    # Copy img folder to html output directory
    img_source = input_path / "img"
    img_dest = output_path / "img"
    if img_source.exists():
        shutil.copytree(img_source, img_dest, dirs_exist_ok=True)
        print(f"Copied images to {img_dest}")


def main():
    # Get optional filename argument
    filename = sys.argv[1] if len(sys.argv) > 1 else None
    convert_md_to_html(filename)


if __name__ == "__main__":
    main()
