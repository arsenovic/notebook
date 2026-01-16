#!/usr/bin/env python3

import os
import sys
import shutil
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
    "exclude_input": False,
    "exclude_output": False,
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

    # Extract nbconvert styles from the first conversion
    nbconvert_styles = None

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

            # Extract nbconvert styles on first run
            if nbconvert_styles is None:
                import re
                style_match = re.search(r'<style[^>]*>(.*?)</style>', body, re.DOTALL)
                if style_match:
                    nbconvert_styles = style_match.group(1)
            
            # Remove embedded styles and add external CSS references
            import re
            html_content = re.sub(r'<style[^>]*>.*?</style>', '', body, flags=re.DOTALL)
            html_content = html_content.replace("</head>", 
                '<link rel="stylesheet" href="nbconvert.css">\n' +
                '<link rel="stylesheet" href="style.css">\n</head>')

            # Add reference to external helpers script
            script_ref = '<script src="helpers.js"></script>'
            html_content = html_content.replace("</html>", script_ref + "\n</html>")

            html_file = output_path / md_file.with_suffix(".html").name
            html_file.write_text(html_content)

            print(f" {html_file}")

        except Exception as e:
            print(f"Error converting {md_file.name}: {e}", file=sys.stderr)
    
    # Save extracted nbconvert styles to nbconvert.css
    if nbconvert_styles:
        nbconvert_css_file = output_path / "nbconvert.css"
        nbconvert_css_file.write_text(nbconvert_styles)
        print(f"Saved nbconvert styles to {nbconvert_css_file}")

    # Copy img folder to html output directory
    img_source = input_path / "img"
    img_dest = output_path / "img"
    if img_source.exists():
        shutil.copytree(img_source, img_dest, dirs_exist_ok=True)
        print(f"Copied images to {img_dest}")


def main():
    convert_md_to_html()


if __name__ == "__main__":
    main()
