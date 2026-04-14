#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LATEX_DIR="$ROOT_DIR/latex"
PDF_DIR="$ROOT_DIR/pdfs"
TEXMFVAR_DIR="$ROOT_DIR/.texmf-var"

mkdir -p "$PDF_DIR"
mkdir -p "$TEXMFVAR_DIR"

shopt -s nullglob
tex_files=("$LATEX_DIR"/*.tex)

if [ ${#tex_files[@]} -eq 0 ]; then
  echo "No LaTeX files found in $LATEX_DIR"
  exit 0
fi

for tex_file in "${tex_files[@]}"; do
  base_name="$(basename "$tex_file" .tex)"
  echo "Compiling $base_name.tex"
  TEXMFVAR="$TEXMFVAR_DIR" pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$LATEX_DIR" "$tex_file" >/dev/null
  TEXMFVAR="$TEXMFVAR_DIR" pdflatex -interaction=nonstopmode -halt-on-error -output-directory "$LATEX_DIR" "$tex_file" >/dev/null
  mv "$LATEX_DIR/$base_name.pdf" "$PDF_DIR/$base_name.pdf"
done

find "$LATEX_DIR" -maxdepth 1 -type f \
  \( -name '*.aux' -o -name '*.log' -o -name '*.out' -o -name '*.toc' -o -name '*.nav' -o -name '*.snm' -o -name '*.fls' -o -name '*.fdb_latexmk' -o -name '*.brf' \) \
  -delete

echo "PDFs written to $PDF_DIR"
