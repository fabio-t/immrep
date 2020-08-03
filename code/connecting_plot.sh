#!/bin/bash

usage() { echo "Usage: $0 <n1> <n2> <n3> <n4> <n5> [--n12 <mean>,<stdev>] [--n13 <mean>,<stdev>] [--n14 <mean>,<stdev>] [--n15 <mean>,<stdev>] [--n23 <mean>,<stdev>] [--n24 <mean>,<stdev>] [--n25 <mean>,<stdev>] [--n34 <mean>,<stdev>] [--n35 <mean>,<stdev>] [--n45 <mean>,<stdev>]" 1>&2; exit 1; }

ARGUMENT_LIST=(
    "n12"
    "n13"
    "n14"
    "n15"
    "n23"
    "n24"
    "n25"
    "n34"
    "n35"
    "n45"
)

# read arguments
opts=$(getopt \
    --longoptions "$(printf "%s:," "${ARGUMENT_LIST[@]}")" \
    --name "$(basename "$0")" \
    --options "" \
    -- "$@"
)

M11=$1
M13=$2
M31=$3
M33=$4
M52=$5
shift 5

if [ "$M11" == "." ]
then
    unset M11
fi
if [ "$M13" == "." ]
then
    unset M13
fi
if [ "$M31" == "." ]
then
    unset M31
fi
if [ "$M33" == "." ]
then
    unset M33
fi
if [ "$M52" == "." ]
then
    unset M52
fi

eval set --$opts

while [[ $# -gt 0 ]]; do
    case "$1" in
        --n12)
            set -f # disable glob
            IFS="," read -r M11_M13_mean M11_M13_stdev <<< "$2"
            shift 2
            ;;

        --n13)
            set -f # disable glob
            IFS="," read -r M11_M31_mean M11_M31_stdev <<< "$2"
            shift 2
            ;;

        --n14)
            set -f # disable glob
            IFS="," read -r M11_M33_mean M11_M33_stdev <<< "$2"
            shift 2
            ;;

        --n15)
            set -f # disable glob
            IFS="," read -r M11_M52_mean M11_M52_stdev <<< "$2"
            shift 2
            ;;

        --n23)
            set -f # disable glob
            IFS="," read -r M13_M31_mean M13_M31_stdev <<< "$2"
            shift 2
            ;;

        --n24)
            set -f # disable glob
            IFS="," read -r M13_M33_mean M13_M33_stdev <<< "$2"
            shift 2
            ;;

        --n25)
            set -f # disable glob
            IFS="," read -r M13_M52_mean M13_M52_stdev <<< "$2"
            shift 2
            ;;

        --n34)
            set -f # disable glob
            IFS="," read -r M31_M33_mean M31_M33_stdev <<< "$2"
            shift 2
            ;;

        --n35)
            set -f # disable glob
            IFS="," read -r M31_M52_mean M31_M52_stdev <<< "$2"
            shift 2
            ;;

        --n45)
            set -f # disable glob
            IFS="," read -r M33_M52_mean M33_M52_stdev <<< "$2"
            shift 2
            ;;

        *)
            break
            ;;
    esac
done

# Create a temporary directory
curdir=$( pwd )
tmpdir=$( mktemp -dt "latex.XXXXXXXX" )

# Set up a trap to clean up and return back when the script ends
# for some reason
clean_up () {
    cd "$curdir"
    [ -d "$tmpdir" ] && rm -rf "$tmpdir"
    exit
}
trap 'clean_up' EXIT SIGHUP SIGINT SIGQUIT SIGTERM

_f() {
    tmp=$( bc -l <<<"5*$tmp" )
}

if [ -n "$M11" ] && [ -n "$M13" ]
then
    tmp=$M11_M13_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M11_M13="(m-1-1) edge [${line},line width=${tmp}pt,-] node [above] {\\tiny \$$M11_M13_mean \\pm $M11_M13_stdev\$} (m-1-3)"
fi
if [ -n "$M11" ] && [ -n "$M31" ]
then
    tmp=$M11_M31_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M11_M31="(m-1-1) edge [${line},line width=${tmp}pt,-] node [left] {\\tiny \$$M11_M31_mean \\pm $M11_M31_stdev\$} (m-3-1)"
fi
if [ -n "$M11" ] && [ -n "$M33" ]
then
    tmp=$M11_M33_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M11_M33="(m-1-1) edge [${line},line width=${tmp}pt,-] node [below right=4pt] {\\tiny \$$M11_M33_mean \\pm $M11_M33_stdev\$} (m-3-3)"
fi
if [ -n "$M11" ] && [ -n "$M52" ]
then
    tmp=$M11_M52_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M11_M52="(m-1-1) edge [${line},line width=${tmp}pt,-] node [below=11pt] {\\tiny \$$M11_M52_mean \\pm $M11_M52_stdev\$} (m-5-2)"
fi
if [ -n "$M13" ] && [ -n "$M31" ]
then
    tmp=$M13_M31_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M13_M31="(m-1-3) edge [${line},line width=${tmp}pt,-] node [below left=4pt] {\\tiny \$$M13_M31_mean \\pm $M13_M31_stdev\$} (m-3-1)"
fi
if [ -n "$M13" ] && [ -n "$M33" ]
then
    tmp=$M13_M33_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M13_M33="(m-1-3) edge [${line},line width=${tmp}pt,-] node [right] {\\tiny \$$M13_M33_mean \\pm $M13_M33_stdev\$} (m-3-3)"
fi
if [ -n "$M13" ] && [ -n "$M52" ]
then
    tmp=$M13_M52_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M13_M52="(m-1-3) edge [${line},line width=${tmp}pt,-] node [below=11pt] {\\tiny \$$M13_M52_mean \\pm $M13_M52_stdev\$} (m-5-2)"
fi
if [ -n "$M31" ] && [ -n "$M33" ]
then
    tmp=$M31_M33_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M31_M33="(m-3-1) edge [${line},line width=${tmp}pt,-] node [below] {\\tiny \$$M31_M33_mean \\pm $M31_M33_stdev\$} (m-3-3)"
fi
if [ -n "$M31" ] && [ -n "$M52" ]
then
    tmp=$M31_M52_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M31_M52="(m-3-1) edge [${line},line width=${tmp}pt,-] node [left=6pt] {\\tiny \$$M31_M52_mean \\pm $M31_M52_stdev\$} (m-5-2)"
fi
if [ -n "$M33" ] && [ -n "$M52" ]
then
    tmp=$M33_M52_mean
    line="solid"
    if [ 1 -eq "$(echo "${tmp} < 0.4" | bc)" ]
    then
        line="dashed"
    fi
    _f
    M33_M52="(m-3-3) edge [${line},line width=${tmp}pt,-] node [right=6pt] {\\tiny \$$M33_M52_mean \\pm $M33_M52_stdev\$} (m-5-2)"
fi

# Switch to the temp. directory and extract the .tex file
cd $tmpdir
# Quoting the 'THEEND' string prevents $-expansion.
cat > tmp.tex <<EOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Welcome to Overleaf --- just edit your LaTeX on the left,
% and we'll compile it for you on the right. If you open the
% 'Share' menu, you can invite other users to edit at the same
% time. See www.overleaf.com/learn for more info. Enjoy!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A simple commutative diagram
% Stefan Kottwitz
\documentclass{article}
\usepackage{tikz}
%%%<
\usepackage{verbatim}
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{tikzpicture}
\setlength\PreviewBorder{5pt}%
%%%>
\usetikzlibrary{matrix}
\begin{document}
\tikzset{every node/.style={fill=white}}
\begin{tikzpicture}
  \matrix (m) [matrix of math nodes,row sep=3em,column sep=4em,minimum width=2em]
  {
     $M11 & & $M13 \\\\
     & & \\\\
     $M31 & & $M33 \\\\
     & & \\\\
     & $M52 & \\\\
  };
  \path[-stealth]
    $M11_M13
    $M11_M31
    $M11_M33
    $M11_M52
    $M13_M31
    $M13_M33
    $M13_M52
    $M31_M33
    $M31_M52
    $M33_M52
    ;
\end{tikzpicture}
\end{document}
EOF

cat tmp.tex

# If the file extracts succesfully, try to run pdflatex 3 times.
# If something fails, print a warning and exit
if [[ -f 'tmp.tex' ]]
then
   for i in {1..3}
   do
      if pdflatex tmp.tex
      then
         echo "Pdflatex run $i finished."
      else
         echo "Pdflatex run $i failed."
         exit 2
      fi
   done
else
   echo "Error extracting .tex file"
   exit 1
fi

# Copy the resulting .pdf file to original directory and exit
cp tmp.pdf $curdir
exit 0
