\NeedsTeXFormat{LaTeX2e}[1994/06/01]
\ProvidesClass{slides}
[2015/3/26, Mark Andrews's Beamer slides class.]

\LoadClass[11pt,xcolor=dvipsnames,professionalfonts]{beamer}

\usefonttheme{serif}

\RequirePackage{pifont,manfnt}
\RequirePackage{booktabs}
\RequirePackage[T1]{fontenc}
\RequirePackage{mathpazo}
\RequirePackage{eulervm}
\RequirePackage{tikz}
\usetikzlibrary{mindmap,trees,decorations.pathreplacing,positioning, arrows.meta}


\linespread{1.05}
\RequirePackage{xspace}


\usepackage[]{hyperref}
\hypersetup{
	pdfauthor={AUTHORS},
	pdftitle={TITLE},
	pdfkeywords={KEYWORDS},
	colorlinks,breaklinks,
	linkcolor={black},
	urlcolor={blue},
	anchorcolor={blue},
	citecolor={blue}
	}


\urlstyle{same}

\RequirePackage{rotating}
\RequirePackage{multirow}
\RequirePackage{multicol}

\usepackage{colortbl}


\usepackage{subcaption}
\captionsetup{labelformat=empty}

\newcommand\iidsim{\mathrel{\overset{\makebox[0pt]{\mbox{\normalfont\tiny iid}}}{\sim}}}
\newcommand\defeq{\mathrel{\overset{\makebox[0pt]{\mbox{\normalfont\tiny def}}}{=}}}
\newcommand{\hpd}{\textsc{hpd}\xspace}

\setbeamerfont{title}{family=\it}
\setbeamerfont{frametitle}{family=\it}
\usetikzlibrary{trees}

\usepackage{tikz-qtree-compat}
\usetikzlibrary{arrows,decorations.pathmorphing,backgrounds,positioning,fit,petri}
\tikzset{every tree node/.style={align=center, anchor=north}}
\tikzset{level 1+/.style={sibling distance=1.5\baselineskip}}

\usetikzlibrary{decorations.pathmorphing}
\usetikzlibrary{positioning}
\usetikzlibrary{fit,shapes.geometric}
\usetikzlibrary{shadings}

\newcommand{\tikzmark}[2]{%
	\tikz[remember picture,baseline=-3pt]
	\node[circle,black,draw,text=black,anchor=center,inner sep=1pt] (#1) {$#2$};}


\RequirePackage{amssymb,latexsym,amsmath,amsfonts,amscd}

\usepackage[absolute,overlay]{textpos}
\usepackage{mathtools}
\usepackage{changepage}
\usepackage{listings}


\usepackage{xspace,graphicx,nth}
\usepackage{fontawesome}
\usepackage{comment}
\usepackage{forest}
\usepackage{flexisym}
\usetikzlibrary{calc}

\usepackage{kantlipsum,setspace}
\usepackage{caption}
\captionsetup[table]{font={stretch=1.3}} 
\newcommand{\updel}{\mathop{}\!\upDelta}

\usecolortheme[named=Gray]{structure} 
\setbeamercolor{titlelike}{fg=black!60!red}
\definecolor{Mygrey}{gray}{0.55}

\let\Tiny=\tiny


\DeclareSymbolFont{legacymaths}{OT1}{cmr}{m}{n}
\DeclareMathAccent{\dot}     {\mathalpha}{legacymaths}{95}
\DeclareMathAccent{\bar}     {\mathalpha}{legacymaths}{22}
\DeclareMathAccent{\tilde}     {\mathalpha}{legacymaths}{126}


\setbeamercovered{transparent}
\mode<presentation>
\beamertemplatenavigationsymbolsempty


% Footnote with astrix
\renewcommand*{\thefootnote}{\fnsymbol{footnote}}


% Poss cite, e.g., Roeser et al.'s (2015)
\newcommand\posscite[1]{\citeauthor{#1}'s (\citeyear{#1})}

% Change figure name
\renewcommand{\figurename}{}

% Background box
\makeatletter
\newbox\@backgroundblock
\newenvironment{backgroundblock}[2]{%
	\global\setbox\@backgroundblock=\vbox\bgroup%
	\unvbox\@backgroundblock%
	\vbox to0pt\bgroup\vskip#2\hbox to0pt\bgroup\hskip#1\relax%
}{\egroup\egroup\egroup}
\addtobeamertemplate{background}{\box\@backgroundblock}{}
\makeatother

% Footnote without index
\newcommand\blfootnote[1]{%
	\begingroup
	\renewcommand\thefootnote{}\footnote{#1}%
	\addtocounter{footnote}{-1}%
	\endgroup
}

% Footnote with astrix
\renewcommand*{\thefootnote}{\fnsymbol{footnote}}

%\setbeamertemplate{itemize items}[ball]

% Bright gray
\definecolor{brightgray}{gray}{0.8}

\usepackage{linguex}

% Draw a cloud
\newcommand{\asymcloud}[2][.1]{%
	\begin{scope}[#2]
		\pgftransformscale{#1}%    
		\pgfpathmoveto{\pgfpoint{261 pt}{115 pt}} 
		\pgfpathcurveto{\pgfqpoint{70 pt}{107 pt}}
		{\pgfqpoint{137 pt}{291 pt}}
		{\pgfqpoint{260 pt}{273 pt}} 
		\pgfpathcurveto{\pgfqpoint{78 pt}{382 pt}}
		{\pgfqpoint{381 pt}{445 pt}}
		{\pgfqpoint{412 pt}{410 pt}}
		\pgfpathcurveto{\pgfqpoint{577 pt}{587 pt}}
		{\pgfqpoint{698 pt}{488 pt}}
		{\pgfqpoint{685 pt}{366 pt}}
		\pgfpathcurveto{\pgfqpoint{840 pt}{192 pt}}
		{\pgfqpoint{610 pt}{157 pt}}
		{\pgfqpoint{610 pt}{157 pt}}
		\pgfpathcurveto{\pgfqpoint{531 pt}{39 pt}}
		{\pgfqpoint{298 pt}{51 pt}}
		{\pgfqpoint{261 pt}{115 pt}}
		\pgfusepath{fill,stroke}         
	\end{scope}}  
	
	
 \tikzset{
 	invisible/.style={opacity=0},
 	visible on/.style={alt={#1{}{invisible}}},
 	alt/.code args={<#1>#2#3}{%
 		\alt<#1>{\pgfkeysalso{#2}}{\pgfkeysalso{#3}} % \pgfkeysalso doesn't change the path
 	},
 }
 