#!/bin/bash

for f in *.m
do
	echo "removing 'clear;'s from file $f"
	sed s/"clear;"// $f > "_$f" && mv -f "_$f" $f
done

echo "xlabel(\"\$x\$\"); ylabel(\"\$y\$\"); print -dpdflatexstandalone 'unicycle_x1x2_sherlock.pdf';" >> Unicycle12.m
octave Unicycle12.m
pdflatex unicycle_x1x2_sherlock
echo "xlabel(\"\$\\\theta\$\"); ylabel(\"\$v\$\"); print -dpdflatexstandalone 'unicycle_x3x4_sherlock.pdf';" >> Unicycle34.m
octave Unicycle34.m
pdflatex unicycle_x3x4_sherlock

echo "xlabel(\"\$x_1\$\"); ylabel(\"\$x_2\$\"); print -dpdflatexstandalone 'tora_x1x2_sherlock.pdf';" >> TORA12.m
octave TORA12.m
pdflatex tora_x1x2_sherlock
echo "xlabel(\"\$x_3\$\"); ylabel(\"\$x_4\$\"); print -dpdflatexstandalone 'tora_x3x4_sherlock.pdf';" >> TORA34.m
octave TORA34.m
pdflatex tora_x3x4_sherlock

echo "xlabel(\"\$x_{lead}\$\"); ylabel(\"\$x_{ego}\$\"); print -dpdflatexstandalone 'acc_sherlock.pdf';" >> ACC.m
octave ACC.m
pdflatex acc_sherlock

echo "xlabel(\"\$\\\theta\$\"); ylabel(\"\$\\\theta'\$\"); print -dpdflatexstandalone 'singlependulum_sherlock.pdf';" >> SinglePendulum.m
octave SinglePendulum.m
pdflatex singlependulum_sherlock

echo "xlabel(\"\$\\\theta\$\"); ylabel(\"\$\\\theta'\$\"); print -dpdflatexstandalone 'doublependulum_more_robust_x1x2_sherlock.pdf';" >> DoublePendulum_mr.m
octave DoublePendulum_mr.m
pdflatex doublependulum_more_robust_x1x2_sherlock

echo "xlabel(\"\$y\$\"); ylabel(\"\$\\\phi\$\"); print -dpdflatexstandalone 'airplane_yphi_sherlock.pdf';" >> Airplane_yphi.m
octave Airplane_yphi.m
pdflatex airplane_yphi_sherlock
echo "xlabel(\"\$\\\theta\$\"); ylabel(\"\$\\\psi\$\"); print -dpdflatexstandalone 'airplane_thetapsi_sherlock.pdf';" >> Airplane_thetapsi.m
octave Airplane_thetapsi.m
pdflatex airplane_thetapsi_sherlock
